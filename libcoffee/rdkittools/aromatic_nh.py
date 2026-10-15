from rdkit import Chem, rdBase


class AromaticNHAddedWarning(UserWarning):
    """
    Warning issued when hydrogens are added to aromatic nitrogens to make a molecule kekulizable.
    """


def _unkekulizable_atoms(mol: Chem.Mol) -> set[int] | None:
    """
    Returns the atoms of a ring system that cannot be kekulized, or None. Raises ValueError for other problems.
    """
    for problem in Chem.DetectChemistryProblems(mol):
        if problem.GetType() != "KekulizeException":
            raise ValueError(f"{problem.GetType()}: {problem.Message()}")
        return set(problem.GetAtomIndices())
    return None


def _is_nh_candidate(atom: Chem.Atom) -> bool:
    """
    Returns whether the atom is a pyridine-type aromatic nitrogen (two heavy-atom neighbors, no hydrogen, no charge).
    """
    neighbors = [nb.GetAtomicNum() for nb in atom.GetNeighbors()]
    # before sanitization, molfile input marks only bonds (not atoms) as aromatic
    aromatic = atom.GetIsAromatic() or any(b.GetIsAromatic() for b in atom.GetBonds())
    return (
        atom.GetAtomicNum() == 7
        and aromatic
        and atom.GetFormalCharge() == 0
        and atom.GetNumExplicitHs() == 0
        and 1 not in neighbors
        and len(neighbors) == 2
    )


def _set_h(atom: Chem.Atom, has_h: bool) -> None:
    atom.SetNumExplicitHs(int(has_h))
    atom.SetNoImplicit(has_h)


def add_missing_aromatic_nh(mol: Chem.Mol) -> tuple[Chem.Mol, list[int]]:
    """
    Adds hydrogens to aromatic nitrogens of an unsanitized molecule (e.g. "c1ccnc1") until it can be kekulized.
    Only one feasible assignment is searched, so the result may be a tautomer of the intended molecule.
    Returns the sanitized molecule and the indices of the nitrogens hydrogens were added to.
    """
    rwmol = Chem.RWMol(mol)
    added: list[int] = []
    with rdBase.BlockLogs():  # suppress RDKit logs of the trial-and-error kekulization
        while (bad := _unkekulizable_atoms(rwmol)) is not None:
            for idx in sorted(bad):
                atom = rwmol.GetAtomWithIdx(idx)
                if not _is_nh_candidate(atom):
                    continue
                _set_h(atom, True)
                try:
                    after = _unkekulizable_atoms(rwmol)
                except ValueError:
                    after = bad
                if after is None or not after & bad:
                    added.append(idx)
                    break
                _set_h(atom, False)
            else:
                raise ValueError(f"No feasible hydrogen assignment on aromatic nitrogens for atoms {sorted(bad)}")
    Chem.SanitizeMol(rwmol)
    return rwmol.GetMol(), added
