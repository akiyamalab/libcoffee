from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt
from rdkit import Chem
from rdkit.Chem import AllChem

from libcoffee.core.molbase import MolBase


class RDKitMol(MolBase):
    """
    rdkit.Chem.Mol wrapper class
    """

    def __init__(self, mol: Chem.Mol):
        if not isinstance(mol, Chem.Mol):
            raise ValueError(
                f"The argument of RDKitMol constructor should be an instance of rdkit.Chem.Mol, not {type(mol)}"
            )
        super().__init__(mol)

    @property
    def _atoms(self) -> tuple[Chem.Atom, ...]:
        return tuple(m for m in self._mol.GetAtoms())

    @_atoms.setter
    def _atoms(self, atoms: tuple[Chem.Atom, ...]) -> None:
        if len(atoms) != len(self._atoms):
            raise ValueError("Length of atoms should be equal to the number of atoms")
        for i in range(len(atoms)):
            self._mol.ReplaceAtom(i, atoms[i])

    @property
    def bonds(self) -> tuple[Chem.Bond, ...]:
        return tuple(b for b in self._mol.GetBonds())

    @property
    def isotopes(self) -> npt.NDArray[np.int32]:
        # a.GetIsotope() should return int but typing says return Any
        return np.array([a.GetIsotope() for a in self._atoms], dtype=np.int32)

    @isotopes.setter
    def isotopes(self, isotopes: npt.NDArray[np.int32]) -> None:
        if len(isotopes) != len(self._atoms):
            raise ValueError("Length of isotopes should be equal to the number of atoms")
        for i in range(len(self._atoms)):
            self._atoms[i].SetIsotope(int(isotopes[i]))

    @property
    def name(self) -> str:
        return self._mol.GetProp("_Name")  # type: ignore[no-any-return]

    @name.setter
    def name(self, name: str) -> None:
        self._mol.SetProp("_Name", name)

    @property
    def heavy_atom_indices(self) -> npt.NDArray[np.intp]:
        atomic_nums = [a.GetAtomicNum() for a in self._atoms]
        return np.where(np.array(atomic_nums) > 1)[0]

    def get_smiles(self, kekulize: bool = False) -> str:
        return Chem.MolToSmiles(self._mol, kekuleSmiles=kekulize)

    def get_attr(self, attr_name: str) -> Any:
        return self._mol.GetProp(attr_name)

    def set_attr(self, attr_name: str, value: Any) -> None:
        self._mol.SetProp(attr_name, str(value))

    def has_attr(self, attr_name: str) -> bool:
        return self._mol.HasProp(attr_name)  # type: ignore[no-any-return]

    def get_coordinates(self, only_heavy_atom: bool = False) -> npt.NDArray[np.float64]:
        conf = self._mol.GetConformer()
        coords = np.array([conf.GetAtomPosition(i) for i in range(self._mol.GetNumAtoms())])
        if only_heavy_atom:
            coords = coords[self.heavy_atom_indices]
        return coords

    def generate_coordinates(self, temporary_add_hydrogens: bool = False) -> None:
        if temporary_add_hydrogens:
            self._mol = Chem.AddHs(self._mol)

        try:
            AllChem.EmbedMolecule(self.raw_mol, useRandomCoords=True)  # type: ignore[attr-defined]
            AllChem.UFFOptimizeMolecule(self.raw_mol, maxIters=100)  # type: ignore[attr-defined]
        except:
            pass

        params = AllChem.ETKDGv2()  # type: ignore[attr-defined]
        params.randomSeed = 1
        params.numThreads = 50
        params.pruneRmsThresh = 0.1
        params.useRandomCoords = True
        params.maxAttempts = 1000

        count = 0
        while count < 10:
            conformers = AllChem.EmbedMultipleConfs(self.raw_mol, numConfs=10, params=params)  # type: ignore[attr-defined]

            if len(conformers) == 0:
                count += 1
                params.randomSeed += 1   
            break

        if len(conformers) == 0:
            raise ValueError(f"Failed to generate coordinates. {self.name}")

        AllChem.UFFOptimizeMolecule(self.raw_mol, maxIters=100)  # type: ignore[attr-defined]

        if temporary_add_hydrogens:
            self._mol = Chem.RemoveHs(self._mol)

    def has_coordinates(self) -> bool:
        return self._mol.GetNumConformers() > 0  # type: ignore[no-any-return]

    def add_hydrogens(self) -> RDKitMol:
        self._mol = Chem.AddHs(self._mol)
        return self

    def remove_hydrogens(self) -> RDKitMol:
        self._mol = Chem.RemoveHs(self._mol)
        return self

    def extract_submol(self, atom_idxs: npt.NDArray[np.intp]) -> RDKitMol:
        rw_mol = Chem.RWMol(self.raw_mol)
        idx_remove_atoms = set(range(self.raw_mol.GetNumAtoms())) - set(atom_idxs)
        atomidxs = sorted(idx_remove_atoms)[::-1]
        for idx in atomidxs:
            rw_mol.RemoveAtom(idx)
        return RDKitMol(rw_mol.GetMol())

    def merge(self, mol: RDKitMol, aps: tuple[int, int] | None = None) -> RDKitMol:
        rwmol = Chem.RWMol(self.raw_mol)

        merged_mol = Chem.RWMol(Chem.CombineMols(rwmol, mol.raw_mol))

        if aps is not None:
            idx1, idx2 = aps
            offset = self.raw_mol.GetNumAtoms()
            merged_mol.AddBond(idx1, idx2 + offset, Chem.BondType.SINGLE)

        Chem.SanitizeMol(merged_mol)
        return RDKitMol(merged_mol.GetMol())

    @classmethod
    def reconstruct_from_fragments(cls, frags: tuple[RDKitMol, ...]) -> RDKitMol:  # type: ignore[override]
        raise NotImplementedError

    @classmethod
    def from_smiles(cls, smiles: str) -> RDKitMol:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Failed to generate a molecule from SMILES: {smiles}")
        return RDKitMol(mol)

    @classmethod
    def read_sdf(cls, file_path: Path) -> tuple[RDKitMol, ...]:
        suppl = Chem.SDMolSupplier(str(file_path))
        return tuple(RDKitMol(m) for m in suppl if m is not None)

    @classmethod
    def write_sdf(cls, file_path: Path, mols: tuple[MolBase, ...]) -> None:
        """
        Writes the given molecules to an SDF file
        """
        writer = Chem.SDWriter(str(file_path))
        for mol in mols:
            writer.write(mol.raw_mol)
        writer.close()
