from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt
from rdkit import Chem

from libcoffee.molecule.molbase import MolBase


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

    @atoms.setter
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

    def has_coordinates(self) -> bool:
        return self._mol.GetNumConformers() > 0

    def add_hydrogens(self) -> "RDKitMol":
        self._mol = Chem.AddHs(self._mol)
        return self

    def remove_hydrogens(self) -> "RDKitMol":
        self._mol = Chem.RemoveHs(self._mol)
        return self

    def extract_submol(self, atom_idxs: npt.NDArray[np.intp]) -> "RDKitMol":
        rw_mol = Chem.RWMol(self.raw_mol)
        idx_remove_atoms = set(range(self.raw_mol.GetNumAtoms())) - set(atom_idxs)
        atomidxs = sorted(idx_remove_atoms)[::-1]
        for idx in atomidxs:
            rw_mol.RemoveAtom(idx)
        return RDKitMol(rw_mol.GetMol())

    def merge(self, mol: "RDKitMol", aps: tuple[int, int] | None = None) -> "RDKitMol":
        raise NotImplementedError

    @classmethod
    def read_sdf(cls, file_path: Path) -> tuple["RDKitMol", ...]:
        suppl = Chem.SDMolSupplier(str(file_path))
        return tuple(RDKitMol(m) for m in suppl if m is not None)

    @classmethod
    def write_sdf(cls, file_path: Path, mols: tuple["MolBase", ...]) -> None:
        """
        Writes the given molecules to an SDF file
        """
        writer = Chem.SDWriter(str(file_path))
        for mol in mols:
            writer.write(mol.raw_mol)
        writer.close()
