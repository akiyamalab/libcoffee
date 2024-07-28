from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt
from rdkit import Chem

from libcoffee.common.molbase import MolBase


class Mol(MolBase):
    """
    rdkit.Chem.Mol wrapper class
    """

    def __init__(self, mol: Chem.Mol):
        super().__init__(mol)

    @property
    def _atoms(self) -> tuple[Chem.Atom, ...]:
        return tuple(m for m in self._mol.GetAtoms())

    @property
    def isotopes(self) -> npt.NDArray[np.int_]:
        # a.GetIsotope() should return int but typing says return Any
        return np.array([a.GetIsotope() for a in self._atoms], dtype=np.int_)  # type: ignore  # Argument is not needed for GetIsotope()

    @isotopes.setter
    def isotopes(self, isotopes: npt.NDArray[np.int_]) -> None:
        if len(isotopes) != len(self._atoms):
            raise ValueError("Length of isotopes should be equal to the number of atoms")
        for i in range(len(self._atoms)):
            self._atoms[i].SetIsotope(int(isotopes[i]))

    @property
    def name(self) -> str:
        return self._mol.GetProp("_Name")

    @property
    def heavy_atom_indices(self) -> npt.NDArray[np.int_]:
        atomic_nums = [a.GetAtomicNum() for a in self._atoms]
        return np.where(np.array(atomic_nums) > 1)[0]

    def get_smiles(self, kekulize: bool = False) -> str:
        return Chem.MolToSmiles(self._mol, kekuleSmiles=kekulize)

    def get_attr(self, attr_name: str) -> Any:
        return self._mol.GetProp(attr_name)

    def has_attr(self, attr_name: str) -> bool:
        return self._mol.HasProp(attr_name)

    def get_coordinates(self, only_heavy_atom: bool = False) -> npt.NDArray[np.float_]:
        conf = self._mol.GetConformer()
        coords = np.array([conf.GetAtomPosition(i) for i in range(self._mol.GetNumAtoms())])
        if only_heavy_atom:
            coords = coords[self.heavy_atom_indices]
        return coords

    def extract_submol(self, atom_idxs: npt.NDArray[np.int_]) -> "Mol":
        rw_mol = Chem.RWMol(self.raw_mol)
        idx_remove_atoms = set(range(self.raw_mol.GetNumAtoms())) - set(atom_idxs)
        atomidxs = sorted(idx_remove_atoms)[::-1]
        for idx in atomidxs:
            rw_mol.RemoveAtom(idx)
        return Mol(rw_mol.GetMol())  # type: ignore  # Argument of GetMol() is not needed

    def merge(self, mol: "Mol", aps: tuple[int, int] | None = None) -> "Mol":
        raise NotImplementedError

    @classmethod
    def read_sdf(cls, file_path: Path) -> tuple["Mol", ...]:
        suppl = Chem.SDMolSupplier(str(file_path))
        return tuple(Mol(m) for m in suppl if m is not None)
