from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt
from openbabel import pybel
from openbabel.openbabel import OBConversion

from libcoffee.common.molbase import MolBase

from ._util import combine_two_mols


class Mol(MolBase):
    """
    openbabel.pybel.Molecule wrapper class
    """

    def __init__(self, mol: pybel.Molecule):
        super().__init__(mol)

    @property
    def _atoms(self) -> tuple[pybel.Atom, ...]:
        return tuple(self.raw_mol.atoms)

    @property
    def isotopes(self) -> npt.NDArray[np.int_]:
        return np.array([a.isotope for a in self._atoms], dtype=np.int_)

    @isotopes.setter
    def isotopes(self, isotopes: npt.NDArray[np.int_]) -> None:
        if len(isotopes) != len(self._atoms):
            raise ValueError("Length of isotopes should be equal to the number of atoms")
        raise NotImplementedError

    @property
    def name(self) -> str:
        return self.raw_mol.title

    @property
    def heavy_atom_indices(self) -> npt.NDArray[np.int_]:
        unsorted_indices = list(set(np.where(np.array([a.atomicnum for a in self.raw_mol.atoms]) > 1)[0]))
        return np.array(sorted(unsorted_indices), dtype=np.int_)

    def get_coordinates(self, only_heavy_atom: bool = False) -> npt.NDArray[np.float_]:
        coords = np.array([a.coords for a in self.raw_mol.atoms])
        if only_heavy_atom:
            coords = coords[self.heavy_atom_indices]
        return coords

    def get_smiles(self, kekulize: bool = False) -> str:
        obConversion = OBConversion()
        obConversion.SetOutFormat("can")
        if kekulize:
            obConversion.SetOptions("k", obConversion.OUTOPTIONS)
        return obConversion.WriteString(self.raw_mol.OBMol).split()[0]

    def get_attr(self, attr_name: str) -> Any:
        return self.raw_mol.data[attr_name]

    def has_attr(self, attr_name: str) -> bool:
        return attr_name in self.raw_mol.data

    def extract_submol(self, atom_idxs: npt.NDArray[np.int_]) -> "MolBase":
        raise NotImplementedError

    def merge(self, mol: "Mol", aps: tuple[int, int] | None = None) -> "Mol":
        natoms = len(self._atoms)
        ret = combine_two_mols(self.raw_mol, mol.raw_mol)
        if aps is not None:
            ap1, ap2 = aps
            ap1, ap2 = ap1 + 1, ap2 + natoms + 1  # +1: atom index starts from 1
            ret.OBMol.AddBond(ap1, ap2, 1)  # ap1, ap2の間に単結合を追加
        return Mol(ret)

    @classmethod
    def read_sdf(cls, file_path: Path) -> tuple["Mol", ...]:
        """
        Reads molecules from an SDF file and returns the molecule objects
        """
        molecules = list(pybel.readfile("sdf", str(file_path)))
        return tuple(cls(mol) for mol in molecules if mol is not None)
