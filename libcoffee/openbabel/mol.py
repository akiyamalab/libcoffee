from typing import Any
from openbabel import pybel
from libcoffee.common.molbase import MolBase
from openbabel.openbabel import OBConversion
import numpy as np
import numpy.typing as npt


class Molecule(MolBase):
    """
    openbabel.pybel.Molecule wrapper class
    """

    def __init__(self, mol: pybel.Molecule):
        super().__init__(mol)

    @property
    def atoms(self) -> list[pybel.Atom]:
        return self.raw_mol.atoms

    @property
    def isotopes(self) -> npt.NDArray[np.int_]:
        return np.array([a.isotope for a in self.atoms], dtype=np.int_)

    @isotopes.setter
    def isotopes(self, isotopes: npt.NDArray[np.int_]) -> None:
        if len(isotopes) != len(self.atoms):
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

    def extract_submol(self, atom_idxs: list[int]) -> "MolBase":
        raise NotImplementedError