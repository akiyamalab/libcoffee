from typing import Any
from openbabel import pybel
from libcoffee.common.molbase import MolBase
from openbabel.openbabel import OBConversion
import numpy as np
import numpy.typing as npt

class Molecule(MolBase):

    def __init__(self, mol: pybel.Molecule):
        self._mol = mol

    @property
    def atoms(self) -> list[pybel.Atom]:
        return self._mol.atoms

    @property
    def isotopes(self) -> npt.NDArray[np.int_]:
        return np.array([a.isotope for a in self.atoms], dtype=np.int_)

    @property
    def name(self):
        return self._mol.title

    @property
    def smiles(self, kekulize: bool = False) -> str:
        obConversion = OBConversion()
        obConversion.SetOutFormat("can")
        if kekulize:
            obConversion.SetOptions("k", obConversion.OUTOPTIONS)
        return obConversion.WriteString(self._mol.OBMol).split()[0]

    @property
    def coordinates(self) -> npt.NDArray[np.float_]:
        return np.array([a.coords for a in self.atoms], dtype=np.float_)

    def get_attr(self, attr_name: str) -> Any:
        return self._mol.data[attr_name]
    
    def has_attr(self, attr_name: str) -> bool:
        return attr_name in self._mol.data
    
    def set_isotopes(self, isotopes: npt.NDArray[np.int_]) -> None:
        if len(isotopes) != len(self.atoms):
            raise ValueError("Length of isotopes should be equal to the number of atoms")
        raise NotImplementedError
    