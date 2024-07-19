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
        self._mol = mol

    @property
    def atoms(self) -> list[pybel.Atom]:
        """
        Returns a list of atoms in the molecule
        """
        return self._mol.atoms

    @property
    def isotopes(self) -> npt.NDArray[np.int_]:
        """
        Returns a list of isotope numbers of atoms in the molecule
        """
        return np.array([a.isotope for a in self.atoms], dtype=np.int_)

    @isotopes.setter
    def isotopes(self, isotopes: npt.NDArray[np.int_]) -> None:
        """
        Sets the isotope numbers of atoms in the molecule
        Now this method is not implemented
        """
        if len(isotopes) != len(self.atoms):
            raise ValueError("Length of isotopes should be equal to the number of atoms")
        raise NotImplementedError

    @property
    def name(self):
        """
        Returns the name of the molecule
        """
        return self._mol.title

    @property
    def coordinates(self) -> npt.NDArray[np.float_]:
        return np.array([a.coords for a in self.atoms], dtype=np.float_)

    def get_smiles(self, kekulize: bool = False) -> str:
        """
        Returns the SMILES representation of the molecule
        """
        obConversion = OBConversion()
        obConversion.SetOutFormat("can")
        if kekulize:
            obConversion.SetOptions("k", obConversion.OUTOPTIONS)
        return obConversion.WriteString(self._mol.OBMol).split()[0]

    def get_attr(self, attr_name: str) -> Any:
        """
        Returns the value of the attribute with the given name
        """
        return self._mol.data[attr_name]

    def has_attr(self, attr_name: str) -> bool:
        """
        Returns True if the molecule has an attribute with the given name
        """
        return attr_name in self._mol.data
