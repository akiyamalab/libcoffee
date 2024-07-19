from typing import Any
from libcoffee.common.molbase import MolBase
import numpy as np
import numpy.typing as npt
from rdkit import Chem


class Mol(MolBase):
    """
    rdkit.Chem.Mol wrapper class
    """

    def __init__(self, mol: Chem.Mol):
        self._mol = mol

    @property
    def atoms(self) -> list[Chem.Atom]:
        """
        Returns a list of atoms in the molecule
        """
        return [m for m in self._mol.GetAtoms()]

    @property
    def isotopes(self) -> npt.NDArray[np.int_]:
        """
        Returns a list of isotope numbers of atoms in the molecule
        """
        # a.GetIsotope() should return int but typing says return Any
        return np.array([a.GetIsotope() for a in self.atoms], dtype=np.int_)  # type: ignore

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
        return self._mol.GetProp("_Name")

    @property
    def coordinates(self) -> npt.NDArray[np.float_]:
        """
        Returns the coordinates of atoms in the molecule
        """
        # a.GetPos() should return Tuple[float, float, float] but typing says return Any
        return np.array([a.GetPos() for a in self.atoms], dtype=np.float_)  # type: ignore

    def get_smiles(self, kekulize: bool = False) -> str:
        """
        Returns the SMILES representation of the molecule
        """
        return Chem.MolToSmiles(self._mol, kekuleSmiles=kekulize)

    def get_attr(self, attr_name: str) -> Any:
        """
        Returns the value of the attribute with the given name
        """
        return self._mol.GetProp(attr_name)

    def has_attr(self, attr_name: str) -> bool:
        """
        Returns True if the molecule has an attribute with the given name
        """
        return self._mol.HasProp(attr_name)
