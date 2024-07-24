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
        super().__init__(mol)


    @property
    def atoms(self) -> list[Chem.Atom]:
        return [m for m in self._mol.GetAtoms()]

    @property
    def isotopes(self) -> npt.NDArray[np.int_]:
        # a.GetIsotope() should return int but typing says return Any
        return np.array([a.GetIsotope() for a in self.atoms], dtype=np.int_)  # type: ignore

    @isotopes.setter
    def isotopes(self, isotopes: npt.NDArray[np.int_]) -> None:
        if len(isotopes) != len(self.atoms):
            raise ValueError("Length of isotopes should be equal to the number of atoms")
        raise NotImplementedError

    @property
    def name(self) -> str:
        return self._mol.GetProp("_Name")

    @property
    def heavy_atom_indices(self) -> npt.NDArray[np.int_]:
        raise NotImplementedError

    def get_smiles(self, kekulize: bool = False) -> str:
        return Chem.MolToSmiles(self._mol, kekuleSmiles=kekulize)

    def get_attr(self, attr_name: str) -> Any:
        return self._mol.GetProp(attr_name)

    def has_attr(self, attr_name: str) -> bool:
        return self._mol.HasProp(attr_name)

    def get_coordinates(self, only_heavy_atom: bool = False) -> npt.NDArray[np.float_]:
        raise NotImplementedError
