from typing import Any
from libcoffee.common.molecule_base import MoleculeBase
import numpy as np
import numpy.typing as npt
from rdkit import Chem
from rdkit.Chem import Mol, Atom

class Molecule(MoleculeBase):

    def __init__(self, mol: Mol):
        self._mol = mol

    @property
    def atoms(self) -> list[Atom]:
        return [m for m in self._mol.GetAtoms()]

    @property
    def isotopes(self) -> npt.NDArray[np.int_]:
        # a.GetIsotope() should return int but typing says return Any
        return np.array([a.GetIsotope() for a in self.atoms], dtype=np.int_) # type: ignore

    @property
    def name(self):
        return self._mol.GetProp("_Name")

    @property
    def smiles(self, kekulize: bool = False) -> str:
        return Chem.MolToSmiles(self._mol, kekuleSmiles=kekulize)

    @property
    def coordinates(self) -> npt.NDArray[np.float_]:
        # a.GetPos() should return Tuple[float, float, float] but typing says return Any
        return np.array([a.GetPos() for a in self.atoms], dtype=np.float_) # type: ignore

    def get_attr(self, attr_name: str) -> Any:
        return self._mol.GetProp(attr_name)
    
    def has_attr(self, attr_name: str) -> bool:
        return self._mol.HasProp(attr_name)

    def set_isotopes(self, isotopes: npt.NDArray[np.int_]) -> None:
        if len(isotopes) != len(self.atoms):
            raise ValueError("Length of isotopes should be equal to the number of atoms")
        raise NotImplementedError