from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt


class MolBase(ABC):
    """
    Wrapper class for molecule objects of RDKit and OpenBabel.
    The aim is to capture common attributes and methods of the two classes,
    and reduce the amount of code duplication in the two classes.
    """

    def __init__(self, mol):
        self._mol = mol

    @property
    def raw_mol(self):
        """
        Returns the raw molecule object
        """
        return self._mol

    @property
    @abstractmethod
    def _atoms(self) -> tuple:
        """
        Returns a list of atoms in the molecule.
        This method should be private method because output is not consistent between RDKit and OpenBabel.
        """
        pass

    @property
    @abstractmethod
    def isotopes(self) -> npt.NDArray[np.int_]:
        """
        Returns a list of isotope numbers of atoms in the molecule
        """
        pass

    @isotopes.setter
    @abstractmethod
    def isotopes(self, isotopes: npt.NDArray[np.int_]) -> None:
        """
        Sets the isotope numbers of atoms in the molecule
        """
        pass

    @property
    @abstractmethod
    def name(self) -> str:
        """
        Returns the name of the molecule
        """
        pass

    @property
    @abstractmethod
    def heavy_atom_indices(self) -> npt.NDArray[np.int_]:
        """
        Returns the indices of heavy atoms in the molecule
        """
        pass

    @abstractmethod
    def get_smiles(self, kekulize: bool = False) -> str:
        """
        Returns the SMILES representation of the molecule
        """
        pass

    @abstractmethod
    def get_attr(self, attr_name: str) -> Any:
        """
        Returns the value of the attribute with the given name
        """
        pass

    @abstractmethod
    def has_attr(self, attr_name: str) -> bool:
        """
        Returns True if the molecule has an attribute with the given name
        """
        pass

    @abstractmethod
    def get_coordinates(self, only_heavy_atom: bool = False) -> npt.NDArray[np.float_]:
        """
        Returns the coordinates of all atoms that make up the molecule.
        If only_heavy_atom is True, return the coordinates of heavy atoms only.
        """
        pass

    @abstractmethod
    def extract_submol(self, atom_idxs: npt.NDArray[np.int_]) -> "MolBase":
        """
        Extracts a substructure molecule from the original molecule with the given atom indices.
        """
        pass

    def center(self, only_heavy_atom: bool = False) -> npt.NDArray[np.float_]:
        """
        Returns the center of the molecule.
        If only_heavy_atom is True, return the center of heavy atoms only.
        """
        return np.mean(self.get_coordinates(only_heavy_atom), axis=0)

    @abstractmethod
    def merge(self, mol, aps: tuple[int, int] | None = None):
        """
        Merges the current molecule with another molecule.
        If aps (attachment points) is given, the two molecules are bonded at the given attachment points.
        Isotope numbers of atoms in the given molecule are updated to avoid conflicts.
        """
        pass

    def split(self) -> tuple:
        """
        Splits the molecule into fragments based on isotopes information
        and returns the fragment molecules
        """
        frag_idxs = np.unique(self.isotopes)
        frags = []
        for idx in frag_idxs:
            atom_idxs = np.where(self.isotopes == idx)[0]
            frags.append(self.extract_submol(atom_idxs))
        return tuple(frags)

    @classmethod
    @abstractmethod
    def read_sdf(cls, file_path: Path) -> tuple:
        """
        Reads molecules from an SDF file and returns the molecule objects
        """
        pass
