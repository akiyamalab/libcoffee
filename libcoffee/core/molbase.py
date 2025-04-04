from __future__ import annotations

import copy
import random
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Self

import numpy as np
import numpy.typing as npt


class MolBase(ABC):
    """
    Wrapper class for molecule objects of RDKit and OpenBabel.
    The aim is to capture common attributes and methods of the two classes,
    and reduce the amount of code duplication in the two classes.
    """

    def __init__(self: Self, mol: Any):
        self._mol = mol

    @property
    def raw_mol(self: Self) -> Any:
        """
        Returns the raw molecule object
        """
        return self._mol

    @property
    @abstractmethod
    def _atoms(self: Self) -> tuple[Any, ...]:
        """
        Returns a list of atoms in the molecule.
        This method should be private method because output is not consistent between RDKit and OpenBabel.
        """
        pass

    @_atoms.setter
    @abstractmethod
    def _atoms(self: Self, atoms: tuple[Any, ...]) -> None:
        """
        Sets the atoms of the molecule
        """
        pass

    @property
    @abstractmethod
    def bonds(self: Self) -> tuple[Any, ...]:
        """
        Returns a list of bonds in the molecule.
        """
        pass

    @property
    @abstractmethod
    def isotopes(self: Self) -> npt.NDArray[np.int32]:
        """
        Returns a list of isotope numbers of atoms in the molecule
        """
        pass

    @isotopes.setter
    @abstractmethod
    def isotopes(self: Self, isotopes: npt.NDArray[np.int32]) -> None:
        """
        Sets the isotope numbers of atoms in the molecule
        """
        pass

    @property
    @abstractmethod
    def name(self: Self) -> str:
        """
        Returns the name of the molecule
        """
        pass

    @name.setter
    @abstractmethod
    def name(self: Self, name: str) -> None:
        """
        Sets the name of the molecule
        """
        pass

    @property
    @abstractmethod
    def heavy_atom_indices(self: Self) -> npt.NDArray[np.intp]:
        """
        Returns the indices of heavy atoms in the molecule
        """
        pass

    @abstractmethod
    def get_smiles(self: Self, kekulize: bool = False) -> str:
        """
        Returns the SMILES representation of the molecule
        """
        pass

    @abstractmethod
    def get_attr(self: Self, attr_name: str) -> Any:
        """
        Returns the value of the attribute with the given name
        """
        pass

    @abstractmethod
    def set_attr(self: Self, attr_name: str, value: Any) -> None:
        """
        Sets the value of the attribute with the given name
        """
        pass

    @abstractmethod
    def has_attr(self, attr_name: str) -> bool:
        """
        Returns True if the molecule has an attribute with the given name
        """
        pass

    @abstractmethod
    def get_coordinates(self: Self, only_heavy_atom: bool = False) -> npt.NDArray[np.float64]:
        """
        Returns the coordinates of all atoms that make up the molecule.
        If only_heavy_atom is True, return the coordinates of heavy atoms only.
        """
        pass

    @abstractmethod
    def generate_coordinates(self: Self, temporary_add_hydrogens: bool = False) -> None:
        """
        Generates 3D coordinates of the molecule.
        """
        pass

    @abstractmethod
    def has_coordinates(self: Self) -> bool:
        """
        Returns True if the molecule has 3D coordinates.
        """
        pass

    @abstractmethod
    def remove_hydrogens(self: Self) -> MolBase:
        """
        Removes hydrogen atoms from the molecule.
        """
        pass

    @abstractmethod
    def add_hydrogens(self: Self) -> MolBase:
        """
        Adds hydrogen atoms to the molecule.
        """
        pass

    @abstractmethod
    def extract_submol(self: Self, atom_idxs: npt.NDArray[np.intp]) -> MolBase:
        """
        Extracts a substructure molecule from the original molecule with the given atom indices.
        """
        pass

    def center(self: Self, only_heavy_atom: bool = False) -> npt.NDArray[np.float64]:
        """
        Returns the center of the molecule.
        If only_heavy_atom is True, return the center of heavy atoms only.
        """
        return np.mean(self.get_coordinates(only_heavy_atom), axis=0)  # type: ignore[no-any-return]

    @abstractmethod
    def merge(self: Self, mol: Any, aps: tuple[int, int] | None = None) -> MolBase:
        """
        Merges the current molecule with another molecule.
        If aps (attachment points) is given, the two molecules are bonded at the given attachment points.
        Isotope numbers of atoms in the given molecule are updated to avoid conflicts.
        """
        pass

    def annotate_bonds_infomation_to_atoms(self: Self) -> None:
        """
        Annotates bond information to atoms
        """
        isotopes = self.isotopes
        for bond in self.bonds:
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()

            if isotopes[atom1.GetIdx()] != isotopes[atom2.GetIdx()]:
                rand = random.randint(0, 100000)
                if atom1.HasProp("attachment_point"):
                    atom1.SetProp("attachment_point", f"{atom1.GetProp("attachment_point")};{str(rand)}")
                else:
                    atom1.SetProp("attachment_point", str(rand))

                if atom2.HasProp("attachment_point"):
                    atom2.SetProp("attachment_point", f"{atom2.GetProp("attachment_point")};{str(rand)}")
                else:
                    atom2.SetProp("attachment_point", str(rand))

    def split(self: Self) -> tuple[MolBase, ...]:
        """
        Splits the molecule into fragments based on isotopes information
        and returns the fragment molecules
        """
        frag_idxs = np.unique(self.isotopes)
        frags: list[MolBase] = []
        for idx in frag_idxs:
            atom_idxs = np.where(self.isotopes == idx)[0]
            frags.append(self.extract_submol(atom_idxs))
        return tuple(frags)

    def deep_copy(self: Self) -> MolBase:
        """
        Returns a deep copy of the molecule object
        """
        new_mol = copy.deepcopy(self.raw_mol)
        return self.__class__(new_mol)

    @classmethod
    @abstractmethod
    def reconstruct_from_fragments(cls, frags: tuple[MolBase, ...]) -> MolBase:
        """
        Reconstructs a molecule from fragments
        """
        pass

    @classmethod
    @abstractmethod
    def from_smiles(cls, smiles: str) -> MolBase:
        """
        Creates a molecule object from a SMILES string
        """
        pass

    @classmethod
    @abstractmethod
    def read_sdf(cls, file_path: Path) -> tuple[MolBase, ...]:
        """
        Reads molecules from an SDF file and returns the molecule objects
        """
        pass

    @classmethod
    @abstractmethod
    def write_sdf(cls, file_path: Path, mols: tuple[MolBase, ...]) -> None:
        """
        Writes the given molecules to an SDF file
        """
        pass
