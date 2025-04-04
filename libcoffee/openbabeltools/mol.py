from __future__ import annotations

from pathlib import Path
from typing import Any, Self

import numpy as np
import numpy.typing as npt
from openbabel import pybel
from openbabel.openbabel import OBConversion

from libcoffee.core.molbase import MolBase

from ._util import combine_two_mols


class PybelMol(MolBase):
    """
    openbabel.pybel.Molecule wrapper class
    """

    def __init__(self, mol: pybel.Molecule):
        if not isinstance(mol, pybel.Molecule):
            raise ValueError(
                f"The argument of PybelMol constructor should be an instance of openbabel.pybel.Molecule, not {type(mol)}"
            )
        super().__init__(mol)

    @property
    def _atoms(self) -> tuple[pybel.Atom, ...]:
        return tuple(self.raw_mol.atoms)

    @_atoms.setter
    def _atoms(self, atoms: tuple[pybel.Atom, ...]) -> None:
        if len(atoms) != len(self._atoms):
            raise ValueError("Length of atoms should be equal to the number of atoms")
        self.raw_mol.atoms = atoms

    @property
    def bonds(self) -> tuple[Any, ...]:
        raise NotImplementedError

    @property
    def isotopes(self) -> npt.NDArray[np.int32]:
        return np.array([a.isotope for a in self._atoms], dtype=np.int32)

    @isotopes.setter
    def isotopes(self, isotopes: npt.NDArray[np.int32]) -> None:
        if len(isotopes) != len(self._atoms):
            raise ValueError("Length of isotopes should be equal to the number of atoms")
        raise NotImplementedError

    @property
    def name(self) -> str:
        return self.raw_mol.title  # type: ignore[no-any-return]

    @name.setter
    def name(self, name: str) -> None:
        self.raw_mol.title = name

    @property
    def heavy_atom_indices(self) -> npt.NDArray[np.intp]:
        unsorted_indices: list[int] = list(set(np.where(np.array([a.atomicnum for a in self.raw_mol.atoms]) > 1)[0]))
        return np.array(sorted(unsorted_indices), dtype=np.intp)

    def get_coordinates(self, only_heavy_atom: bool = False) -> npt.NDArray[np.float64]:
        coords = np.array([a.coords for a in self.raw_mol.atoms])
        if only_heavy_atom:
            coords = coords[self.heavy_atom_indices]
        return coords

    def generate_coordinates(self, temporary_add_hydrogens: bool = False) -> None:
        raise NotImplementedError

    def has_coordinates(self) -> bool:
        return all(a.coords is not None for a in self.raw_mol.atoms)

    def add_hydrogens(self) -> PybelMol:
        self.raw_mol.addh()
        return self

    def remove_hydrogens(self) -> PybelMol:
        self.raw_mol.removeh()
        return self

    def get_smiles(self, kekulize: bool = False) -> str:
        obConversion = OBConversion()
        obConversion.SetOutFormat("can")
        if kekulize:
            obConversion.SetOptions("k", obConversion.OUTOPTIONS)
        return obConversion.WriteString(self.raw_mol.OBMol).split()[0]  # type: ignore[no-any-return]

    def get_attr(self, attr_name: str) -> Any:
        return self.raw_mol.data[attr_name]

    def set_attr(self, attr_name: str, value: Any) -> None:
        self.raw_mol.data[attr_name] = value

    def has_attr(self, attr_name: str) -> bool:
        return attr_name in self.raw_mol.data

    def extract_submol(self, atom_idxs: npt.NDArray[np.intp]) -> MolBase:
        raise NotImplementedError

    def merge(self, mol: PybelMol, aps: tuple[int, int] | None = None) -> PybelMol:
        natoms = len(self._atoms)
        ret = combine_two_mols(self.raw_mol, mol.raw_mol)
        if aps is not None:
            ap1, ap2 = aps
            ap1, ap2 = ap1 + 1, ap2 + natoms + 1  # +1: atom index starts from 1
            ret.OBMol.AddBond(ap1, ap2, 1)  # ap1, ap2の間に単結合を追加
        return PybelMol(ret)

    @classmethod
    def reconstruct_from_fragments(cls, frags: tuple[PybelMol, ...]) -> PybelMol:  # type: ignore[override]
        raise NotImplementedError

    @classmethod
    def from_smiles(cls, smiles: str) -> PybelMol:
        """
        Generates a molecule object from SMILES
        """
        return cls(pybel.readstring("smi", smiles))

    @classmethod
    def read_sdf(cls, file_path: Path) -> tuple[PybelMol, ...]:
        """
        Reads molecules from an SDF file and returns the molecule objects
        """
        molecules = list(pybel.readfile("sdf", str(file_path)))
        return tuple(cls(mol) for mol in molecules if mol is not None)

    @classmethod
    def write_sdf(cls, file_path: Path, mols: tuple[MolBase, ...]) -> None:
        """
        Writes the given molecules to an SDF file
        """
        with open(file_path, "w") as f:
            for mol in mols:
                f.write(mol.raw_mol.write("sdf"))
