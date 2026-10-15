from abc import ABC, abstractmethod
from dataclasses import dataclass
from functools import lru_cache

import numpy as np
import numpy.typing as npt
from rdkit import Chem, DataStructs
from rdkit.Chem import Crippen, rdFingerprintGenerator, rdMolAlign, rdMolDescriptors, rdShapeHelpers

from libcoffee.rdkittools.mol import RDKitMol


class SimilarityBase(ABC):
    """
    Base class for similarity calculation methods.
    Method-specific options are given to the constructor, and calc() has the same signature for all methods,
    in the same manner as RDKit fingerprint generators (e.g. GetMorganGenerator(radius=2).GetFingerprint(mol)).
    """

    @abstractmethod
    def calc(self, prb_mol: RDKitMol, ref_mol: RDKitMol, prealigned: bool = False) -> float:
        """
        Calculate the similarity between a probe molecule and a reference molecule.
        The argument order follows the RDKit alignment functions (e.g. rdMolAlign.AlignMol(prbMol, refMol)).

        Parameters
        ----------
        prb_mol : RDKitMol
            Probe molecule, which is aligned to ref_mol if needed. It is never modified.
        ref_mol : RDKitMol
            Reference molecule. It is never modified.
        prealigned : bool, default=False
            Whether prb_mol is already aligned to ref_mol.
            If False, methods that depend on 3D poses align a copy of prb_mol to ref_mol before the calculation.
            Methods that do not depend on 3D poses ignore this argument.
        """
        pass


@dataclass(frozen=True)
class SimilarityNunobe2024(SimilarityBase):
    """
    Calculate the similarity between two mols using the Nunobe 2024 similarity calculation method.

    Parameters
    ----------
    radius : int, default=2
        Radius of the Morgan fingerprint.
    bits : int, default=2048
        Number of bits of the Morgan fingerprint.
    """

    radius: int = 2
    bits: int = 2048

    # Fingerprints are cached with a binary snapshot of a molecule (Chem.Mol.ToBinary) as the key,
    # so a molecule modified after the first calculation is never served stale fingerprints.
    # This relies on ToBinary containing all information that MQN and Morgan fingerprints depend on.

    @classmethod
    @lru_cache(maxsize=100_000)
    def __calc_mqn_fps(cls, mol_binary: bytes) -> npt.NDArray[np.float64]:
        mol = Chem.Mol(mol_binary)  # type: ignore[call-overload]  # rdkit-stubs lacks the bytes overload
        return np.array(rdMolDescriptors.MQNs_(mol), dtype=np.float64)

    @classmethod
    def __calc_mqn_sim(cls, prb_binary: bytes, ref_binary: bytes) -> float:
        prb_mqn = cls.__calc_mqn_fps(prb_binary)
        ref_mqn = cls.__calc_mqn_fps(ref_binary)
        similarity: float = 1 / (1 + np.sum(np.abs(prb_mqn - ref_mqn)) / 42)
        return similarity

    # self is part of the cache key; instances with the same options are equal, so they share cached results.
    @lru_cache(maxsize=100_000)
    def __calc_morgan_fps(self, mol_binary: bytes) -> DataStructs.ExplicitBitVect:
        mol = Chem.Mol(mol_binary)  # type: ignore[call-overload]  # rdkit-stubs lacks the bytes overload
        mfgen = rdFingerprintGenerator.GetMorganGenerator(radius=self.radius, fpSize=self.bits)
        fp: DataStructs.ExplicitBitVect = mfgen.GetFingerprint(mol)
        return fp

    def __calc_morgan_sim(self, prb_binary: bytes, ref_binary: bytes) -> float:
        prb_morgan = self.__calc_morgan_fps(prb_binary)
        ref_morgan = self.__calc_morgan_fps(ref_binary)
        return DataStructs.TanimotoSimilarity(prb_morgan, ref_morgan)

    def calc(self, prb_mol: RDKitMol, ref_mol: RDKitMol, prealigned: bool = False) -> float:
        """
        Calculate the similarity as the average of MQN similarity and Morgan fingerprint Tanimoto similarity.
        The result is symmetric with respect to prb_mol and ref_mol.

        Parameters
        ----------
        prb_mol : RDKitMol
            Probe molecule.
        ref_mol : RDKitMol
            Reference molecule.
        prealigned : bool, default=False
            Ignored, because this method does not depend on 3D poses.
        """
        prb_binary = prb_mol.raw_mol.ToBinary()
        ref_binary = ref_mol.raw_mol.ToBinary()
        mqn_sim = self.__calc_mqn_sim(prb_binary, ref_binary)
        morgan_sim = self.__calc_morgan_sim(prb_binary, ref_binary)
        return (mqn_sim + morgan_sim) / 2


@dataclass(frozen=True)
class SimilarityYoneyama2025(SimilarityBase):
    """
    Calculate the similarity between two mols using the Yoneyama 2025 similarity calculation method.
    """

    @classmethod
    def align_mols(cls, prb_mol: RDKitMol, ref_mol: RDKitMol) -> None:
        """
        Align prb_mol to ref_mol in place with Crippen O3A,
        in the same manner as rdMolAlign.GetCrippenO3A(prbMol, refMol).
        The coordinates of prb_mol are modified; ref_mol is not modified.
        """
        prb_contribs = Crippen.rdMolDescriptors._CalcCrippenContribs(prb_mol.raw_mol)
        ref_contribs = Crippen.rdMolDescriptors._CalcCrippenContribs(ref_mol.raw_mol)

        o3a = rdMolAlign.GetCrippenO3A(prb_mol.raw_mol, ref_mol.raw_mol, prb_contribs, ref_contribs, maxIters=100)
        o3a.Align()

        match = o3a.Matches()
        if len(match) == 0:
            if len(prb_mol.raw_mol.GetAtoms()) == 1 or len(ref_mol.raw_mol.GetAtoms()) == 1:
                rdMolAlign.AlignMol(prb_mol.raw_mol, ref_mol.raw_mol, atomMap=[(0, 0)])
            else:
                rdMolAlign.AlignMol(prb_mol.raw_mol, ref_mol.raw_mol, atomMap=[(0, 0), (1, 1)])

    @classmethod
    def __calc_shape_tanimoto(cls, mol1: RDKitMol, mol2: RDKitMol) -> float:
        return 1 - rdShapeHelpers.ShapeTanimotoDist(mol1.raw_mol, mol2.raw_mol)

    def calc(self, prb_mol: RDKitMol, ref_mol: RDKitMol, prealigned: bool = False) -> float:
        """
        Calculate the shape Tanimoto similarity between a probe molecule and a reference molecule.
        Unless prealigned is True, the result is not symmetric, because prb_mol is aligned to ref_mol.

        Parameters
        ----------
        prb_mol : RDKitMol
            Probe molecule with 3D coordinates. It is never modified.
        ref_mol : RDKitMol
            Reference molecule with 3D coordinates. It is never modified.
        prealigned : bool, default=False
            Whether prb_mol is already aligned to ref_mol.
            If False, a copy of prb_mol is aligned to ref_mol with Crippen O3A before the calculation.
        """
        if not prb_mol.has_coordinates() or not ref_mol.has_coordinates():
            raise ValueError("Molecules must have coordinates to calculate the similarity.")

        aligned_prb = prb_mol
        if not prealigned:
            aligned_prb = prb_mol.deep_copy()  # type: ignore[assignment]
            self.align_mols(aligned_prb, ref_mol)
        return self.__calc_shape_tanimoto(aligned_prb, ref_mol)


def calc_similarities(
    prb_mols: tuple[RDKitMol, ...],
    ref_mol: RDKitMol,
    criteria: SimilarityBase,
    attribution: int,
    prealigned: bool = False,
) -> npt.NDArray[np.float64]:
    """
    Calculate the similarities between probe molecules and a reference molecule.

    Parameters
    ----------
    prb_mols : tuple[RDKitMol, ...]
        Probe molecules. They are never modified.
    ref_mol : RDKitMol
        Reference molecule. It is never modified.
    criteria : SimilarityBase
        Similarity calculation method with its options, e.g. SimilarityNunobe2024(radius=3).
    attribution : int
        If 0, only the maximum similarity is kept and the others are set to 0.
    prealigned : bool, default=False
        Whether prb_mols are already aligned to ref_mol. See SimilarityBase.calc.
    """
    similarities = np.array([criteria.calc(prb_mol, ref_mol, prealigned) for prb_mol in prb_mols], dtype=np.float64)
    if attribution == 0:
        max_sim_index = np.argmax(similarities)
        similarities = np.array(
            [sim if i == max_sim_index else 0 for i, sim in enumerate(similarities)], dtype=np.float64
        )
    return similarities
