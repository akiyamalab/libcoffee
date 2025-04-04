from abc import ABC, abstractmethod
from functools import lru_cache

import numpy as np
import numpy.typing as npt
from rdkit import Chem, DataStructs
from rdkit.Chem import Crippen, rdFingerprintGenerator, rdMolAlign, rdMolDescriptors, rdShapeHelpers

from libcoffee.rdkittools.mol import RDKitMol

# attributetion all max liratery


class SimilarityBase(ABC):
    """
    Base class for similarity calculation methods.
    """

    @classmethod
    @abstractmethod
    def calc(cls, mol1: RDKitMol, mol2: RDKitMol, prealigned: bool = True) -> float:
        """
        Calculate the similarity between two mols.
        """
        pass


class SimilarityNunobe2024(SimilarityBase):
    """
    Calculate the similarity between two mols using the Nunobe 2024 similarity calculation method.
    """

    @classmethod
    @lru_cache(maxsize=10000)
    def __calc_mqn_fps(cls, mol: RDKitMol) -> npt.NDArray[np.float64]:
        return np.array(rdMolDescriptors.MQNs_(mol.raw_mol), dtype=np.float64)

    @classmethod
    def __calc_mqn_sim(cls, mol1: RDKitMol, mol2: RDKitMol) -> float:
        mqn1 = cls.__calc_mqn_fps(mol1)
        mqn2 = cls.__calc_mqn_fps(mol2)
        similarity: float = 1 / (1 + np.sum(np.abs(mqn1 - mqn2)) / 42)
        return similarity

    @classmethod
    @lru_cache(maxsize=10000)
    def __calc_morgan_fps(cls, mol: RDKitMol, radius=2, bits=2048):  # type: ignore
        mfgen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=bits)
        return mfgen.GetFingerprint(mol.raw_mol)

    @classmethod
    def __calc_morgan_sim(cls, mol1: RDKitMol, mol2: RDKitMol, radius: int = 2, bits: int = 2048) -> float:
        morgan_fps1 = cls.__calc_morgan_fps(mol1, radius, bits)
        morgan_fps2 = cls.__calc_morgan_fps(mol2, radius, bits)
        return DataStructs.TanimotoSimilarity(morgan_fps1, morgan_fps2)

    @classmethod
    def calc(cls, mol1: RDKitMol, mol2: RDKitMol, prealigned: bool = True) -> float:
        return (cls.__calc_mqn_sim(mol1, mol2) + cls.__calc_morgan_sim(mol1, mol2)) / 2


class SimilarityYoneyama2025(SimilarityBase):
    """
    Calculate the similarity between two mols using the Yoneyama 2025 similarity calculation method.
    """

    @classmethod
    def align_mols(cls, ref_mol: RDKitMol, target_mol: RDKitMol) -> None:
        contribs1 = Crippen.rdMolDescriptors._CalcCrippenContribs(ref_mol.raw_mol)
        contribs2 = Crippen.rdMolDescriptors._CalcCrippenContribs(target_mol.raw_mol)

        o3a = rdMolAlign.GetCrippenO3A(ref_mol.raw_mol, target_mol.raw_mol, contribs1, contribs2, maxIters=100)
        o3a.Align()

        match = o3a.Matches()
        if len(match) == 0:
            if len(ref_mol.raw_mol.GetAtoms()) == 1 or len(target_mol.raw_mol.GetAtoms()) == 1:
                rdMolAlign.AlignMol(ref_mol.raw_mol, target_mol.raw_mol, atomMap=[(0, 0)])
            else:
                rdMolAlign.AlignMol(ref_mol.raw_mol, target_mol.raw_mol, atomMap=[(0, 0), (1, 1)])

    @classmethod
    @lru_cache(maxsize=10000)
    def __calc_shape_tanimoto(cls, mol1: RDKitMol, mol2: RDKitMol) -> float:
        return 1 - rdShapeHelpers.ShapeTanimotoDist(mol1.raw_mol, mol2.raw_mol)

    @classmethod
    def calc(cls, mol1: RDKitMol, mol2: RDKitMol, prealigned: bool = True) -> float:
        if not mol1.has_coordinates() or not mol2.has_coordinates():
            raise ValueError("Molecules must have coordinates to calculate the similarity.")

        mol2_deepcopy = mol2.deep_copy()

        if prealigned:
            cls.align_mols(mol1, mol2_deepcopy)  # type: ignore[arg-type]

        return cls.__calc_shape_tanimoto(mol1, mol2_deepcopy)


def calc_similarities(
    ref_mol: RDKitMol, target_mols: tuple[RDKitMol], criteria_name: str, attribution: int, prealigned: bool = True
) -> npt.NDArray[np.float64]:
    """
    Calculate the similarities between a reference molecule and target molecules.
    """

    criteria: type[SimilarityBase]
    if criteria_name == "nunobe2024":
        criteria = SimilarityNunobe2024
    elif criteria_name == "yoneyama2025":
        criteria = SimilarityYoneyama2025
    else:
        raise ValueError(f"Invalid criteria name: {criteria_name}")

    similarities = np.array(
        [criteria.calc(ref_mol, target_mol, prealigned) for target_mol in target_mols], dtype=np.float64
    )
    if attribution == 0:
        max_sim_index = np.argmax(similarities)
        similarities = np.array(
            [sim if i == max_sim_index else 0 for i, sim in enumerate(similarities)], dtype=np.float64
        )
    return similarities
