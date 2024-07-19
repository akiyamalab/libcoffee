from abc import ABCMeta, abstractmethod
from typing import Any

import numpy as np
import numpy.typing as npt


class MolBase(metaclass=ABCMeta):
    """
    Wrapper class for molecule objects of RDKit and OpenBabel.
    The aim is to capture common attributes and methods of the two classes,
    and reduce the amount of code duplication in the two classes.
    """

    @abstractmethod
    def __init__(self, mol): ...

    @property
    @abstractmethod
    def atoms(self) -> list: ...

    @property
    @abstractmethod
    def isotopes(self) -> npt.NDArray[np.int_]: ...

    @isotopes.setter
    @abstractmethod
    def isotopes(self, isotopes: npt.NDArray[np.int_]) -> None: ...

    @property
    @abstractmethod
    def name(self) -> str: ...

    @property
    @abstractmethod
    def coordinates(self) -> npt.NDArray[np.float_]: ...

    @abstractmethod
    def get_smiles(self, kekulize: bool = False) -> str: ...

    @abstractmethod
    def get_attr(self, attr_name: str) -> Any: ...

    @abstractmethod
    def has_attr(self, attr_name: str) -> bool: ...

    @property
    def center(self) -> npt.NDArray[np.float_]:
        return np.mean(self.coordinates, axis=0)
