from abc import ABCMeta, abstractmethod
from typing import Any

import numpy as np
import numpy.typing as npt

class MoleculeBase(metaclass=ABCMeta):
    @abstractmethod
    def __init__(self, mol): ...

    @property
    @abstractmethod
    def atoms(self) -> list: ...

    @property
    @abstractmethod
    def isotopes(self) -> npt.NDArray[np.int_]: ...
    
    @property
    @abstractmethod
    def name(self) -> str: ...

    @property
    @abstractmethod
    def smiles(self, kekulize: bool = False) -> str: ...

    @property
    @abstractmethod
    def coordinates(self) -> npt.NDArray[np.float_]: ...

    @abstractmethod
    def get_attr(self, attr_name: str) -> Any: ...
    
    @abstractmethod
    def has_attr(self, attr_name: str) -> bool: ...
    
    @abstractmethod
    def set_isotopes(self, isotopes: npt.NDArray[np.int_]) -> None: ...

    @property
    def center(self) -> npt.NDArray[np.float_]:
        return np.mean(self.coordinates, axis=0)
