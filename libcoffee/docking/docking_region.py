from dataclasses import dataclass
from typing import Any, Final

import numpy as np
import numpy.typing as npt

from libcoffee.core.molbase import MolBase


@dataclass(frozen=True, init=False)
class DockingRegion:
    center: tuple[float, float, float]
    size: tuple[int, int, int]
    pitch: tuple[float, float, float]

    def __init__(
        self,
        center: npt.ArrayLike = (0.0, 0.0, 0.0),
        size: npt.ArrayLike = (10, 10, 10),
        pitch: npt.ArrayLike = (1.0, 1.0, 1.0),
    ) -> None:
        # accept any array-like input, but always store immutable tuples
        object.__setattr__(self, "center", _to_float3(center))
        object.__setattr__(self, "size", _to_int3(size))
        object.__setattr__(self, "pitch", _to_float3(pitch))

    @staticmethod
    def from_molecule(mol: MolBase, pitch: npt.ArrayLike = (1.0, 1.0, 1.0)) -> "DockingRegion":
        size = (_eboxsize(mol),) * 3
        return DockingRegion(center=mol.center(only_heavy_atom=True), size=size, pitch=pitch)


def _check_3d_vector(arr: npt.NDArray[Any], v: npt.ArrayLike) -> None:
    if arr.shape != (3,):
        raise ValueError(f"Must be a sequence of 3 numbers, but got shape {arr.shape}: {v!r}")


def _to_float3(v: npt.ArrayLike) -> tuple[float, float, float]:
    arr = np.asarray(v, dtype=np.float64)
    _check_3d_vector(arr, v)
    x, y, z = arr
    return float(x), float(y), float(z)


def _to_int3(v: npt.ArrayLike) -> tuple[int, int, int]:
    arr = np.asarray(v)
    _check_3d_vector(arr, v)
    if not np.issubdtype(arr.dtype, np.integer):
        raise ValueError(f"Must consist of integers, but got dtype {arr.dtype}: {v!r}")
    x, y, z = arr
    return int(x), int(y), int(z)


def _eboxsize(mol: MolBase) -> int:
    """
    Calculate the box size for ligand with eBoxSize algorithm
    https://www.brylinski.org/eboxsize
    """
    _GY_BOX_RATIO: Final[float] = 0.23
    center = mol.get_coordinates(only_heavy_atom=True).mean(axis=0)
    sq_gyration = ((mol.get_coordinates(only_heavy_atom=True) - center) ** 2).sum(axis=1).mean()
    size = sq_gyration**0.5 / _GY_BOX_RATIO
    return int((size + 1) / 2) * 2  # round up to the nearest even number
