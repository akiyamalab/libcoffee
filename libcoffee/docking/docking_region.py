from dataclasses import dataclass
from typing import Final

from libcoffee.core.molbase import MolBase


@dataclass(frozen=True)
class DockingRegion:
    center: tuple[float, float, float] = (0.0, 0.0, 0.0)
    size: tuple[int, int, int] = (10, 10, 10)
    pitch: tuple[float, float, float] = (1.0, 1.0, 1.0)

    @staticmethod
    def from_molecule(mol: MolBase, pitch: tuple[float, float, float] = (1.0, 1.0, 1.0)) -> "DockingRegion":
        center = mol.get_coordinates(only_heavy_atom=True).mean(axis=0)
        size = (_eboxsize(mol),) * 3
        return DockingRegion(center=center, size=size, pitch=pitch)


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
