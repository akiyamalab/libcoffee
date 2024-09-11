from typing import Final

from libcoffee.molecule.molbase import MolBase


def _eboxsize(mol: MolBase) -> int:
    """
    Calculate the box size for ligand with eBoxSize algorithm
    https://www.brylinski.org/eboxsize
    """
    _GY_BOX_RATIO: Final[float] = 0.23
    center = mol.get_coordinates().mean(axis=0)
    sq_gyration = ((mol.get_coordinates() - center) ** 2).sum(axis=1).mean()
    size = sq_gyration**0.5 / _GY_BOX_RATIO
    return int((size + 1) / 2) * 2  # round up to the nearest even number
