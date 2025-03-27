from typing import Any, Self

import numpy as np
import pytest
from openbabel import pybel

from libcoffee.openbabeltools import PybelMol


class TestMolFromSmiles:
    @pytest.fixture
    def init(self: Self) -> None:
        self.mol = PybelMol(pybel.readstring("smi", "c1ccccc1 benzene"))

    def test_isotopes(self: Self, init: Any) -> None:
        assert np.all(self.mol.isotopes == 0)

    @pytest.mark.skip(reason="isotopes.setter is not implemented")
    def test_isotopes_setter(self: Self, init: Any) -> None:
        self.mol.isotopes = np.array([1, 2, 3, 4, 5, 6])
        assert np.all(self.mol.isotopes == np.array([1, 2, 3, 4, 5, 6]))

    def test_name(self: Self, init: Any) -> None:
        assert self.mol.name == "benzene"

    def test_heavy_atom_indices(self: Self, init: Any) -> None:
        assert np.all(self.mol.heavy_atom_indices == np.array([0, 1, 2, 3, 4, 5]))

    def test_get_smiles(self: Self, init: Any) -> None:
        assert self.mol.get_smiles() == "c1ccccc1"
