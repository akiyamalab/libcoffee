from libcoffee.rdkit.mol import Mol
from rdkit import Chem
import pytest

class TestMol:
  @pytest.fixture
  def init(self):
    self.mol = Mol(Chem.MolFromSmiles('c1ccccc1'))

  def test_atoms(self, init):
    assert len(self.mol.atoms) == 6
    
    