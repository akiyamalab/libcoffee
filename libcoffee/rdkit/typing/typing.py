from typing import TypeVar, NewType
from rdkit import Chem

T_RDKitMol = TypeVar("T_RDKitMol", bound=Chem.rdchem.Mol | Chem.rdchem.RWMol)
RDKitMol = NewType("RDKitMol", type[T_RDKitMol])
