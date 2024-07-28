from openbabel import pybel


def combine_two_mols(mol1: pybel.Molecule, mol2: pybel.Molecule) -> pybel.Molecule:
    mol = mol1.clone
    cnt = len(mol.atoms)
    obmol = mol.OBMol
    mol2 = mol2.clone
    for i in range(len(mol2.atoms)):
        obatom = mol2.atoms[i].OBAtom
        obmol.AddAtom(obatom)
    for i in range(mol2.OBMol.NumBonds()):
        bond = mol2.OBMol.GetBond(i)
        obmol.AddBond(bond.GetBeginAtomIdx() + cnt, bond.GetEndAtomIdx() + cnt, bond.GetBondOrder(), bond.GetFlags())
    return mol
