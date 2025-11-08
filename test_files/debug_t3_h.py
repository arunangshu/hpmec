"""Check hydrogen degrees in t3 molecule"""
from calculator import load_xyz, infer_topology

mol = load_xyz('t3.xyz')
infer_topology(mol)
rdkit_mol = mol.rdkit_mol

print("=" * 70)
print("HYDROGEN ANALYSIS FOR T3")
print("=" * 70)
print()

from collections import Counter

if rdkit_mol:
    h_degrees = []
    for atom in rdkit_mol.GetAtoms():
        if atom.GetSymbol() == 'H':
            idx = atom.GetIdx()
            degree = atom.GetDegree()
            neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
            h_degrees.append(degree)
            print(f"Atom {idx:2d} (H): Degree={degree}, Neighbors={neighbors}")
    
    print()
    print("Degree distribution:")
    for deg, count in Counter(h_degrees).items():
        print(f"  Degree {deg}: {count} hydrogens")
    
    print()
    print("Testing [#1;D1] pattern:")
    from rdkit import Chem
    patt = Chem.MolFromSmarts('[#1;D1]')
    matches = rdkit_mol.GetSubstructMatches(patt)
    print(f"  [#1;D1] matches: {len(matches)} hydrogens")
    print(f"  Expected: {len(h_degrees)} hydrogens")
    
    if len(matches) < len(h_degrees):
        print(f"  ⚠️  MISSING: {len(h_degrees) - len(matches)} hydrogens NOT matched!")
else:
    print("ERROR: RDKit molecule not created!")
    print("Check for RDKit errors in molecule construction")
