"""Test environment-based grouping"""
from calculator import load_xyz, infer_topology

mol = load_xyz('ethanol_sample.xyz')
infer_topology(mol)
rdkit_mol = mol.rdkit_mol

print("=" * 70)
print("ATOM GROUPING BY ELEMENT + NEIGHBORS")
print("=" * 70)
print()

from collections import defaultdict
groups = defaultdict(list)

for atom in rdkit_mol.GetAtoms():
    idx = atom.GetIdx()
    symbol = atom.GetSymbol()
    neighbors = atom.GetNeighbors()
    neighbor_symbols = sorted([n.GetSymbol() for n in neighbors])
    
    # Create environment key
    neighbor_key = "_".join(neighbor_symbols) if neighbor_symbols else "isolated"
    environment_key = f"{symbol}_{neighbor_key}"
    
    groups[environment_key].append(idx)
    
    print(f"Atom {idx:2d} ({symbol}): neighbors={neighbor_symbols} → key='{environment_key}'")

print()
print("=" * 70)
print("UNIQUE ENVIRONMENTS:")
print("=" * 70)
for env_key, atom_indices in sorted(groups.items()):
    print(f"{env_key:30s}: {len(atom_indices)} atoms → {atom_indices}")

print()
print(f"Total unique environments: {len(groups)}")
print(f"Expected for ethanol: 6 (CH3, CH2, OH, 3 types of H)")
