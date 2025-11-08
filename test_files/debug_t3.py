"""
Debug script to see which atoms are not getting typed in t3 molecule
"""
from calculator import load_xyz, load_force_field, assign_parameters, infer_topology

# Load t3 molecule
mol = load_xyz('t3.xyz')
infer_topology(mol)

# Load force field
ff = load_force_field('t3.yaml')

# Assign types
assign_parameters(mol, ff)

print("=" * 70)
print("ATOM TYPE ASSIGNMENTS FOR T3")
print("=" * 70)
print()

# Group by type
from collections import defaultdict
by_type = defaultdict(list)
for atom in mol.atoms:
    by_type[atom.atom_type].append(f"{atom.index}:{atom.element}")

for atom_type in sorted(by_type.keys()):
    atoms_str = ", ".join(by_type[atom_type])
    print(f"{atom_type:20s}: {atoms_str}")

print()
print("=" * 70)
print("UNTYPED ATOMS (using generic fallback):")
print("=" * 70)
for atom in mol.atoms:
    if atom.atom_type and atom.atom_type.startswith('generic_'):
        neighbors = [mol.atoms[n].element for n in mol.bonds.get(atom.index, [])]
        print(f"  Atom {atom.index:2d} ({atom.element}): {atom.atom_type:15s} - neighbors: {neighbors}")
