"""
Show which atoms in t3 are not being typed and why
"""
from calculator import load_xyz, load_force_field, assign_parameters, infer_topology
from rdkit import Chem

# Load
mol = load_xyz('t3.xyz')
infer_topology(mol)
ff = load_force_field('t3_forcefield.yaml')

# Assign
assign_parameters(mol, ff)

print("=" * 80)
print("T3 ATOM TYPING ANALYSIS")
print("=" * 80)
print()

# Show which atoms got typed
untyped = []
for atom in mol.atoms:
    if atom.atom_type and atom.atom_type.startswith('generic_'):
        untyped.append(atom)

print(f"Untyped atoms: {len(untyped)}/{len(mol.atoms)}")
print()

if untyped:
    print("UNTYPED ATOMS (using generic fallback):")
    print()
    
    rdkit_mol = mol.rdkit_mol
    if rdkit_mol:
        for atom in untyped:
            idx = atom.index
            rdk_atom = rdkit_mol.GetAtomWithIdx(idx)
            degree = rdk_atom.GetDegree()
            neighbors = [n.GetSymbol() for n in rdk_atom.GetNeighbors()]
            h_count = sum(1 for n in neighbors if n == 'H')
            
            print(f"Atom {idx:2d}: {atom.element} → {atom.atom_type}")
            print(f"  Degree: {degree}, H count: {h_count}")
            print(f"  Neighbors: {neighbors}")
            print(f"  Expected SMARTS: [{atom.element};D{degree};H{h_count}]" if h_count > 0 else f"[{atom.element};D{degree}]")
            print()

print("=" * 80)
print("SOLUTION:")
print("=" * 80)
print("Add these SMARTS patterns to your YAML:")
print()

# Generate missing patterns
if rdkit_mol:
    missing_patterns = set()
    for atom in untyped:
        idx = atom.index
        rdk_atom = rdkit_mol.GetAtomWithIdx(idx)
        degree = rdk_atom.GetDegree()
        neighbors = [n.GetSymbol() for n in rdk_atom.GetNeighbors()]
        h_count = sum(1 for n in neighbors if n == 'H')
        
        symbol = '#1' if atom.element == 'H' else atom.element
        if h_count > 0 and atom.element != 'H':
            pattern = f"[{symbol};D{degree};H{h_count}]"
        else:
            pattern = f"[{symbol};D{degree}]"
        
        missing_patterns.add((pattern, atom.element, degree, h_count))
    
    for pattern, elem, deg, hcount in sorted(missing_patterns):
        type_name = f"{elem}_deg{deg}_h{hcount}"
        print(f"- smarts: '{pattern}'")
        print(f"  type_name: {type_name}")
        print(f"  charge: 0.0")
        print(f"  sigma: 0.35")
        print(f"  epsilon: 0.276")
        print()
