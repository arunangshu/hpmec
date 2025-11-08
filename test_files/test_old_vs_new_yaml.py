"""
Test the OLD yaml (with extended SMARTS) to confirm it causes generic fallbacks.
"""

from calculator import load_molecule, infer_topology, assign_parameters, load_force_field

print("=" * 60)
print("TESTING OLD YAML (Extended SMARTS)")
print("=" * 60)
print()

try:
    mol = load_molecule('ethanol_sample.xyz')
    ff = load_force_field('ethanol.yaml')  # OLD yaml
    infer_topology(mol)
    assign_parameters(mol, ff)
    
    print("Atom typing results:")
    for atom in mol.atoms:
        fallback = "❌ FALLBACK" if 'generic_' in atom.atom_type else "✅"
        print(f"  {fallback} Atom {atom.index} ({atom.element}): {atom.atom_type}")
    
    generic_count = sum(1 for a in mol.atoms if 'generic_' in a.atom_type)
    
    print()
    if generic_count > 0:
        print(f"❌ PROBLEM CONFIRMED: {generic_count}/{len(mol.atoms)} atoms using generic fallbacks")
    else:
        print("✅ All atoms typed correctly")
        
except Exception as e:
    print(f"❌ Error: {e}")
    import traceback
    traceback.print_exc()

print()
print("=" * 60)
print("Now testing CORRECTED YAML (Base SMARTS)")
print("=" * 60)
print()

try:
    mol = load_molecule('ethanol_sample.xyz')
    ff = load_force_field('ethanol_corrected.yaml')  # CORRECTED yaml
    infer_topology(mol)
    assign_parameters(mol, ff)
    
    print("Atom typing results:")
    for atom in mol.atoms:
        fallback = "❌ FALLBACK" if 'generic_' in atom.atom_type else "✅"
        print(f"  {fallback} Atom {atom.index} ({atom.element}): {atom.atom_type}")
    
    generic_count = sum(1 for a in mol.atoms if 'generic_' in a.atom_type)
    
    print()
    if generic_count == 0:
        print(f"✅ SOLUTION CONFIRMED: {len(mol.atoms)}/{len(mol.atoms)} atoms typed correctly!")
    else:
        print(f"❌ Still has problems: {generic_count}/{len(mol.atoms)} atoms using fallbacks")
        
except Exception as e:
    print(f"❌ Error: {e}")
    import traceback
    traceback.print_exc()

print()
print("=" * 60)
