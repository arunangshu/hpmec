"""
Test scenario: Using ethanol.yaml with T3 molecule (mismatch).
This is likely the user's actual scenario.
"""

from calculator import load_molecule, infer_topology, assign_parameters, load_force_field

print("=" * 60)
print("SCENARIO: Using ETHANOL force field on T3 molecule")
print("(This simulates using the wrong YAML for a molecule)")
print("=" * 60)
print()

try:
    mol = load_molecule('t3.xyz')
    ff = load_force_field('ethanol.yaml')  # WRONG FF for this molecule!
    infer_topology(mol)
    print(f"Loaded T3: {len(mol.atoms)} atoms")
    print()
    
    assign_parameters(mol, ff)
    
    print("Atom typing results:")
    generic_count = 0
    for i, atom in enumerate(mol.atoms):
        if 'generic_' in atom.atom_type:
            print(f"  ❌ Atom {atom.index} ({atom.element}): {atom.atom_type}")
            generic_count += 1
    
    print()
    print(f"Result: {generic_count}/{len(mol.atoms)} atoms using generic fallbacks")
    print()
    print("💡 This is EXPECTED behavior:")
    print("   - ethanol.yaml only has patterns for ethanol atoms")
    print("   - T3 has different atom types (aromatic C, iodine, etc.)")
    print("   - Unmatched atoms get generic fallbacks")
    print()
    print("✅ SOLUTION: Generate force field YAML for each specific molecule")
        
except FileNotFoundError as e:
    print(f"⚠️ T3 file not found - skipping this test")
    print(f"   (This test shows what happens when using wrong YAML)")
except Exception as e:
    print(f"❌ Error: {e}")
    import traceback
    traceback.print_exc()

print()
print("=" * 60)
print("CORRECT SCENARIO: Using ethanol YAML with ethanol molecule")
print("=" * 60)
print()

try:
    mol = load_molecule('ethanol_sample.xyz')
    ff = load_force_field('ethanol.yaml')
    infer_topology(mol)
    assign_parameters(mol, ff)
    
    generic_count = sum(1 for a in mol.atoms if 'generic_' in a.atom_type)
    
    if generic_count == 0:
        print(f"✅ CORRECT: {len(mol.atoms)}/{len(mol.atoms)} atoms typed successfully")
        print("   All ethanol atoms match ethanol force field patterns")
    else:
        print(f"❌ PROBLEM: {generic_count} atoms using fallbacks")
        
except Exception as e:
    print(f"❌ Error: {e}")
    import traceback
    traceback.print_exc()

print()
print("=" * 60)
