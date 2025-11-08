"""
Test if corrected YAML with base SMARTS works without generic fallbacks.
"""

from calculator import calculate_energy_with_breakdown

print("=" * 60)
print("TESTING CORRECTED YAML WITH BASE SMARTS")
print("=" * 60)
print()

try:
    result = calculate_energy_with_breakdown('ethanol_sample.xyz', 'ethanol_corrected.yaml', n_cores=4)
    
    print("✅ Energy calculation successful!")
    print()
    print("Energy Results:")
    print(f"  Total Energy: {result['total']:.4f} kJ/mol")
    print(f"  Bonded:       {result['bonded']:.4f} kJ/mol")
    print(f"  Non-bonded:   {result['nonbonded']:.4f} kJ/mol")
    print()
    print("Performance:")
    print(f"  Molecule size: {result['timing']['n_atoms']} atoms")
    print(f"  Brute force O(N²): {result['timing']['brute_force_time']*1000:.2f} ms")
    print(f"  Optimized (k-d tree): {result['timing']['single_core_time']*1000:.2f} ms")
    print(f"  Speedup: {result['timing']['speedup_optimized']:.1f}x")
    print()
    
    # Check for generic fallbacks
    print("Checking for generic fallback types...")
    from calculator import load_molecule, infer_topology, assign_parameters, load_force_field
    
    mol = load_molecule('ethanol_sample.xyz')
    ff = load_force_field('ethanol_corrected.yaml')
    infer_topology(mol)
    assign_parameters(mol, ff)
    
    generic_count = 0
    for atom in mol.atoms:
        if 'generic_' in atom.atom_type:
            print(f"  ❌ Atom {atom.index} ({atom.element}): {atom.atom_type}")
            generic_count += 1
    
    if generic_count == 0:
        print("  ✅ NO GENERIC FALLBACKS - All atoms typed correctly!")
    else:
        print(f"  ❌ Found {generic_count} atoms using generic fallbacks")
    
except Exception as e:
    print(f"❌ Error: {e}")
    import traceback
    traceback.print_exc()

print()
print("=" * 60)
print("TEST COMPLETE")
print("=" * 60)
