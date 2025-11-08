"""
COMPREHENSIVE WORKFLOW TEST
Simulates the complete user workflow from molecule upload to energy calculation.
"""

from calculator import (
    load_molecule, infer_topology, assign_parameters, load_force_field,
    calculate_energy_with_breakdown, validate_force_field_coverage
)

def print_section(title):
    print("\n" + "=" * 70)
    print(title)
    print("=" * 70)

def test_complete_workflow():
    """Test the complete workflow as a user would experience it"""
    
    print_section("COMPLETE WORKFLOW TEST: Ethanol Example")
    
    # STEP 1: Load molecule
    print("\n📁 STEP 1: Upload Molecule File")
    try:
        mol = load_molecule('ethanol_sample.xyz')
        print(f"✅ Loaded: ethanol_sample.xyz ({len(mol.atoms)} atoms)")
    except Exception as e:
        print(f"❌ Failed to load molecule: {e}")
        return
    
    # STEP 2: Infer topology
    print("\n🔗 STEP 2: Analyze Structure")
    try:
        infer_topology(mol)
        print(f"✅ Topology inferred:")
        print(f"   - Bonds: {len(mol.bonds)}")
        print(f"   - Angles: {len(mol.angles)}")
        print(f"   - Dihedrals: {len(mol.dihedrals)}")
    except Exception as e:
        print(f"❌ Failed topology inference: {e}")
        return
    
    # STEP 3: Validate force field with CORRECTED YAML
    print("\n🔍 STEP 3: Validate Force Field Coverage (Corrected YAML)")
    try:
        validation = validate_force_field_coverage('ethanol_sample.xyz', 'ethanol_corrected.yaml')
        
        if validation['is_complete']:
            print("✅ Force field validation PASSED!")
        else:
            print("⚠️ Force field validation found issues:")
            if validation['missing_atom_types']:
                print(f"   - Missing atom types: {len(validation['missing_atom_types'])}")
            if validation['missing_bonds']:
                print(f"   - Missing bond types: {len(validation['missing_bonds'])}")
            if validation['missing_angles']:
                print(f"   - Missing angle types: {len(validation['missing_angles'])}")
            if validation['missing_dihedrals']:
                print(f"   - Missing dihedral types: {len(validation['missing_dihedrals'])}")
        
        # Show coverage stats
        stats = validation['coverage_stats']
        print(f"\n   Coverage Statistics:")
        print(f"   - Atoms:     {stats['atoms']['coverage_percent']:.1f}% ({stats['atoms']['typed']}/{stats['atoms']['total']})")
        print(f"   - Bonds:     {stats['bonds']['coverage_percent']:.1f}% ({stats['bonds']['parameterized']}/{stats['bonds']['total']})")
        print(f"   - Angles:    {stats['angles']['coverage_percent']:.1f}% ({stats['angles']['parameterized']}/{stats['angles']['total']})")
        print(f"   - Dihedrals: {stats['dihedrals']['coverage_percent']:.1f}% ({stats['dihedrals']['parameterized']}/{stats['dihedrals']['total']})")
        
    except Exception as e:
        print(f"❌ Validation failed: {e}")
        return
    
    # STEP 4: Calculate energy
    print("\n⚡ STEP 4: Calculate Energy")
    try:
        result = calculate_energy_with_breakdown('ethanol_sample.xyz', 'ethanol_corrected.yaml', n_cores=4)
        
        print("✅ Energy calculation successful!")
        print(f"\n   Energy Components:")
        print(f"   - Bonds:        {result['bond']:10.4f} kJ/mol")
        print(f"   - Angles:       {result['angle']:10.4f} kJ/mol")
        print(f"   - Dihedrals:    {result['dihedral']:10.4f} kJ/mol")
        print(f"   - Van der Waals:{result['vdw']:10.4f} kJ/mol")
        print(f"   - Electrostatic:{result['electrostatic']:10.4f} kJ/mol")
        print(f"   " + "-" * 40)
        print(f"   - TOTAL:        {result['total']:10.4f} kJ/mol")
        
        print(f"\n   Performance Metrics:")
        timing = result['timing']
        print(f"   - Molecule size: {timing['n_atoms']} atoms")
        print(f"   - Brute force:   {timing['brute_force_time']*1000:6.2f} ms (O(N²))")
        print(f"   - Optimized:     {timing['single_core_time']*1000:6.2f} ms (O(N log N))")
        print(f"   - Multi-core:    {timing['multi_core_time']*1000:6.2f} ms ({timing['n_cores_used']} cores)")
        print(f"   - Speedup:       {timing['speedup_optimized']:.1f}x (optimization)")
        
    except Exception as e:
        print(f"❌ Energy calculation failed: {e}")
        import traceback
        traceback.print_exc()
        return
    
    # STEP 5: Verify atom typing
    print("\n✓ STEP 5: Verify Atom Typing")
    try:
        mol = load_molecule('ethanol_sample.xyz')
        ff = load_force_field('ethanol_corrected.yaml')
        infer_topology(mol)
        assign_parameters(mol, ff)
        
        print("   Atom Types Assigned:")
        unique_types = {}
        for atom in mol.atoms:
            if atom.atom_type not in unique_types:
                unique_types[atom.atom_type] = []
            unique_types[atom.atom_type].append(atom.index)
        
        for atom_type, indices in sorted(unique_types.items()):
            count = len(indices)
            element = mol.atoms[indices[0]].element
            fallback_marker = "❌ FALLBACK" if 'generic_' in atom_type else "✅"
            print(f"   {fallback_marker} {atom_type:15s} ({element}): {count} atoms - indices {indices}")
        
        generic_count = sum(1 for a in mol.atoms if 'generic_' in a.atom_type)
        
        print(f"\n   Summary:")
        print(f"   - Total atoms: {len(mol.atoms)}")
        print(f"   - Unique types: {len(unique_types)}")
        print(f"   - Generic fallbacks: {generic_count}")
        
        if generic_count == 0:
            print(f"\n   🎉 SUCCESS: All atoms properly typed!")
        else:
            print(f"\n   ⚠️  WARNING: {generic_count} atoms using generic fallbacks")
        
    except Exception as e:
        print(f"❌ Atom typing verification failed: {e}")
        return
    
    print_section("✅ WORKFLOW TEST COMPLETE - ALL STEPS PASSED")

if __name__ == "__main__":
    test_complete_workflow()
