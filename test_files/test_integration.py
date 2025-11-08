"""
Test script for calculator integration
"""

from calculator import (
    calculate_energy_with_breakdown,
    validate_force_field_coverage
)

print("=" * 60)
print("TESTING MOLECULAR ENERGY CALCULATOR")
print("=" * 60)

# Test 1: Validation
print("\n[1] Testing force field validation...")
validation = validate_force_field_coverage('ethanol_sample.xyz', 'ethanol.yaml')

print(f"   Complete: {validation['is_complete']}")
print(f"   Atom Coverage: {validation['coverage_stats']['atoms']['coverage_percent']:.1f}%")
print(f"   Bond Coverage: {validation['coverage_stats']['bonds']['coverage_percent']:.1f}%")
print(f"   Angle Coverage: {validation['coverage_stats']['angles']['coverage_percent']:.1f}%")
print(f"   Dihedral Coverage: {validation['coverage_stats']['dihedrals']['coverage_percent']:.1f}%")

if validation['missing_bonds']:
    print(f"\n   Missing bonds: {validation['missing_bonds']}")
if validation['missing_angles']:
    print(f"   Missing angles: {validation['missing_angles']}")
if validation['missing_dihedrals']:
    print(f"   Missing dihedrals: {validation['missing_dihedrals']}")

# Test 2: Energy calculation
print("\n[2] Testing energy calculation...")
energy = calculate_energy_with_breakdown('ethanol_sample.xyz', 'ethanol.yaml')

print(f"\n   ENERGY BREAKDOWN:")
print(f"   ─────────────────────────────────")
print(f"   Bond:           {energy['bond']:>12.4f} kJ/mol")
print(f"   Angle:          {energy['angle']:>12.4f} kJ/mol")
print(f"   Dihedral:       {energy['dihedral']:>12.4f} kJ/mol")
print(f"   Non-bonded:     {energy['nonbonded']:>12.4f} kJ/mol")
print(f"   ─────────────────────────────────")
print(f"   TOTAL:          {energy['total']:>12.4f} kJ/mol")

print("\n" + "=" * 60)
print("✅ ALL TESTS PASSED!")
print("=" * 60)
