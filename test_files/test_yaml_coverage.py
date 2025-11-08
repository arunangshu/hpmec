"""
Test script to verify YAML Builder generates complete coverage
Run this after regenerating ethanol.yaml with the YAML Builder
"""
from calculator import load_xyz, load_force_field, validate_force_field_coverage

print("=" * 60)
print("TESTING YAML BUILDER OUTPUT")
print("=" * 60)
print()

# Load molecule
mol_obj = load_xyz('ethanol_sample.xyz')

# Load YAML (user should regenerate this first)
ff_params = load_force_field('ethanol.yaml')

print("[1] Validating force field coverage...")
result = validate_force_field_coverage(mol_obj, ff_params)

print()
print("Coverage Results:")
print(f"  Atoms: {result['coverage']['atoms']['typed']}/{result['coverage']['atoms']['total']} ({result['coverage']['atoms']['coverage_percent']:.1f}%)")
print(f"  Bonds: {result['coverage']['bonds']['parameterized']}/{result['coverage']['bonds']['total']} ({result['coverage']['bonds']['coverage_percent']:.1f}%)")
print(f"  Angles: {result['coverage']['angles']['parameterized']}/{result['coverage']['angles']['total']} ({result['coverage']['angles']['coverage_percent']:.1f}%)")
print(f"  Dihedrals: {result['coverage']['dihedrals']['parameterized']}/{result['coverage']['dihedrals']['total']} ({result['coverage']['dihedrals']['coverage_percent']:.1f}%)")
print()

if result['complete']:
    print("✅ SUCCESS: 100% coverage achieved!")
    print("   No generic fallback types used.")
else:
    print("⚠️  INCOMPLETE COVERAGE")
    if result['missing']['atom_types']:
        print(f"   Missing atom types: {result['missing']['atom_types']}")
    if result['missing']['bond_params']:
        print(f"   Missing bond types: {result['missing']['bond_params'][:5]}")  # Show first 5
    if result['missing']['angle_params']:
        print(f"   Missing angle types: {result['missing']['angle_params'][:5]}")
    if result['missing']['dihedral_params']:
        print(f"   Missing dihedral types: {result['missing']['dihedral_params'][:5]}")

print()
print("=" * 60)
print("Atom Type Assignments:")
print("=" * 60)
for atom in mol_obj.atoms:
    print(f"  Atom {atom.index}: {atom.element:2s} → {atom.atom_type}")

print()
print("NEXT STEPS:")
print("1. Open Streamlit GUI: streamlit run app.py")
print("2. Upload ethanol_sample.xyz")
print("3. Click 'Generate YAML Builder Parameters'")
print("4. Download and save as ethanol.yaml")
print("5. Run this test again to verify 100% coverage")
print("=" * 60)
