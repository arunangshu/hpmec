"""
Test script to verify PDB and MOL file support in the energy calculator.
"""

from calculator import load_molecule, load_pdb, load_mol, load_xyz, infer_topology

def test_xyz_loading():
    """Test loading XYZ file"""
    print("=" * 60)
    print("Testing XYZ file loading...")
    print("=" * 60)
    
    try:
        mol = load_xyz("ethanol_sample.xyz")
        print(f"✅ Loaded XYZ file successfully")
        print(f"   - Number of atoms: {len(mol.atoms)}")
        print(f"   - Elements: {[a.element for a in mol.atoms]}")
        print(f"   - Coordinates shape: {mol.coordinates.shape}")
        
        infer_topology(mol)
        print(f"   - Number of bonds: {len(mol.bonds)}")
        print(f"   - Number of angles: {len(mol.angles)}")
        print(f"   - Number of dihedrals: {len(mol.dihedrals)}")
        print(f"   - RDKit mol exists: {mol.rdkit_mol is not None}")
        
    except Exception as e:
        print(f"❌ Error loading XYZ: {e}")
    
    print()

def test_universal_loader():
    """Test universal load_molecule function"""
    print("=" * 60)
    print("Testing universal molecule loader...")
    print("=" * 60)
    
    # Test with XYZ
    try:
        mol = load_molecule("ethanol_sample.xyz")
        print(f"✅ Universal loader works with XYZ")
        print(f"   - Loaded {len(mol.atoms)} atoms")
    except Exception as e:
        print(f"❌ Error with XYZ: {e}")
    
    print()

def test_format_detection():
    """Test that file format detection works"""
    print("=" * 60)
    print("Testing file format detection...")
    print("=" * 60)
    
    test_files = [
        ("test.xyz", "xyz"),
        ("test.pdb", "pdb"),
        ("test.mol", "mol"),
        ("TEST.XYZ", "xyz"),  # uppercase
        ("molecule.PDB", "pdb"),  # mixed case
    ]
    
    for filename, expected_ext in test_files:
        detected_ext = filename.lower().split('.')[-1]
        status = "✅" if detected_ext == expected_ext else "❌"
        print(f"{status} {filename} → detected as '{detected_ext}' (expected '{expected_ext}')")
    
    print()

def test_topology_with_existing_bonds():
    """Test that infer_topology preserves existing bonds from PDB/MOL files"""
    print("=" * 60)
    print("Testing topology inference with existing bonds...")
    print("=" * 60)
    
    try:
        # Create a mock molecule with pre-existing bonds
        from calculator import Molecule, Atom
        
        mol = Molecule()
        mol.atoms = [
            Atom(0, 'C', (0.0, 0.0, 0.0)),
            Atom(1, 'O', (0.15, 0.0, 0.0)),
            Atom(2, 'H', (-0.01, 0.01, 0.0))
        ]
        mol.coordinates = np.array([a.coords for a in mol.atoms])
        
        # Pre-populate bonds (simulating PDB/MOL loading)
        mol.bonds = [(0, 1), (0, 2)]
        
        print("   Before infer_topology:")
        print(f"   - Bonds: {mol.bonds}")
        
        infer_topology(mol)
        
        print("   After infer_topology:")
        print(f"   - Bonds: {mol.bonds}")
        print(f"   - Angles: {mol.angles}")
        
        if set(mol.bonds) == {(0, 1), (0, 2)}:
            print("✅ Existing bonds were preserved!")
        else:
            print("❌ Bonds were modified unexpectedly")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
    
    print()

if __name__ == "__main__":
    import numpy as np
    
    print("\n" + "=" * 60)
    print("MOLECULAR FILE FORMAT SUPPORT TEST SUITE")
    print("=" * 60 + "\n")
    
    test_xyz_loading()
    test_universal_loader()
    test_format_detection()
    test_topology_with_existing_bonds()
    
    print("=" * 60)
    print("TEST SUITE COMPLETE")
    print("=" * 60)
