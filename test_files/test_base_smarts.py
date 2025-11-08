"""
Test script to verify SMARTS pattern generation produces BASE patterns only.
"""

from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
import yaml

def test_smarts_generation():
    """Test that we generate base SMARTS correctly for ethanol"""
    
    # Load ethanol
    mol = Chem.MolFromXYZFile("ethanol_sample.xyz")
    if mol is None:
        print("❌ Failed to load ethanol")
        return
    
    # Determine bonds
    rdDetermineBonds.DetermineBonds(mol, charge=0)
    
    print("=" * 60)
    print("ETHANOL SMARTS PATTERN GENERATION TEST")
    print("=" * 60)
    print(f"Number of atoms: {mol.GetNumAtoms()}")
    print()
    
    # Generate SMARTS for each atom
    smarts_data = []
    
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        symbol = atom.GetSymbol()
        degree = atom.GetDegree()
        
        # Count explicit hydrogens
        neighbors = atom.GetNeighbors()
        explicit_h_count = sum(1 for n in neighbors if n.GetSymbol() == 'H')
        
        # Build base SMARTS (single atom pattern)
        smarts_parts = []
        
        if symbol == 'H':
            smarts_parts.append('#1')  # Use atomic number for H
        else:
            smarts_parts.append(symbol)
        
        smarts_parts.append(f'D{degree}')
        smarts_parts.append(f'H{explicit_h_count}')
        
        base_smarts = f"[{';'.join(smarts_parts)}]"
        
        smarts_data.append({
            'atom_idx': idx,
            'element': symbol,
            'degree': degree,
            'h_count': explicit_h_count,
            'base_smarts': base_smarts
        })
        
        print(f"Atom {idx:2d} ({symbol:2s}): degree={degree}, H={explicit_h_count} → {base_smarts}")
    
    print()
    print("=" * 60)
    print("UNIQUE PATTERNS (grouped by base SMARTS):")
    print("=" * 60)
    
    # Group by base SMARTS
    groups = {}
    for item in smarts_data:
        key = item['base_smarts']
        if key not in groups:
            groups[key] = {
                'smarts': key,
                'element': item['element'],
                'count': 0,
                'indices': []
            }
        groups[key]['count'] += 1
        groups[key]['indices'].append(item['atom_idx'])
    
    for i, (pattern, data) in enumerate(groups.items(), 1):
        print(f"\n{i}. Pattern: {pattern}")
        print(f"   Element: {data['element']}")
        print(f"   Count: {data['count']} atoms")
        print(f"   Atom indices: {data['indices']}")
    
    print()
    print("=" * 60)
    print("VERIFICATION:")
    print("=" * 60)
    
    # Verify patterns match
    all_typed = True
    for pattern_data in groups.values():
        pattern = pattern_data['smarts']
        
        try:
            patt = Chem.MolFromSmarts(pattern)
            if patt is None:
                print(f"❌ Invalid SMARTS: {pattern}")
                continue
            
            matches = mol.GetSubstructMatches(patt)
            matched_atoms = set()
            for match in matches:
                if len(match) > 0:
                    matched_atoms.add(match[0])  # Single-atom pattern
            
            expected = set(pattern_data['indices'])
            if matched_atoms == expected:
                print(f"✅ {pattern} matches {len(matched_atoms)} atoms correctly")
            else:
                print(f"❌ {pattern} mismatch!")
                print(f"   Expected: {expected}")
                print(f"   Got: {matched_atoms}")
                all_typed = False
                
        except Exception as e:
            print(f"❌ Error matching {pattern}: {e}")
            all_typed = False
    
    print()
    if all_typed:
        print("✅ ALL PATTERNS VERIFIED - All atoms can be typed!")
    else:
        print("❌ PATTERN ERRORS DETECTED")
    
    return groups

def generate_corrected_yaml(groups):
    """Generate corrected ethanol.yaml with base SMARTS"""
    
    print()
    print("=" * 60)
    print("GENERATING CORRECTED YAML:")
    print("=" * 60)
    
    atom_types = []
    for i, (pattern, data) in enumerate(groups.items(), 1):
        atom_types.append({
            'smarts': pattern,  # BASE SMARTS
            'type_name': f"{data['element']}_type_{i}",
            'charge': 0.0,
            'sigma': 0.35,
            'epsilon': 0.276
        })
    
    yaml_content = {
        'atom_types': atom_types,
        'bond_types': {
            'C_type_1-C_type_2': [224262.4, 0.1529],
            'C_type_2-O_type_3': [224262.4, 0.1529],
            'C_type_1-H_type_4': [224262.4, 0.1529],
            'C_type_2-H_type_4': [224262.4, 0.1529],
            'H_type_4-O_type_3': [224262.4, 0.1529]
        },
        'angle_types': {
            'C_type_2-C_type_1-H_type_4': [418.4, 1.911],
            'H_type_4-C_type_1-H_type_4': [418.4, 1.911],
            'C_type_1-C_type_2-O_type_3': [418.4, 1.911],
            'C_type_1-C_type_2-H_type_4': [418.4, 1.911],
            'H_type_4-C_type_2-O_type_3': [418.4, 1.911],
            'H_type_4-C_type_2-H_type_4': [418.4, 1.911],
            'C_type_2-O_type_3-H_type_4': [418.4, 1.911]
        },
        'dihedral_types': {
            'O_type_3-C_type_2-C_type_1-H_type_4': [0.0, 0.0, 0.0, 0.0],
            'H_type_4-C_type_2-C_type_1-H_type_4': [0.0, 0.0, 0.0, 0.0],
            'H_type_4-O_type_3-C_type_2-C_type_1': [0.0, 0.0, 0.0, 0.0],
            'H_type_4-O_type_3-C_type_2-H_type_4': [0.0, 0.0, 0.0, 0.0]
        }
    }
    
    with open('ethanol_corrected.yaml', 'w') as f:
        f.write("# Generated Force Field Parameters (CORRECTED)\n")
        f.write("# Created for molecule: ethanol_sample.xyz\n")
        f.write("# Number of atoms: 9\n")
        f.write("# UNITS: kJ/mol, nm, radians, elementary charge (e)\n")
        f.write("# BASE SMARTS patterns (single-atom)\n")
        f.write("# -----------------------------------------------------------------\n\n")
        yaml.dump(yaml_content, f, default_flow_style=False, sort_keys=False)
    
    print("✅ Saved to: ethanol_corrected.yaml")
    print()
    print("YAML Preview:")
    print(yaml.dump({'atom_types': atom_types}, default_flow_style=False))

if __name__ == "__main__":
    groups = test_smarts_generation()
    if groups:
        generate_corrected_yaml(groups)
