from calculator import load_xyz, infer_topology
from rdkit import Chem

mol_obj = load_xyz('ethanol_sample.xyz')
infer_topology(mol_obj)
rdkit_mol = mol_obj.rdkit_mol

print("Testing new SMARTS patterns with explicit H count:")
print()

test_patterns = [
    ('[C;D4;H3][C;D4;H2]', 'CH3 carbon'),
    ('[C;D4;H2][C;D4;H3]', 'CH2 carbon (reverse)'),
    ('[C;D4;H2][O;D2]', 'CH2 carbon (with O neighbor)'),
    ('[O;D2;H1][C;D4]', 'OH oxygen'),
    ('[#1;D1][C;D4;H3]', 'H on CH3'),
    ('[#1;D1][C;D4;H2]', 'H on CH2'),
    ('[#1;D1][O;D2;H1]', 'H on OH'),
]

for smarts, desc in test_patterns:
    patt = Chem.MolFromSmarts(smarts)
    if patt:
        matches = rdkit_mol.GetSubstructMatches(patt)
        print(f"{smarts:25s} ({desc:25s}): {len(matches)} matches")
        if len(matches) > 0 and len(matches) <= 10:
            print(f"  Matches: {matches}")
    else:
        print(f"{smarts:25s} ({desc:25s}): INVALID PATTERN")
    print()
