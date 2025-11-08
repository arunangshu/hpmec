from calculator import load_xyz, infer_topology
from rdkit import Chem

mol_obj = load_xyz('ethanol_sample.xyz')
infer_topology(mol_obj)
rdkit_mol = mol_obj.rdkit_mol

print("Testing SMARTS patterns:")
print(f"Molecule has {rdkit_mol.GetNumAtoms()} atoms")
print()

# Test simple patterns first
test_patterns = [
    ('[H]', 'any hydrogen'),
    ('[C]', 'any carbon'),  
    ('[O]', 'any oxygen'),
    ('[H;X1]', 'hydrogen degree 1'),
    ('[H;D1]', 'hydrogen degree 1 (D notation)'),
    ('[H;x1]', 'hydrogen degree 1 (lowercase x)'),
    ('[C;X4]', 'carbon degree 4'),
    ('[O;X2]', 'oxygen degree 2'),
    ('[H][C]', 'H-C bond'),
    ('[H][O]', 'H-O bond'),
    ('[H;D1][C]', 'H(D1)-C bond'),
    ('[H;D1][C;X4]', 'H(D1)-C(X4) bond'),
    ('[H;D1][O;X2]', 'H(D1)-O(X2) bond'),
    ('[#1;X1]', 'hydrogen by atomic num'),
    ('[#1;D1]', 'hydrogen (#1) with D1'),
    ('[#1][C]', 'H(#1)-C bond'),
    ('[#1;D1][C;X4]', 'H(#1,D1)-C(X4)'),
]

# Check hydrogen atom properties in detail
print("\nDetailed hydrogen properties:")
for atom in rdkit_mol.GetAtoms():
    if atom.GetSymbol() == 'H':
        idx = atom.GetIdx()
        print(f"  Atom {idx}: Symbol={atom.GetSymbol()}, "
              f"AtomicNum={atom.GetAtomicNum()}, "
              f"Degree={atom.GetDegree()}, "
              f"TotalDegree={atom.GetTotalDegree()}, "
              f"TotalValence={atom.GetTotalValence()}, "
              f"ExplicitValence={atom.GetExplicitValence()}, "
              f"ImplicitValence={atom.GetImplicitValence()}, "
              f"NumExplicitHs={atom.GetNumExplicitHs()}, "
              f"NumImplicitHs={atom.GetNumImplicitHs()}")

print("\nTesting patterns:")
for smarts, desc in test_patterns:
    patt = Chem.MolFromSmarts(smarts)
    matches = rdkit_mol.GetSubstructMatches(patt)
    print(f"{smarts:20s} ({desc:25s}): {len(matches)} matches")
    if len(matches) > 0 and len(matches) <= 10:
        print(f"  Matches: {matches}")
