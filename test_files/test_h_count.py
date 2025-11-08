from calculator import load_xyz, infer_topology

mol_obj = load_xyz('ethanol_sample.xyz')
infer_topology(mol_obj)
rdkit_mol = mol_obj.rdkit_mol

print("Checking hydrogen counts on carbon atoms:")
print()
for atom in rdkit_mol.GetAtoms():
    if atom.GetSymbol() == 'C':
        idx = atom.GetIdx()
        neighbors = atom.GetNeighbors()
        h_neighbors = [n for n in neighbors if n.GetSymbol() == 'H']
        
        print(f"Atom {idx} (C):")
        print(f"  GetTotalNumHs() = {atom.GetTotalNumHs()}")
        print(f"  GetNumExplicitHs() = {atom.GetNumExplicitHs()}")
        print(f"  GetNumImplicitHs() = {atom.GetNumImplicitHs()}")
        print(f"  Actual H neighbors = {len(h_neighbors)}")
        print(f"  All neighbors: {[n.GetSymbol() for n in neighbors]}")
        print()
