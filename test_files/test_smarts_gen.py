from calculator import load_xyz, infer_topology

mol_obj = load_xyz('ethanol_sample.xyz')
infer_topology(mol_obj)
rdkit_mol = mol_obj.rdkit_mol

print("Testing SMARTS generation with explicit H count:")
print()

for atom in rdkit_mol.GetAtoms():
    if atom.GetSymbol() in ['C', 'O']:
        idx = atom.GetIdx()
        symbol = atom.GetSymbol()
        degree = atom.GetDegree()
        neighbors = atom.GetNeighbors()
        
        # Count explicit H neighbors (the fix)
        explicit_h_count = sum(1 for n in neighbors if n.GetSymbol() == 'H')
        
        # Find most electronegative non-H neighbor for extended SMARTS
        ELECTRONEGATIVITY = {
            'F': 4.0, 'O': 3.5, 'N': 3.0, 'Cl': 3.0, 'Br': 2.8,
            'S': 2.5, 'C': 2.5, 'P': 2.1, 'H': 2.1, 'Si': 1.8
        }
        
        heavy_neighbors = [n for n in neighbors if n.GetSymbol() != 'H']
        
        if heavy_neighbors:
            sig_neighbor = max(heavy_neighbors, 
                             key=lambda n: ELECTRONEGATIVITY.get(n.GetSymbol(), 2.0))
        else:
            sig_neighbor = None
        
        # Build SMARTS
        smarts_parts = [symbol, f"D{degree}"]
        if explicit_h_count > 0:
            smarts_parts.append(f"H{explicit_h_count}")
        
        base_smarts = f"[{';'.join(smarts_parts)}]"
        
        if sig_neighbor:
            n_degree = sig_neighbor.GetDegree()
            n_symbol = sig_neighbor.GetSymbol()
            n_neighbors = sig_neighbor.GetNeighbors()
            n_h_count = sum(1 for nn in n_neighbors if nn.GetSymbol() == 'H')
            
            n_parts = [n_symbol, f"D{n_degree}"]
            if n_h_count > 0 and n_symbol != 'H':
                n_parts.append(f"H{n_h_count}")
            
            neighbor_smarts = f"[{';'.join(n_parts)}]"
            extended_smarts = f"{base_smarts}{neighbor_smarts}"
        else:
            extended_smarts = base_smarts
        
        print(f"Atom {idx} ({symbol}):")
        print(f"  Explicit H count: {explicit_h_count}")
        print(f"  Neighbors: {[n.GetSymbol() for n in neighbors]}")
        print(f"  Base SMARTS: {base_smarts}")
        print(f"  Extended SMARTS: {extended_smarts}")
        print()
