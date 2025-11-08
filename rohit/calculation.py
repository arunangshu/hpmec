import numpy as np
from itertools import combinations

# ========================
# Utility: Covalent Radii
# ========================
COVALENT_RADII = {
    'H': 0.37, 'C': 0.77, 'O': 0.73, 'N': 0.75, 'S': 1.02
}

# ========================
# Step 1: Read XYZ file
# ========================
def read_xyz(filename):
    atoms, coords = [], []
    with open(filename, 'r') as f:
        lines = f.readlines()[2:]  # skip atom count + comment line
        for line in lines:
            parts = line.split()
            atoms.append(parts[0])
            coords.append(list(map(float, parts[1:4])))
    return np.array(atoms), np.array(coords)

# ========================
# Step 2: Infer Bonds
# ========================
def infer_bonds(atoms, coords, scale=1.2):
    bonds = []
    n = len(atoms)
    for i in range(n):
        for j in range(i+1, n):
            r = np.linalg.norm(coords[i] - coords[j])
            cutoff = scale * (COVALENT_RADII.get(atoms[i], 0.7) +
                              COVALENT_RADII.get(atoms[j], 0.7))
            if r < cutoff:
                bonds.append((i, j))
    return bonds

# ========================
# Step 3: Generate Angles
# ========================
def generate_angles(bonds):
    angles = []
    for i, j in bonds:
        for k, l in bonds:
            if j == k and i != l:
                angles.append((i, j, l))
    return list(set(angles))  # remove duplicates

# ========================
# Step 4: Generate Dihedrals
# ========================
def generate_dihedrals(bonds):
    dihedrals = set()
    for (i, j), (k, l) in combinations(bonds, 2):
        if j == k:
            for m, n in bonds:
                if l == m and n != i:
                    dihedrals.add((i, j, l, n))
    return list(dihedrals)
 #=============================
 #get non bonded atoms
 #==============================
def generate_nonbondedAtoms(atoms,bonds):
  nonbonded_pairs = []
  n = len(atoms)
  bonded_set = set(bonds)
  for i in range(n):
      for j in range(i+1, n):
          if (i,j) in bonded_set or (j,i) in bonded_set:
              continue
          nonbonded_pairs.append((i,j))


# ========================
# Step 5: Force Field (Dummy) to be written
# ========================
FF = {
    "bond_k": 300.0,        # kJ/mol·Å²
    "r0": 1.0,              # Å
    "angle_k": 40.0,        # kJ/mol·rad²
    "theta0": np.deg2rad(109.5),  # radians
    "dihedral_Vn": 2.0,     # kJ/mol
    "n": 3,
    "gamma": 0.0,
    "epsilon": 0.2,         # kJ/mol
    "sigma": 3.5,           # Å
}

# ========================
# Step 6: Compute Energies
# ========================

def bond_energy(bonds, coords):
    E = 0
    for i, j in bonds:
        r = np.linalg.norm(coords[i] - coords[j])
        E += 0.5 * FF["bond_k"] * (r - FF["r0"])**2
    return E

def angle_energy(angles, coords):
    E = 0
    for i, j, k in angles:
        v1 = coords[i] - coords[j]
        v2 = coords[k] - coords[j]
        theta = np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1)*np.linalg.norm(v2)))
        E += 0.5 * FF["angle_k"] * (theta - FF["theta0"])**2
    return E

def dihedral_energy(dihedrals, coords):
    E = 0
    for i, j, k, l in dihedrals:
        b1 = coords[j] - coords[i]
        b2 = coords[k] - coords[j]
        b3 = coords[l] - coords[k]

        n1 = np.cross(b1, b2)
        n2 = np.cross(b2, b3)

        n1 /= np.linalg.norm(n1)
        n2 /= np.linalg.norm(n2)

        m1 = np.cross(n1, b2/np.linalg.norm(b2))
        phi = np.arctan2(np.dot(m1, n2), np.dot(n1, n2))

        E += 0.5 * FF["dihedral_Vn"] * (1 + np.cos(FF["n"] * phi - FF["gamma"]))
    return E

def nonbonded_energy(atoms, coords, bonds):
    E = 0
    bonded_pairs = set(bonds)
    n = len(atoms)
    for i in range(n):
        for j in range(i+1, n):
            if (i, j) in bonded_pairs or (j, i) in bonded_pairs:
                continue  # skip bonded
            r = np.linalg.norm(coords[i] - coords[j])
            E_lj = 4 * FF["epsilon"] * ((FF["sigma"]/r)**12 - (FF["sigma"]/r)**6)
            E += E_lj
    return E

# ========================
# Step 7: Run Everything
# ========================
def compute_energy(filename):
    atoms, coords = read_xyz(filename)
    bonds = infer_bonds(atoms, coords)
    angles = generate_angles(bonds)
    dihedrals = generate_dihedrals(bonds)

    Eb = bond_energy(bonds, coords)
    Ea = angle_energy(angles, coords)
    Ed = dihedral_energy(dihedrals, coords)
    Enb = nonbonded_energy(atoms, coords, bonds)
    Etot = Eb + Ea + Ed + Enb

    print(f"--- Energy Breakdown for {filename} ---")
    print(f"Bonds:      {Eb:10.4f} kJ/mol")
    print(f"Angles:     {Ea:10.4f} kJ/mol")
    print(f"Dihedrals:  {Ed:10.4f} kJ/mol")
    print(f"Nonbonded:  {Enb:10.4f} kJ/mol")
    print(f"--------------------------------------")
    print(f"Total:      {Etot:10.4f} kJ/mol")
    print()
    return Etot

# Example: Run with ethanol.xyz

compute_energy("ethanol.xyz")

