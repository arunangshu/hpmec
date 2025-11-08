"""
High-Performance Molecular Energy Calculator
Module 4: Energy Calculation Engine

This module implements molecular mechanics energy calculations using:
- OPLS-AA force field format
- Harmonic bond stretching
- Harmonic angle bending  
- Fourier series dihedral torsions
- Lennard-Jones (12-6) van der Waals
- Coulombic electrostatics

UNITS CONVENTION (strictly enforced):
- Energy: kJ/mol
- Distance: nanometers (nm)
- Angles: radians
- Charge: elementary charge (e)
- LJ sigma: nanometers (nm)
- LJ epsilon: kJ/mol

NOTE: XYZ input files should have coordinates in Ångströms (standard XYZ format).
      The code automatically converts to nm internally.

Authors: Team Member 4 (Ankana) - Core calculations
         Integration with GUI by Team Member 5
"""

import numpy as np
import yaml
import itertools
import multiprocessing as mp
import time
from functools import partial
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from rdkit import Chem
from rdkit.Chem import AllChem

# Define global constants (UNITS: kJ/mol, nm, radians, elementary charge (e))
COULOMB_CONST = 1389.35458  # 1 / (4 * pi * epsilon0) * (e^2 / nm) to kJ/mol

# --- 1. CORE DATA STRUCTURES ---

class Atom:
    """Holds all data for a single atom."""
    def __init__(self, index, element, coords):
        self.index = index          # int
        self.element = element      # str (e.g., 'C', 'H')
        self.coords = coords        # tuple (x, y, z)
        self.atom_type = None       # str (e.g., 'opls_157')
        self.charge = 0.0           # float
        self.sigma = 0.0            # float (nm)
        self.epsilon = 0.0          # float (kJ/mol)

class Molecule:
    """Holds all data for the entire molecule."""
    def __init__(self):
        self.atoms = []             # List of Atom objects
        self.coordinates = None     # Nx3 NumPy array of all coords
        self.rdkit_mol = None       # RDKit Mol object

        # Topology and Exclusions
        self.bonds = []             # List of (i, j) tuples (original indexing)
        self.angles = []            # List of (i, j, k) tuples (j=center)
        self.dihedrals = []         # List of (i, j, k, l) tuples
        self.non_bonded_exclusions = set() # Set of (i, j) tuples for 1-2, 1-3, and 1-4 pairs
        self.orig_to_rdkit = {}     # mapping from original atom index -> rdkit atom index

# --- 2. INPUT/OUTPUT AND TOPOLOGY (Refactored) ---

def load_xyz(xyz_file_path):
    """
    Parses a .xyz file and returns a Molecule object.
    
    XYZ files are expected to have coordinates in Ångströms (standard XYZ format).
    Coordinates are automatically converted to nanometers for internal calculations.
    """
    with open(xyz_file_path, 'r') as f:
        lines = f.readlines()

    # Handle case where file might be empty or malformed
    if not lines or not lines[0].strip().isdigit():
        raise ValueError("Invalid or empty XYZ file.")

    num_atoms = int(lines[0].strip())
    coords_list = []
    mol = Molecule()

    atom_lines = lines[2:2 + num_atoms]
    for i, line in enumerate(atom_lines):
        parts = line.split()
        if len(parts) < 4: 
            continue  # Skip malformed lines

        element = parts[0]
        try:
            # Read coordinates in Ångströms (standard XYZ format)
            x_angstrom = float(parts[1])
            y_angstrom = float(parts[2])
            z_angstrom = float(parts[3])
            # Convert to nanometers for internal use (1 Å = 0.1 nm)
            coords = (x_angstrom * 0.1, y_angstrom * 0.1, z_angstrom * 0.1)
        except ValueError:
            raise ValueError(f"Invalid coordinate format in line {i+3}: {line.strip()}")

        mol.atoms.append(Atom(index=i, element=element, coords=coords))
        coords_list.append(coords)

    # Validate atom count
    if len(coords_list) != num_atoms:
        raise ValueError(f"XYZ file header says {num_atoms} atoms but found {len(coords_list)} coordinate lines")
    
    mol.coordinates = np.array(coords_list)
    return mol


def load_force_field(yaml_file_path):
    """Parses a.yaml force field file and returns a dict."""
    with open(yaml_file_path, 'r') as f:
        ff_parameters = yaml.safe_load(f)
    return ff_parameters


# --- Helper functions for robust RDKit mapping ---

def get_mol_with_mapping(atom_symbols, coords, add_hs=True):
    """Build an RDKit Mol from atom_symbols and coords, return (rdkit_mol, original_to_rdkit_idx).

    original_to_rdkit_idx maps indices in the input list to atom indices in rdkit_mol.
    If add_hs=True, hydrogens will be added but mapping will reflect the heavy-atom indices correctly
    by assuming original atoms occupy the first N indices.
    """
    rw = Chem.RWMol()
    for sym in atom_symbols:
        rw.AddAtom(Chem.Atom(sym))
    # Create conformer and set positions
    conf = Chem.Conformer(rw.GetNumAtoms())
    for i, p in enumerate(coords):
        conf.SetAtomPosition(i, p)
    rw.AddConformer(conf, assignId=True)

    mol = rw.GetMol()

    if add_hs:
        # Chem.AddHs appends hydrogens AFTER existing atoms, so original indices remain 0..N-1
        mol = Chem.AddHs(mol, addCoords=True)

    # Try to embed and optimize (best-effort, failures are caught)
    try:
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        AllChem.UFFOptimizeMolecule(mol)
    except Exception:
        # embedding/optimization failed - continue with the molecule as-is
        pass

    # Build mapping: we assume original atoms are first N atoms in rdkit mol (true if AddHs appended H)
    original_to_rdkit = {i: i for i in range(len(atom_symbols))}
    return mol, original_to_rdkit


def xyz_to_mol(xyz_block_lines, charge=0):
    '''
    Convert an xyz block (as a list of lines) to an RDKit molecule and mapping.
    Returns: (rdkit_mol, orig_to_rdkit_map)
    '''
    num_atoms = int(xyz_block_lines[0].strip())
    atom_symbols = []
    coords = []
    for line in xyz_block_lines[2:2 + num_atoms]:
        parts = line.split()
        atom_symbols.append(parts[0])
        coords.append((float(parts[1]), float(parts[2]), float(parts[3])))

    rdkit_mol, mapping = get_mol_with_mapping(atom_symbols, coords, add_hs=True)
    return rdkit_mol, mapping


# ---------- Replacement: bond inference via covalent radii and RDKit construction ----------

# Simple covalent radii (Å). We'll convert to nm inside code if your coords are nm.
COVALENT_RADII_ANGSTROM = {
    'H': 0.31, 'He': 0.28,
    'Li': 1.28, 'Be': 0.96, 'B': 0.84, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57, 'Ne': 0.58,
    'Na': 1.66, 'Mg': 1.41, 'Al': 1.21, 'Si': 1.11, 'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Ar': 1.06,
}

def infer_bonds_by_distance(molecule, factor=1.2, coords_in_nm=True):
    """
    Infer bonds from coordinates using covalent radii.
    factor: multiplier for sum of radii to use as bond cutoff (1.1-1.25 typical).
    coords_in_nm: if True assumes molecule.coordinates are in nm and converts radii to nm.
    """
    coords = np.array([a.coords for a in molecule.atoms])
    n = len(molecule.atoms)
    mol_coords = coords
    # Convert radii to nm if coords are in nm (radii are given in angstrom above)
    radius_factor = 0.1 if coords_in_nm else 1.0
    radii = []
    for a in molecule.atoms:
        rA = COVALENT_RADII_ANGSTROM.get(a.element, 0.77)  # fallback approx (C-like)
        radii.append(rA * radius_factor)
    radii = np.array(radii)

    for i in range(n):
        for j in range(i+1, n):
            dij = np.linalg.norm(mol_coords[i] - mol_coords[j])
            cutoff = factor * (radii[i] + radii[j])
            if dij > 1e-6 and dij <= cutoff:
                molecule.bonds.append((i, j))
                # FIX: Always store exclusions as sorted tuples for consistent lookup
                molecule.non_bonded_exclusions.add(tuple(sorted((i, j))))

def build_rdkit_from_bonds(molecule):
    """
    Construct an RDKit Mol from the molecule atoms and inferred bonds (no AddHs).
    This enables SMARTS matching without running AddHs on an atom-only mol.
    """
    rw = Chem.RWMol()
    for atom in molecule.atoms:
        a = Chem.Atom(atom.element)
        rw.AddAtom(a)

    # Add single bonds for inferred connectivity (this is a heuristic)
    for i, j in molecule.bonds:
        try:
            rw.AddBond(int(i), int(j), Chem.BondType.SINGLE)
        except Exception:
            pass

    mol = rw.GetMol()
    # Attach coordinates as a conformer
    conf = Chem.Conformer(mol.GetNumAtoms())
    for idx, atom in enumerate(molecule.atoms):
        conf.SetAtomPosition(idx, atom.coords)
    mol.AddConformer(conf, assignId=True)

    # Sanitize to calculate implicit valences, aromaticity, etc.
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        # Sanitization can fail for odd valences; still return the mol to allow substructure matching
        pass

    return mol


def infer_topology(molecule):
    """
    New infer_topology: infer bonds with distance heuristics, then angles/dihedrals from connectivity.
    Build an RDKit Mol from those bonds for SMARTS-based atom typing (no AddHs on incomplete mol).
    """
    # ensure coordinates array exists
    if molecule.coordinates is None:
        molecule.coordinates = np.array([a.coords for a in molecule.atoms])

    # Reset topology
    molecule.bonds = []
    molecule.angles = []
    molecule.dihedrals = []
    molecule.non_bonded_exclusions = set()
    molecule.rdkit_mol = None
    molecule.orig_to_rdkit = {}

    # 1) Infer bonds
    infer_bonds_by_distance(molecule, factor=1.2, coords_in_nm=True)

    # 2) Build adjacency (neighbors) list
    neighbors = {i: set() for i in range(len(molecule.atoms))}
    for i, j in molecule.bonds:
        neighbors[i].add(j)
        neighbors[j].add(i)

    # 3) Find angles (i - j - k) where j is center
    for j in range(len(molecule.atoms)):
        nbrs = sorted(neighbors[j])
        for a_idx, b_idx in itertools.combinations(nbrs, 2):
            molecule.angles.append((a_idx, j, b_idx))
            molecule.non_bonded_exclusions.add(tuple(sorted((a_idx, b_idx))))

    # 4) Find dihedrals i-j-k-l by traversing bonds (j-k bond defines central)
    for j, k in molecule.bonds:
        j_nbrs = [n for n in neighbors[j] if n != k]
        k_nbrs = [n for n in neighbors[k] if n != j]
        for i in j_nbrs:
            for l in k_nbrs:
                molecule.dihedrals.append((i, j, k, l))
                molecule.non_bonded_exclusions.add(tuple(sorted((i, l))))

    # 5) Build a sanitized RDKit Mol from inferred bonds to enable SMARTS matching
    try:
        rdkit_mol = build_rdkit_from_bonds(molecule)
        molecule.rdkit_mol = rdkit_mol
        # Map original indices trivially (we didn't add atoms)
        molecule.orig_to_rdkit = {i: i for i in range(len(molecule.atoms))}
    except Exception:
        molecule.rdkit_mol = None
        molecule.orig_to_rdkit = {}

# --- 3. PARAMETER ASSIGNMENT (IMPROVED) (IMPROVED) ---

def find_bond_param(ff_bonds, type_i, type_j):
    key1 = f"{type_i}-{type_j}"
    key2 = f"{type_j}-{type_i}"
    if key1 in ff_bonds:
        return ff_bonds[key1]
    if key2 in ff_bonds:
        return ff_bonds[key2]
    return None


def find_angle_param(ff_angles, type_i, type_j_center, type_k):
    key1 = f"{type_i}-{type_j_center}-{type_k}"
    key2 = f"{type_k}-{type_j_center}-{type_i}"
    if key1 in ff_angles:
        return ff_angles[key1]
    if key2 in ff_angles:
        return ff_angles[key2]
    return None


def assign_parameters(molecule, ff_parameters):
    """Assigns atom types and parameters with robust SMARTS matching and index mapping."""
    rdkit_mol = molecule.rdkit_mol
    orig_to_rdkit = getattr(molecule, 'orig_to_rdkit', {a.index: a.index for a in molecule.atoms})
    atom_type_rules = ff_parameters.get('atom_types', [])
    param_maps = {'bonds': {}, 'angles': {}, 'dihedrals': {}}

    # Reset atom fields
    for atom in molecule.atoms:
        atom.atom_type = None
        atom.charge = 0.0
        atom.sigma = 0.0
        atom.epsilon = 0.0

    # Build reverse mapping: rdkit_idx -> original_idx
    rd_to_orig = {rd: orig for orig, rd in orig_to_rdkit.items()}

    # Atom Typing: iterate rules and apply to all matched RDKit atoms
    for rule in atom_type_rules:
        smarts = rule.get('smarts')
        try:
            patt = Chem.MolFromSmarts(smarts)
            if patt is None:
                print(f"Warning: invalid SMARTS pattern: {smarts}")
                continue
        except Exception as e:
            print(f"Warning: smarts compile error for {smarts}: {e}")
            continue

        matches = rdkit_mol.GetSubstructMatches(patt, useChirality=False)
        for match in matches:
            for rd_idx in match:
                orig_idx = rd_to_orig.get(rd_idx, None)
                if orig_idx is None:
                    continue
                atom = molecule.atoms[orig_idx]
                if atom.atom_type is None:
                    atom.atom_type = rule.get('type_name')
                    atom.charge = rule.get('charge', atom.charge)
                    atom.sigma = rule.get('sigma', atom.sigma)
                    atom.epsilon = rule.get('epsilon', atom.epsilon)

    # FIX: Add fallback typing for unassigned atoms to prevent None-type failures
    # Default LJ parameters are carbon-like values (commonly used as fallback in MD)
    FALLBACK_PARAMS = {
        'H': {'sigma': 0.25, 'epsilon': 0.126},  # Hydrogen-like
        'C': {'sigma': 0.35, 'epsilon': 0.276},  # Carbon-like
        'N': {'sigma': 0.325, 'epsilon': 0.711}, # Nitrogen-like
        'O': {'sigma': 0.296, 'epsilon': 0.879}, # Oxygen-like
        'DEFAULT': {'sigma': 0.35, 'epsilon': 0.3}  # Generic fallback
    }
    
    for atom in molecule.atoms:
        if atom.atom_type is None:
            print(f"Warning: Atom {atom.index} {atom.element} has no assigned atom_type. Using fallback generic_{atom.element}")
            atom.atom_type = f"generic_{atom.element}"
            atom.charge = 0.0  # Neutral fallback
            fallback = FALLBACK_PARAMS.get(atom.element, FALLBACK_PARAMS['DEFAULT'])
            atom.sigma = fallback['sigma']
            atom.epsilon = fallback['epsilon']

    # Map bonds (all atoms now guaranteed to have types via fallback)
    for i, j in molecule.bonds:
        type_i = molecule.atoms[i].atom_type
        type_j = molecule.atoms[j].atom_type
        params = find_bond_param(ff_parameters.get('bond_types', {}), type_i, type_j)
        if params is not None:
            param_maps['bonds'][(i, j)] = params
        else:
            print(f"Warning: Missing bond params for {type_i}-{type_j}")

    # Map angles (all atoms now guaranteed to have types via fallback)
    for i, j_center, k in molecule.angles:
        t_i = molecule.atoms[i].atom_type
        t_j = molecule.atoms[j_center].atom_type
        t_k = molecule.atoms[k].atom_type
        params = find_angle_param(ff_parameters.get('angle_types', {}), t_i, t_j, t_k)
        if params is not None:
            param_maps['angles'][(i, j_center, k)] = params
        else:
            print(f"Warning: Missing angle params for {t_i}-{t_j}-{t_k}")

    # Map dihedrals (checking forward and reverse keys, all atoms guaranteed to have types)
    ff_dihedrals = ff_parameters.get('dihedral_types', {})
    for i, j, k, l in molecule.dihedrals:
        type_i = molecule.atoms[i].atom_type
        type_j = molecule.atoms[j].atom_type
        type_k = molecule.atoms[k].atom_type
        type_l = molecule.atoms[l].atom_type
        key = f"{type_i}-{type_j}-{type_k}-{type_l}"
        key_rev = f"{type_l}-{type_k}-{type_j}-{type_i}"
        if key in ff_dihedrals:
            param_maps['dihedrals'][(i, j, k, l)] = ff_dihedrals[key]
        elif key_rev in ff_dihedrals:
            param_maps['dihedrals'][(i, j, k, l)] = ff_dihedrals[key_rev]
        else:
            # it's common that FF dihedrals are specified in some canonical form; warn if missing
            print(f"Warning: Missing dihedral params for {key} (or reverse)")

    return param_maps

# --- 4. GEOMETRY UTILITIES (IMPROVED) ---

def get_distance(coords, i, j):
    """Calculate the Euclidean distance between two atoms."""
    return np.linalg.norm(coords[i] - coords[j])


def get_angle(coords, i, j, k):
    """Calculate the bond angle (i-j-k) in radians."""
    r_ji = coords[i] - coords[j]
    r_jk = coords[k] - coords[j]
    dot_prod = np.dot(r_ji, r_jk)
    norm_ji = np.linalg.norm(r_ji)
    norm_jk = np.linalg.norm(r_jk)
    cos_theta = dot_prod / (norm_ji * norm_jk + 1e-12)
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    return np.arccos(cos_theta)


def get_dihedral(coords, i, j, k, l):
    """Robust dihedral calculation returning signed angle in radians."""
    p0 = coords[i]
    p1 = coords[j]
    p2 = coords[k]
    p3 = coords[l]
    b0 = -1.0 * (p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    b1_norm = np.linalg.norm(b1)
    if b1_norm < 1e-12:
        return 0.0
    b1 = b1 / (b1_norm)

    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1

    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.arctan2(y, x)

# --- 5. CORE ENERGY CALCULATION (NEW BONDED) ---

def calculate_bond_energy(molecule, param_maps):
    """
    Calculates harmonic bond stretching energy.
    E_bond = sum(0.5 * k_b * (r - b0)^2)  OR as user had k*(r-b0)^2 depending on FF convention.
    Here we will assume parameters are given such that energy = k * (r - b0)^2 as in original demo.
    """
    energy = 0.0
    coords = molecule.coordinates

    for (i, j), params in param_maps['bonds'].items():
        k_b, b0 = params
        r_ij = get_distance(coords, i, j)
        energy += k_b * (r_ij - b0)**2

    return energy


def calculate_angle_energy(molecule, param_maps):
    """
    Calculates harmonic angle bending energy.
    E_angle = sum(k_theta * (theta - theta0)^2)
    """
    energy = 0.0
    coords = molecule.coordinates

    for (i, j, k), params in param_maps['angles'].items():
        k_theta, theta0 = params
        theta_ijk = get_angle(coords, i, j, k)
        energy += k_theta * (theta_ijk - theta0)**2

    return energy


def calculate_dihedral_energy(molecule, param_maps):
    """
    Calculates dihedral (torsion) energy using a multi-term Fourier series (OPLS-like).
    """
    energy = 0.0
    coords = molecule.coordinates

    for (i, j, k, l), params in param_maps['dihedrals'].items():
        # Params are V1, V2, V3, V4 (as defined in the dummy YAML)
        # Ensure params is list-like before concatenation (ff may provide tuples)
        p_list = list(params) if not isinstance(params, list) else params
        V1, V2, V3, V4 = (p_list + [0.0, 0.0, 0.0, 0.0])[:4]
        phi = get_dihedral(coords, i, j, k, l)

        e_torsion = (
            0.5 * V1 * (1 + np.cos(1 * phi)) +
            0.5 * V2 * (1 - np.cos(2 * phi)) +
            0.5 * V3 * (1 + np.cos(3 * phi)) +
            0.5 * V4 * (1 - np.cos(4 * phi))
        )
        energy += e_torsion

    return energy


def calculate_nonbonded_optimized(molecule, cutoff=1.0):
    """
    Calculates non-bonded energy using k-d tree neighbor search (O(N log N)).
    """
    energy = 0.0
    coords = molecule.coordinates
    atoms = molecule.atoms

    # 1. Build the k-d tree from coordinates
    tree = cKDTree(coords)

    # 2. Find all unique pairs (i, j) within the cutoff radius.
    # Use output_type='set' when available (newer SciPy); otherwise fall back and convert.
    try:
        pairs = tree.query_pairs(r=cutoff, output_type='set')
    except TypeError:
        raw_pairs = tree.query_pairs(r=cutoff)
        # query_pairs may return a set of frozensets or a set of tuples; normalize to set of (i,j)
        pairs = set()
        for pair in raw_pairs:
            try:
                i, j = pair
            except Exception:
                # If pair is frozenset-like
                lst = list(pair)
                if len(lst) != 2:
                    continue
                i, j = lst
            pairs.add((min(i, j), max(i, j)))

    for i, j in pairs:
        i, j = min(i, j), max(i, j)
        if (i, j) in molecule.non_bonded_exclusions:
            continue

        atom_i = atoms[i]
        atom_j = atoms[j]

        r_ij = get_distance(coords, i, j)
        if r_ij < 1e-12: continue

        # VDW (Lennard-Jones) — combining rules (Lorentz-Berthelot)
        sigma_ij = (atom_i.sigma + atom_j.sigma) / 2.0
        epsilon_ij = np.sqrt(max(0.0, atom_i.epsilon * atom_j.epsilon))
        if sigma_ij > 0 and epsilon_ij > 0:
            r_ratio = sigma_ij / r_ij
            r6 = r_ratio**6
            r12 = r6**2
            energy += 4.0 * epsilon_ij * (r12 - r6)

        # Coulomb (Electrostatic)
        energy += COULOMB_CONST * (atom_i.charge * atom_j.charge) / r_ij

    return energy


# --- 6. GUI INTEGRATION WRAPPER FUNCTIONS ---

def calculate_single_molecule_energy(xyz_file_path, ff_yaml_path):
    """
    Main entry point for GUI integration.
    
    Orchestrates the complete energy calculation pipeline:
    1. Load XYZ file
    2. Infer molecular topology (bonds, angles, dihedrals)
    3. Assign force field parameters via SMARTS matching
    4. Calculate all energy components
    5. Return total energy
    
    Args:
        xyz_file_path: Path to .xyz molecular geometry file (Ångström coordinates)
        ff_yaml_path: Path to .yaml force field parameters file (OPLS-AA format)
    
    Returns:
        Tuple of (filename, total_energy_in_kJ_per_mol)
    
    Raises:
        ValueError: If XYZ file invalid or force field missing
        Exception: For other calculation errors
    """
    try:
        # Load molecule and force field
        mol = load_xyz(xyz_file_path)
        ff = load_force_field(ff_yaml_path)
        
        # Build topology
        infer_topology(mol)
        
        # Assign parameters
        param_maps = assign_parameters(mol, ff)
        
        # Calculate energy components
        E_bond = calculate_bond_energy(mol, param_maps)
        E_angle = calculate_angle_energy(mol, param_maps)
        E_dihedral = calculate_dihedral_energy(mol, param_maps)
        E_nonbonded = calculate_nonbonded_optimized(mol, cutoff=1.0)
        
        # Total energy
        total_energy = E_bond + E_angle + E_dihedral + E_nonbonded
        
        return (xyz_file_path, total_energy)
    
    except Exception as e:
        print(f"Error calculating energy for {xyz_file_path}: {e}")
        raise


def run_parallel_calculations(list_of_xyz_files, ff_yaml_path):
    """
    Parallel calculation wrapper for GUI benchmarking.
    
    Uses multiprocessing.Pool to distribute energy calculations
    across multiple CPU cores for performance testing.
    
    Args:
        list_of_xyz_files: List of paths to .xyz files
        ff_yaml_path: Path to force field .yaml file
    
    Returns:
        List of (filename, energy) tuples
    """
    # Create partial function with fixed ff_yaml_path
    worker = partial(calculate_single_molecule_energy, ff_yaml_path=ff_yaml_path)
    
    # Use multiprocessing Pool for parallel execution
    with mp.Pool() as pool:
        results = pool.map(worker, list_of_xyz_files)
    
    return results


if __name__ == '__main__':
    # Minimal module sanity check when run directly.
    print('calculator.py loaded. Functions available: load_xyz, infer_topology, assign_parameters, calculate_bond_energy, calculate_angle_energy, calculate_dihedral_energy, calculate_nonbonded_optimized')
    print('GUI integration functions: calculate_single_molecule_energy, run_parallel_calculations')
