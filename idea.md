

# **Project Blueprint: A High-Performance Molecular Energy Calculator**

## **1\. Foundational Theory and System Architecture**

This document provides a comprehensive blueprint for developing a high-performance molecular energy calculator in Python. It begins by establishing the necessary scientific theory, then presents a robust 5-module software architecture designed to separate complex tasks. This structure improves upon the proposed project plan by explicitly isolating the critical (and often overlooked) challenges of topology inference and parameter assignment from the core energy calculation.

### **1.1 The Classical Potential Energy Function (The "Force Field")**

In molecular mechanics (MM), the Born-Oppenheimer approximation is assumed, allowing the total potential energy ($E\_{\\text{Total}}$) of a molecular system to be calculated as a function of the nuclear coordinates using a classical potential, or "force field". The total energy is the sum of two major components: bonded (internal) interactions and non-bonded (external) interactions.2

$$E\_{\\text{Total}} \= E\_{\\text{Bonded}} \+ E\_{\\text{Non-Bonded}}$$

#### **1.1.1 Bonded Interactions ($E\_{\\text{Bonded}}$)**

These terms represent the energy stored in the covalent structure of the molecule.2

* **Bond Stretching:** This term models the energy required to deform a covalent bond (between atoms $i, j$) from its equilibrium length, $b\_0$. It is most commonly modeled as a harmonic oscillator, similar to a simple spring.4  
  $$V\_{\\text{bond}}(b) \= \\sum\_{\\text{bonds}} \\frac{1}{2} k\_b (b \- b\_0)^2$$

  The required parameters are the force constant $k\_b$ and the equilibrium bond length $b\_0$ for each unique bond type (e.g., C-C, C-O).  
* **Angle Bending:** This term models the energy required to deform the angle ($\\theta$) formed by three consecutively bonded atoms ($i-j-k$) from its equilibrium angle, $\\theta\_0$. This is also typically modeled with a harmonic potential.2  
  $$V\_{\\text{angle}}(\\theta) \= \\sum\_{\\text{angles}} \\frac{1}{2} k\_\\theta (\\theta \- \\theta\_0)^2$$

  The required parameters are the angle force constant $k\_\\theta$ and the equilibrium angle $\\theta\_0$ for each unique angle type (e.g., C-C-H, C-O-H).  
* **Dihedral (Torsional) Twisting:** This term models the energy barrier associated with rotation around a central bond (connecting atoms $j-k$ in a $i-j-k-l$ system). This is a periodic function, typically represented by a Fourier series.6 A single dihedral interaction may be the sum of multiple cosine terms with different periodicities ($n$).6  
  $$V\_{\\text{dihedral}}(\\phi) \= \\sum\_{\\text{dihedrals}} \\sum\_{n} K\_{\\phi,n} (1 \+ \\cos(n\\phi \- \\delta\_n))$$

  The required parameters for each dihedral type (e.g., H-C-C-O) are a set of amplitudes $K\_{\\phi,n}$, periodicities $n$, and phase offsets $\\delta\_n$.3

#### **1.1.2 Non-Bonded Interactions ($E\_{\\text{Non-Bonded}}$)**

These terms describe the "through-space" interactions between all pairs of atoms ($i, j$) that are not already part of the bonded structure (typically, 1-2, 1-3, and sometimes 1-4 interactions are excluded or scaled).2 This calculation is the primary computational bottleneck in molecular simulations, as a naïve implementation scales quadratically, $O(N^2)$, with the number of atoms $N$.9

* **Van der Waals (VDW):** This term accounts for two distinct forces: strong, short-range repulsion (Pauli exclusion principle) and weaker, long-range attraction (London dispersion forces).11 The most common model is the Lennard-Jones (LJ) 12-6 potential.12  
  $$V\_{\\text{LJ}}(r\_{ij}) \= \\sum\_{i\<j} 4\\epsilon\_{ij} \\left\[ \\left(\\frac{\\sigma\_{ij}}{r\_{ij}}\\right)^{12} \- \\left(\\frac{\\sigma\_{ij}}{r\_{ij}}\\right)^6 \\right\]$$

  The required parameters are the potential well depth $\\epsilon$ and the collision diameter $\\sigma$ for each atom type. Parameters for mixed pairs ($i, j$) are typically calculated from individual atom parameters using "combining rules," such as the Lorentz-Berthelot rules.14  
* **Electrostatic (Coulomb):** This term models the electrostatic interaction between atoms based on a distribution of fixed, partial atomic charges ($q$).3 It is a direct application of Coulomb's Law.17  
  $$V\_{\\text{Coulomb}}(r\_{ij}) \= \\sum\_{i\<j} \\frac{1}{4\\pi\\epsilon\_0 D} \\frac{q\_i q\_j}{r\_{ij}}$$

  The required parameters are the partial charge $q\_i$ for each atom type.3 The dielectric constant $D$ is typically set to 1.0 for calculations in a vacuum.

### **1.2 A Revised Project Architecture: The 5-Module Pipeline**

The initial project plan correctly identifies the key domains: I/O, core logic, and HPC optimization. However, a robust implementation must first solve a critical *cheminformatics* problem that precedes the "core logic": an input .xyz file contains only atomic elements and coordinates, providing no information on which atoms are bonded, what the angles are, or what the dihedrals are.20

Therefore, the "core logic" must be expanded into a 5-module pipeline that explicitly manages this data-wrangling and inference challenge. This modular design is more scalable and provides a clearer division of labor for the team.

* **Module 1: Input / Output (I/O) & Data Structures:** This module defines the central Python classes (e.g., Molecule, Atom) and contains parsers for the .xyz molecular geometry 21 and .yaml force field 24 input files.  
* **Module 2: Topology Inference Engine:** This is the *first* part of the core logic. It consumes the raw coordinates and element list from Module 1 and uses cheminformatics rules to *infer* the molecular graph, populating the Molecule object with lists of bonds, angles, and dihedrals.25  
* **Module 3: Parameter Assignment Engine:** This is the *second* part of the core logic. It acts as an "expert system" that inspects the inferred topology from Module 2, determines the correct "atom type" for each atom (e.g., distinguishing a methyl carbon CT from an aromatic carbon CA) 28, and then queries the force field data (from Module 1\) to assign the specific physical parameters ($k\_b$, $q\_i$, $\\epsilon\_i$, etc.) to every atom, bond, angle, and dihedral.30  
* **Module 4: Core Energy Calculator (HPC-1):** This module contains the *mathematical* implementation of the energy functions from Section 1.1. It takes the fully parameterized Molecule object and calculates the total energy. This module will include both a naïve $O(N^2)$ non-bonded function and the optimized $O(N)$ single-core solution.32  
* **Module 5: Parallelization & Interface (HPC-2 & GUI):** This module addresses the two scaling-related tasks. It implements the multi-core parallelization (HPC-2) for calculating the energies of *many* molecules simultaneously.34 It also contains the simple graphical user interface (GUI) for system operation.36

### **1.3 Recommended 4-Member Project Division**

This 5-module architecture maps cleanly onto a 4-person team, minimizing dependencies and creating distinct, reportable domains for each member.37

**Table 1: Recommended 4-Member Project Role and Task Division**

| Member | Role | Key Modules | Core Responsibilities & Deliverables |
| :---- | :---- | :---- | :---- |
| **Member 1** | **System Architect & Topology Lead** | Module 1 & 2 | \- Define and code the central Molecule and Atom Python classes. \- Implement the .xyz file parser.21 \- Implement the **Topology Inference Engine:** Use xyz2mol and RDKit libraries 27 to construct the molecular graph from coordinates.39 \- **Deliverable:** A Molecule object populated with atomic coordinates and inferred lists of bonds, angles, and dihedrals. |
| **Member 2** | **Force Field & Parameterization Lead** | Module 1 & 3 | \- Define the full .yaml force field file structure.41 \- Implement the .yaml file parser.24 \- Implement the **Parameter Assignment Engine:** Use SMARTS rules 43 to perform "atom typing".28 \- Implement the parameter lookup logic.31 \- **Deliverable:** The final ethanol.yaml test file and a Molecule object populated with all parameters (charges, VDW, bonded constants). |
| **Member 3** | **Core Logic & HPC-1 Lead** | Module 4 | \- Implement all bonded energy functions (calc\_bonds, calc\_angles, calc\_dihedrals) using NumPy.3 \- Implement the **naïve** $O(N^2)$ non-bonded energy function.10 \- Implement the **optimized** (HPC-1) non-bonded function using scipy.spatial.cKDTree as a highly efficient cell/neighbor list.32 \- **Deliverable:** The core calculate\_energy(molecule) function and benchmark data (naïve vs. optimized). |
| **Member 4** | **Parallelization & GUI Lead** | Module 5 | \- Implement the **HPC-2** scaling solution: Use multiprocessing.Pool.map 35 to distribute a list of molecule calculations across all available CPU cores.34 \- Implement the **GUI:** Use **Streamlit** 36 (recommended for its simplicity over PyQt 50) for file upload and results display. \- **Role:** Act as the primary **Integrator**, managing the main Git repository and combining the modules from all members. |

---

## **2\. Modules 1 & 2: Input, Topology, and Force Field Design**

This section details the specific implementation plan for **Member 1 (Topology)** and **Member 2 (Force Field)**.

### **2.1 Module 1: Input Parsers and Data Structures**

The foundation of the project is a set of robust data classes.

* **Central Data Classes:** The Molecule class will act as the central data container that is passed between all modules.  
  Python  
  import numpy as np

  class Atom:  
      """Holds all data for a single atom."""  
      def \_\_init\_\_(self, index, element, coords):  
          self.index \= index        \# int  
          self.element \= element    \# str (e.g., 'C', 'H')  
          self.coords \= coords      \# tuple (x, y, z)

          \# \--- To be populated by Module 3 (Parameterization) \---  
          self.atom\_type \= None     \# str (e.g., 'opls\_157')  
          self.charge \= 0.0         \# float  
          self.sigma \= 0.0          \# float  
          self.epsilon \= 0.0        \# float

  class Molecule:  
      """Holds all data for the entire molecule."""  
      def \_\_init\_\_(self):  
          self.atoms \=           \# List of Atom objects  
          self.coordinates \= None   \# Nx3 NumPy array of all coords  
          self.rdkit\_mol \= None     \# RDKit Mol object (from Module 2\)

          \# \--- To be populated by Module 2 (Topology) \---  
          self.bonds \=           \# List of (i, j) tuples  
          self.angles \=          \# List of (i, j, k) tuples (j=center)  
          self.dihedrals \=       \# List of (i, j, k, l) tuples

          \# Set of (i, j) tuples for 1-2, 1-3, and 1-4 pairs  
          self.non\_bonded\_exclusions \= set()

* **.xyz Parser (Member 1):** This function reads the standard .xyz file format 21 and populates the Molecule object.  
  Python  
  def load\_xyz(xyz\_file\_path):  
      """Parses a.xyz file and returns a Molecule object."""  
      with open(xyz\_file\_path, 'r') as f:  
          lines \= f.readlines()

      num\_atoms \= int(lines.strip())  
      coords\_list \=  
      mol \= Molecule()

      atom\_lines \= lines\[2:2 \+ num\_atoms\]  
      for i, line in enumerate(atom\_lines):  
          parts \= line.split()  
          element \= parts  
          coords \= (float(parts), float(parts), float(parts))  
          mol.atoms.append(Atom(index=i, element=element, coords=coords))  
          coords\_list.append(coords)

      mol.coordinates \= np.array(coords\_list)  
      return mol

* **.yaml Parser (Member 2):** This function reads the force field file using the PyYAML library.24 The safe\_load function is used to prevent security vulnerabilities.53  
  Python  
  import yaml

  def load\_force\_field(yaml\_file\_path):  
      """Parses a.yaml force field file and returns a dict."""  
      with open(yaml\_file\_path, 'r') as f:  
          ff\_parameters \= yaml.safe\_load(f)  
      return ff\_parameters

### **2.2 Module 2: The Topology Inference Engine (Member 1\)**

This module performs the most complex cheminformatics task: determining the molecular graph (bonds, angles, dihedrals) from coordinates alone. Re-inventing this logic is highly error-prone.55 The recommended solution is to leverage the xyz2mol and RDKit libraries.27 xyz2mol is specifically designed to infer bond orders from .xyz files and produce a high-quality RDKit.Mol object.57

* **Installation:** pip install rdkit-pypi xyz2mol  
* **Core Logic:** The following function takes the partially-filled Molecule object, uses xyz2mol to create an RDKit.Mol object, and then traverses this object to find all bonds, angles, and dihedrals.  
  Python  
  from rdkit import Chem  
  from rdkit.Chem import AllChem  
  from xyz2mol import xyz\_to\_mol  
  import itertools

  def infer\_topology(molecule):  
      """  
      Infers bonds, angles, and dihedrals from atom coordinates.  
      Populates the.bonds,.angles, and.dihedrals lists in the   
      Molecule object.  
      """  
      \# 1\. Create an xyz-formatted string block from our Molecule object  
      xyz\_block \= f"{len(molecule.atoms)}\\n\\n"  
      for atom in molecule.atoms:  
          xyz\_block \+= f"{atom.element} {atom.coords} {atom.coords} {atom.coords}\\n"

      \# 2\. Use xyz2mol to get an RDKit Mol object \[27, 57\]  
      \# We assume a neutral charge (0) for this project.  
      try:  
          \# xyz\_to\_mol returns a list of Mol objects; we take the first  
          rdkit\_mol \= xyz\_to\_mol(xyz\_block.split('\\n'), charge=0)  
      except Exception as e:  
          print(f"Error in xyz2mol topology inference: {e}")  
          return

      \# 3\. Add explicit hydrogens. This is a critical step, as   
      \#    force fields require them for correct energy.\[44\]  
      rdkit\_mol \= Chem.AddHs(rdkit\_mol, addCoords=True)

      \# NOTE: At this point, the atom list in \`molecule\` and \`rdkit\_mol\`  
      \# may differ. The system must be re-synchronized. For this project,  
      \# we will assume the initial.xyz \*contains\* all atoms,   
      \# including hydrogens, and \`xyz2mol\` just finds the bonds.  
      \# A more robust solution would rebuild the \`molecule.atoms\` list  
      \# from the \`rdkit\_mol\` object.  
      molecule.rdkit\_mol \= rdkit\_mol

      \# 4\. Find and store all bonds   
      for bond in rdkit\_mol.GetBonds():  
          i \= bond.GetBeginAtomIdx()  
          j \= bond.GetEndAtomIdx()  
          molecule.bonds.append(tuple(sorted((i, j))))  
          \# Add 1-2 pairs to exclusion list  
          molecule.non\_bonded\_exclusions.add(tuple(sorted((i, j))))

      \# 5\. Find and store all angles  
      \# Angles are found by finding all atoms with \>= 2 neighbors \[58\]  
      for atom in rdkit\_mol.GetAtoms():  
          j \= atom.GetIdx() \# j is the central atom  
          neighbors \= \[n.GetIdx() for n in atom.GetNeighbors()\]  
          if len(neighbors) \< 2:  
              continue

          \# Find all unique combinations of 2 neighbors (i, k)  
          for i, k in itertools.combinations(neighbors, 2):  
              molecule.angles.append((i, j, k)) \# Store as (a1, center, a3)  
              \# Add 1-3 pairs to exclusion list  
              molecule.non\_bonded\_exclusions.add(tuple(sorted((i, k))))

      \# 6\. Find and store all dihedrals (i-j-k-l)  
      \# Iterate over all central bonds (j-k) \[40, 59\]  
      for bond in rdkit\_mol.GetBonds():  
          j \= bond.GetBeginAtomIdx()  
          k \= bond.GetEndAtomIdx()

          j\_neighbors \=  
          k\_neighbors \=

          if not j\_neighbors or not k\_neighbors:  
              continue

          \# Find all i-j-k-l combinations  
          for i in j\_neighbors:  
              for l in k\_neighbors:  
                  molecule.dihedrals.append((i, j, k, l))  
                  \# Add 1-4 pairs to exclusion list  
                  molecule.non\_bonded\_exclusions.add(tuple(sorted((i, l))))

---

## **3\. Module 3: The Parameter Assignment Engine (Member 2\)**

This module *connects* the inferred topology (Module 2\) to the force field parameters (Module 1). This involves two steps: (1) defining the force field data structure and (2) implementing the logic to assign parameters based on that structure.

### **3.1 The Force Field .yaml Structure**

The .yaml file must be structured as a queryable database.41 The structure below is based on the OPLS-AA (Optimized Potentials for Liquid Simulations, All-Atom) force field, which is well-suited for organic molecules like ethanol.61

Proposed ethanol.yaml (Dummy Force Field):  
This file, synthesized from parameters found in OPLS-AA files and documentation 62, serves as the required test file.

YAML

\# \--- ethanol.yaml \---  
\# A complete OPLS-AA-like dummy force field for C2H5OH (Ethanol)  
\# for testing the HPC project.  
\# UNITS: kJ/mol, nm, radians, elementary charge (e)  
\# \-----------------------------------------------------------------

\# Atom typing rules. Processed in order. Most specific first.  
\# SMARTS  define the chemical environment.  
\# Parameters from \[67\] (charges) and representative  
\# VDW params.\[65, 66, 68, 69\]  
atom\_types:  
  \- smarts: '\[H\]\[OX2H1\](\[CX4H2\])'  \# H on an alcohol: HO (opls\_155)  
    type\_name: 'opls\_155'  
    charge: 0.418  
    sigma: 0.0       
    epsilon: 0.0   

  \- smarts: '\[OX2H1\](\[H\])\[CX4H2\]'   \# O in an alcohol: OH (opls\_154)  
    type\_name: 'opls\_154'  
    charge: \-0.683  
    sigma: 0.312  
    epsilon: 0.71132

  \- smarts: '\[CX4H2\](\[OX2H1\])\[CH3X4\]' \# C in CH2-OH: CT (opls\_157)  
    type\_name: 'opls\_157'  
    charge: 0.145  
    sigma: 0.350  
    epsilon: 0.276144

  \- smarts: '\[CH3X4\](\[CX4H2\])'       \# C in CH3: CT (opls\_157)  
    type\_name: 'opls\_157'  
    charge: \-0.180  
    sigma: 0.350  
    epsilon: 0.276144

  \- smarts: '\[H\]\[CX4H2\](\[OX2H1\])'   \# H on CH2 group: HC (opls\_156)  
    type\_name: 'opls\_156'  
    charge: 0.060  
    sigma: 0.250  
    epsilon: 0.12552

  \- smarts: '\[H\]\[CH3X4\]'           \# H on CH3 group: HC (opls\_156)  
    type\_name: 'opls\_156'  
    charge: 0.060  
    sigma: 0.250  
    epsilon: 0.12552

\# Bond parameters \[63, 70\]  
\# Format: 'TYPE1-TYPE2': \[k\_b (kJ/mol/nm^2), b0 (nm)\]  
bond\_types:  
  'opls\_157-opls\_157': \[224262.4, 0.1529\]  \# CT-CT from \[63\]  
  'opls\_157-opls\_154': \[267776.0, 0.1410\]  \# CT-OH  
  'opls\_154-opls\_155': \[462750.4, 0.0945\]  \# OH-HO  
  'opls\_157-opls\_156': \[284512.0, 0.1090\]  \# CT-HC

\# Angle parameters \[64, 70\]  
\# Format: 'TYPE1-CENTER-TYPE3': \[k\_theta (kJ/mol/rad^2), theta0 (radians)\]  
\# 109.5 deg \= 1.911 rad; 107.5 deg \= 1.876 rad; 108.5 deg \= 1.894 rad  
angle\_types:  
  'opls\_157-opls\_157-opls\_154': \[418.4, 1.911\]  \# CT-CT-OH  
  'opls\_156-opls\_157-opls\_157': \[292.8, 1.911\]  \# HC-CT-CT  
  'opls\_156-opls\_157-opls\_154': \[292.8, 1.911\]  \# HC-CT-OH  
  'opls\_156-opls\_157-opls\_156': \[276.1, 1.876\]  \# HC-CT-HC  
  'opls\_157-opls\_154-opls\_155': \[334.7, 1.894\]  \# CT-OH-HO

\# Dihedral parameters \[64, 70, 71, 72\]  
\# OPLS Fourier series (function type 3\)   
\# Format: 'T1-T2-T3-T4': \[V1, V2, V3, V4\] (all in kJ/mol)  
dihedral\_types:  
  'opls\_156-opls\_157-opls\_157-opls\_154': \[0.97905, 2.93716, 0.0, \-3.91622\] \# HC-CT-CT-OH \[64\]  
  'opls\_156-opls\_157-opls\_157-opls\_156': \[0.62760, 1.88280, 0.0, \-2.51040\] \# HC-CT-CT-HC \[64\]  
  'opls\_157-opls\_157-opls\_154-opls\_155': \[\-0.44350, 3.83255, 0.72801, \-4.11705\] \# CT-CT-OH-HO \[64\]  
  'opls\_156-opls\_157-opls\_154-opls\_155': \[0.94140, 2.82420, 0.0, \-3.76560\] \# HC-CT-OH-HO \[64\]

\# \--- End of ethanol.yaml \---

Note on OPLS Dihedral Formula: The calculator in Module 4 must use the OPLS-style Fourier series (function type 3\) corresponding to these $V\_1$ to $V\_4$ parameters 73:

$$V(\\phi) \= \\frac{V\_1}{2}(1+\\cos(\\phi)) \+ \\frac{V\_2}{2}(1-\\cos(2\\phi)) \+ \\frac{V\_3}{2}(1+\\cos(3\\phi)) \+ \\frac{V\_4}{2}(1-\\cos(4\\phi))$$

### **3.2 Module 3: The Parameter Assignment Engine (Member 2\)**

This function implements the "atom typing" logic using the SMARTS rules from the .yaml file and then builds parameter maps for the calculator module.

Python

def assign\_parameters(molecule, ff\_parameters):  
    """  
    Assigns atom types and parameters to the Molecule object.  
    Returns a dictionary of parameter maps for the calculator.  
    """  
    rdkit\_mol \= molecule.rdkit\_mol  \# Get the Mol object from Module 2  
    atom\_type\_rules \= ff\_parameters\['atom\_types'\]  
      
    \# \--- 3.2.1: Atom Typing and Parameter Assignment \---  
    for atom in molecule.atoms:  
        for rule in atom\_type\_rules:  
            smarts\_pattern \= Chem.MolFromSmarts(rule\['smarts'\])  
              
            \# Find all atoms matching this SMARTS pattern \[44\]  
            matches \= rdkit\_mol.GetSubstructMatches(smarts\_pattern)  
              
            \# A match is a tuple of atom indices (e.g., (4, 3, 0)).  
            \# The atom type applies to the \*first\* atom in the SMARTS  
            \# (index 0). We check if our atom.index is in that list.  
            atom\_index\_in\_matches \= \[m for m in matches\]  
              
            if atom.index in atom\_index\_in\_matches:  
                \# Found the first matching (most specific) rule  
                atom.atom\_type \= rule\['type\_name'\]  
                atom.charge \= rule\['charge'\]  
                atom.sigma \= rule\['sigma'\]  
                atom.epsilon \= rule\['epsilon'\]  
                break  \# Stop processing rules for this atom  
          
        if atom.atom\_type is None:  
            raise Exception(f"Atom {atom.index} ({atom.element}) did not match any atom\_type rules.")  
              
    \# \--- 3.2.2: Build Parameter Maps for Calculator \---  
    \# These maps link atom indices to their specific parameters,  
    \# making the calculator (Module 4\) simple and fast.  
    param\_maps \= {  
        'bonds': {},  
        'angles': {},  
        'dihedrals': {}  
    }

    \# Map bonds: (i, j) \-\> \[kb, b0\]  
    for i, j in molecule.bonds:  
        type\_i \= molecule.atoms\[i\].atom\_type  
        type\_j \= molecule.atoms\[j\].atom\_type  
        \# Create a sorted, unique key (e.g., 'CT-HC')   
        key \= "-".join(sorted(\[type\_i, type\_j\]))  
        if key in ff\_parameters\['bond\_types'\]:  
            param\_maps\['bonds'\]\[(i, j)\] \= ff\_parameters\['bond\_types'\]\[key\]  
        else:  
            raise Exception(f"Missing bond parameter for type {key}")

    \# Map angles: (i, j\_center, k) \-\> \[k\_theta, theta0\]  
    for i, j\_center, k in molecule.angles:  
        type\_i \= molecule.atoms\[i\].atom\_type  
        type\_j \= molecule.atoms\[j\_center\].atom\_type  
        type\_k \= molecule.atoms\[k\].atom\_type  
        \# Key sorted on outer atoms (e.g., 'HC-CT-OH')  
        key\_atoms \= sorted(\[type\_i, type\_k\])  
        key \= f"{key\_atoms}\-{type\_j}\-{key\_atoms}"  
          
        if key in ff\_parameters\['angle\_types'\]:  
            param\_maps\['angles'\]\[(i, j\_center, k)\] \= ff\_parameters\['angle\_types'\]\[key\]  
        else:  
            raise Exception(f"Missing angle parameter for type {key}")

    \# Map dihedrals: (i, j, k, l) \-\> \[V1, V2, V3, V4\]  
    for i, j, k, l in molecule.dihedrals:  
        type\_i \= molecule.atoms\[i\].atom\_type  
        type\_j \= molecule.atoms\[j\].atom\_type  
        type\_k \= molecule.atoms\[k\].atom\_type  
        type\_l \= molecule.atoms\[l\].atom\_type  
        key \= f"{type\_i}\-{type\_j}\-{type\_k}\-{type\_l}"  
          
        if key in ff\_parameters\['dihedral\_types'\]:  
            param\_maps\['dihedrals'\]\[(i, j, k, l)\] \= ff\_parameters\['dihedral\_types'\]\[key\]  
        else:  
            \# Check for reversed key (l-k-j-i)  
            key\_rev \= f"{type\_l}\-{type\_k}\-{type\_j}\-{type\_i}"  
            if key\_rev in ff\_parameters\['dihedral\_types'\]:  
                param\_maps\['dihedrals'\]\[(i, j, k, l)\] \= ff\_parameters\['dihedral\_types'\]\[key\_rev\]  
            else:  
                \# OPLS allows wildcards 'X' \[72\], but for this  
                \# project, we will assume all specific types are defined.  
                pass \# Or raise Exception(f"Missing dihedral {key}")  
                  
    return param\_maps  
    \`\`\`

\---

\#\# 4\. Module 4: The Core Energy Calculator (HPC-1) (Member 3\)

This module implements the mathematical formulas from Section 1.1. It is the responsibility of \*\*Member 3\*\* to write these functions, including both the naïve and optimized non-bonded calculators for their final report.

\#\#\# 4.1 Geometry and Bonded Energy Functions

These functions consume the \`Molecule\` object and the \`param\_maps\` dictionary.

\`\`\`python  
\# \--- Helper functions for geometry calculations \---

def get\_distance(coords, i, j):  
    """Calculate distance between atom i and atom j."""  
    return np.linalg.norm(coords\[i\] \- coords\[j\])

def get\_angle(coords, i, j, k):  
    """Calculate angle i-j-k in radians."""  
    v\_ji \= coords\[i\] \- coords\[j\]  
    v\_jk \= coords\[k\] \- coords\[j\]  
    cosine\_angle \= np.dot(v\_ji, v\_jk) / (np.linalg.norm(v\_ji) \* np.linalg.norm(v\_jk))  
    return np.arccos(np.clip(cosine\_angle, \-1.0, 1.0))

def get\_dihedral(coords, i, j, k, l):  
    """Calculate dihedral i-j-k-l in radians."""  
    b1 \= coords\[j\] \- coords\[i\]  
    b2 \= coords\[k\] \- coords\[j\]  
    b3 \= coords\[l\] \- coords\[k\]  
      
    n1 \= np.cross(b1, b2)  
    n1 /= np.linalg.norm(n1)  
      
    n2 \= np.cross(b2, b3)  
    n2 /= np.linalg.norm(n2)  
      
    m1 \= np.cross(n1, b2 / np.linalg.norm(b2))  
      
    x \= np.dot(n1, n2)  
    y \= np.dot(m1, n2)  
      
    return np.arctan2(y, x)

\# \--- Energy calculation functions \[45\] \---

def calculate\_bond\_energy(molecule, param\_maps):  
    energy \= 0.0  
    coords \= molecule.coordinates  
    for (i, j), params in param\_maps\['bonds'\].items():  
        kb, b0 \= params  
        r \= get\_distance(coords, i, j)  
        energy \+= 0.5 \* kb \* (r \- b0)\*\*2  \# \[3, 4\]  
    return energy

def calculate\_angle\_energy(molecule, param\_maps):  
    energy \= 0.0  
    coords \= molecule.coordinates  
    for (i, j\_center, k), params in param\_maps\['angles'\].items():  
        k\_theta, theta0 \= params  
        theta \= get\_angle(coords, i, j\_center, k)  
        energy \+= 0.5 \* k\_theta \* (theta \- theta0)\*\*2  \#   
    return energy

def calculate\_dihedral\_energy(molecule, param\_maps):  
    energy \= 0.0  
    coords \= molecule.coordinates  
    for (i, j, k, l), params in param\_maps\['dihedrals'\].items():  
        V1, V2, V3, V4 \= params  
        phi \= get\_dihedral(coords, i, j, k, l)  
          
        \# OPLS-specific dihedral formula   
        energy\_term \= (  
            0.5 \* V1 \* (1.0 \+ np.cos(phi)) \+  
            0.5 \* V2 \* (1.0 \- np.cos(2.0 \* phi)) \+  
            0.5 \* V3 \* (1.0 \+ np.cos(3.0 \* phi)) \+  
            0.5 \* V4 \* (1.0 \- np.cos(4.0 \* phi))  
        )  
        energy \+= energy\_term  
    return energy

### **4.2 The $O(N^2)$ Bottleneck: Naïve Non-Bonded Calculation**

This function is the "slow" baseline, crucial for demonstrating the HPC-1 optimization. It iterates over all possible pairs of atoms.10

Python

def calculate\_nonbonded\_naive(molecule):  
    energy \= 0.0  
    coords \= molecule.coordinates  
    atoms \= molecule.atoms  
    num\_atoms \= len(atoms)  
      
    \# Coulomb constant for kJ/mol, nm, e units  
    COULOMB\_CONST \= 1389.35458   
      
    for i in range(num\_atoms):  
        for j in range(i \+ 1, num\_atoms):  
            \# Check if pair is a 1-2, 1-3, or 1-4 exclusion  
            if (i, j) in molecule.non\_bonded\_exclusions:  
                continue  
                  
            atom\_i \= atoms\[i\]  
            atom\_j \= atoms\[j\]  
            r\_ij \= get\_distance(coords, i, j)  
              
            if r\_ij \< 1e-6: continue \# Avoid division by zero  
              
            \# 1\. Lennard-Jones (VDW) Energy   
            \# Use Lorentz-Berthelot combining rules   
            sigma\_ij \= (atom\_i.sigma \+ atom\_j.sigma) / 2.0  
            epsilon\_ij \= np.sqrt(atom\_i.epsilon \* atom\_j.epsilon)  
              
            r\_ratio \= sigma\_ij / r\_ij  
            r6 \= r\_ratio\*\*6  
            r12 \= r6\*\*2  
            energy \+= 4.0 \* epsilon\_ij \* (r12 \- r6)  
              
            \# 2\. Coulomb (Electrostatic) Energy \[3, 17\]  
            energy \+= COULOMB\_CONST \* (atom\_i.charge \* atom\_j.charge) / r\_ij  
              
    return energy

### **4.3 HPC-1: Optimizing Non-Bonded with scipy.cKDTree**

This is the "Verlet-Cell List" task. Non-bonded interactions decay rapidly, so they can be ignored beyond a cutoff radius (e.g., 1.0 nm).33 Instead of a naïve $O(N^2)$ loop, a k-d tree spatial data structure can find all pairs within this cutoff in $O(N \\log N)$ or $O(N)$ time.32 scipy.spatial.cKDTree is a highly optimized C implementation of this algorithm.77

* **Installation:** pip install scipy  
* **Core Logic:**  
  Python  
  from scipy.spatial import cKDTree

  def calculate\_nonbonded\_optimized(molecule, cutoff=1.0):  
      """  
      Calculates non-bonded energy using an optimized k-d tree  
      neighbor search, replacing the O(N^2) loop.  
      """  
      \# cutoff \= 1.0 nm (10 Angstroms)  
      energy \= 0.0  
      coords \= molecule.coordinates  
      atoms \= molecule.atoms  
      COULOMB\_CONST \= 1389.35458 

      \# 1\. Build the k-d tree from coordinates \[32, 78\]  
      tree \= cKDTree(coords)

      \# 2\. Find all unique pairs (i, j) within the cutoff radius \[46, 79\]  
      \# This is the O(N) replacement for the O(N^2) double loop.  
      pairs \= tree.query\_pairs(r=cutoff, output\_type='set')

      for (i, j) in pairs:  
          \# 3\. Check for 1-2, 1-3, 1-4 exclusions  
          if tuple(sorted((i, j))) in molecule.non\_bonded\_exclusions:  
              continue

          atom\_i \= atoms\[i\]  
          atom\_j \= atoms\[j\]

          \# 4\. Calculate energy ONLY for these nearby pairs  
          r\_ij \= get\_distance(coords, i, j) \# We know r\_ij \<= cutoff  
          if r\_ij \< 1e-6: continue

          \# VDW  
          sigma\_ij \= (atom\_i.sigma \+ atom\_j.sigma) / 2.0  
          epsilon\_ij \= np.sqrt(atom\_i.epsilon \* atom\_j.epsilon)  
          r\_ratio \= sigma\_ij / r\_ij  
          r6 \= r\_ratio\*\*6  
          r12 \= r6\*\*2  
          energy \+= 4.0 \* epsilon\_ij \* (r12 \- r6)

          \# Coulomb  
          energy \+= COULOMB\_CONST \* (atom\_i.charge \* atom\_j.charge) / r\_ij

      return energy

**Reportable Outcome (Member 3):** A benchmark plot showing execution time of \_naive vs. \_optimized as a function of system size (e.g., $N$ molecules in a box) will clearly demonstrate the $O(N^2)$ vs. $O(N)$ scaling, which is a classic HPC result.

---

## **5\. Module 5: HPC Scaling (HPC-2) & GUI (Member 4\)**

This module addresses the final project task: scaling the calculation to multiple cores for a *batch* of many molecules.

### **5.1 The "Pleasantly Parallel" Problem (HPC-2)**

The task of calculating the energy for 1000 different molecules (or 1000 different conformations of the same molecule) is "pleasantly parallel".34 The calculation for molecule 1 is completely independent of the calculation for molecule 2\.

This is a **CPU-bound** task, dominated by heavy computation, not waiting for I/O.81 In Python, CPU-bound tasks *cannot* be parallelized with threading due to the Global Interpreter Lock (GIL), which ensures only one thread executes Python code at a time.34

The correct solution is the multiprocessing library, which bypasses the GIL by creating separate processes (each with its own GIL), allowing the operating system to schedule them on different CPU cores for true parallel execution.35

### **5.2 Implementation with multiprocessing.Pool**

The multiprocessing.Pool class is the standard tool for this task.35 It manages a pool of worker processes and its .map() function automatically distributes a list of tasks among them.84

* **Core Logic:**  
  Python  
  import multiprocessing as mp  
  import time  
  import os  
  from functools import partial

  \# \--- This MASTER function wraps the entire pipeline for one molecule \---  
  \# It MUST be defined at the top level of the script (not inside  
  \# another function) so it can be "pickled" and sent to workers.\[86\]

  def calculate\_single\_molecule\_energy(xyz\_file\_path, ff\_yaml\_path):  
      """  
      Runs the full pipeline (Modules 1-4) for a single xyz file.  
      This is the function that each worker process will execute.  
      """  
      try:  
          \# Module 1: Load inputs  
          molecule \= load\_xyz(xyz\_file\_path)  
          ff\_params \= load\_force\_field(ff\_yaml\_path)

          \# Module 2: Infer topology  
          infer\_topology(molecule)

          \# Module 3: Assign parameters  
          param\_maps \= assign\_parameters(molecule, ff\_params)

          \# Module 4: Calculate Energy (using the optimized version)  
          energy\_bond \= calculate\_bond\_energy(molecule, param\_maps)  
          energy\_angle \= calculate\_angle\_energy(molecule, param\_maps)  
          energy\_dihedral \= calculate\_dihedral\_energy(molecule, param\_maps)  
          energy\_nonbonded \= calculate\_nonbonded\_optimized(molecule)

          total\_energy \= (energy\_bond \+ energy\_angle \+   
                          energy\_dihedral \+ energy\_nonbonded)

          return (xyz\_file\_path, total\_energy)

      except Exception as e:  
          return (xyz\_file\_path, f"Error: {e}")  
  \# \------------------------------------------------------------

  def run\_parallel\_calculations(list\_of\_xyz\_files, ff\_yaml\_path):  
      """  
      Distributes the energy calculation for a list of molecules  
      across all available CPU cores.  
      """

      \# 1\. Use functools.partial to "freeze" the force field argument.  
      \# pool.map only accepts a function with one iterable argument.  
      worker\_function \= partial(calculate\_single\_molecule\_energy,   
                                ff\_yaml\_path=ff\_yaml\_path)

      \# 2\. Get number of available CPU cores \[47\]  
      n\_cores \= mp.cpu\_count()  
      print(f"Starting parallel calculation on {n\_cores} cores...")

      \# 3\. Create the Pool and run the calculations \[35, 87\]  
      start\_time \= time.time()  
      with mp.Pool(processes=n\_cores) as pool:  
          \# map() distributes the tasks and blocks until all are complete  
          results \= pool.map(worker\_function, list\_of\_xyz\_files)

      end\_time \= time.time()  
      print(f"--- Parallel execution finished in {end\_time \- start\_time:.4f} seconds \---")

      return results

  \# \--- Example main execution block \---  
  if \_\_name\_\_ \== "\_\_main\_\_":  
      \# This block is ESSENTIAL for multiprocessing to work  
      \# correctly on all platforms.\[35\]

      \# 1\. Create a dummy list of tasks (e.g., 1000 molecules)  
      \# We will use the \*same\* ethanol file 1000 times for this benchmark.  
      \# (Assumes 'ethanol.xyz' and 'ethanol.yaml' are in the same directory)  
      dummy\_tasks \= \['ethanol.xyz'\] \* 1000  
      ff\_file \= 'ethanol.yaml'

      \# 2\. Run in parallel  
      parallel\_results \= run\_parallel\_calculations(dummy\_tasks, ff\_file)

      \# 3\. Run in serial for comparison (for the report)  
      print("Starting serial calculation...")  
      serial\_worker \= partial(calculate\_single\_molecule\_energy, ff\_yaml\_path=ff\_file)  
      start\_time \= time.time()  
      serial\_results \= \[serial\_worker(task) for task in dummy\_tasks\]  
      end\_time \= time.time()  
      print(f"--- Serial execution finished in {end\_time \- start\_time:.4f} seconds \---")

**Reportable Outcome (Member 4):** A benchmark plot showing Time vs. Number of Cores (from 1 to n\_cores) for a fixed job size (e.g., 1000 molecules). This will demonstrate the linear speedup provided by the parallel implementation, a key success metric for HPC.

### **5.3 GUI: A Simple Streamlit Dashboard (Member 4\)**

While PyQt was suggested, it is a complex, heavyweight framework for desktop applications.88 For a simple dashboard (upload files, see a result), a modern web-based framework like **Streamlit** is vastly simpler and faster to implement.50 It allows the creation of an interactive UI with pure Python.36

* **Installation:** pip install streamlit  
* **Implementation (app.py):**  
  Python  
  import streamlit as st  
  import os

  \# \--- Import all your functions from other files \---  
  \# (This assumes all functions are in 'calculator.py')  
  \# from calculator import calculate\_single\_molecule\_energy

  st.title("High-Performance Molecular Energy Calculator")

  \# 1\. File Uploaders \[49\]  
  uploaded\_xyz \= st.file\_uploader("Upload your Molecule (.xyz file)")  
  uploaded\_yaml \= st.file\_uploader("Upload your Force Field (.yaml file)")

  if uploaded\_xyz is not None and uploaded\_yaml is not None:

      \# Must save files to disk for the main script to read  
      with open("temp.xyz", "wb") as f:  
          f.write(uploaded\_xyz.getvalue())  
      with open("temp.yaml", "wb") as f:  
          f.write(uploaded\_yaml.getvalue())

      \# 2\. Run Calculation on Button Click \[91\]  
      if st.button("Calculate Total Energy"):  
          with st.spinner("Calculating..."):

              \# This calls the master function  
              \# (Assumes function is imported)  
              \# result \= calculate\_single\_molecule\_energy("temp.xyz", "temp.yaml")

              \# \--- Placeholder for testing GUI \---  
              import time  
              time.sleep(1) \# Simulate work  
              result \= ("temp.xyz", \-123.4567)  
              \# \--- End Placeholder \---

              \# 3\. Display Results  
              st.success("Calculation complete\!")  
              st.metric(label="Total Potential Energy",   
                        value=f"{result:.4f} kJ/mol")

          \# Cleanup temp files  
          os.remove("temp.xyz")  
          os.remove("temp.yaml")

* **To Run:** Save the file as app.py and execute streamlit run app.py in the terminal. A web browser will automatically open with the interactive dashboard. This provides a professional and modern front-end with minimal development effort.

#### **Works cited**

1. Molecular Mechanics Theory in Brief \- GMU, accessed on November 7, 2025, [https://mason.gmu.edu/\~sslayden/Chem350/manual/docs/MM.pdf](https://mason.gmu.edu/~sslayden/Chem350/manual/docs/MM.pdf)  
2. Molecular Mechanics \- PMC \- PubMed Central \- NIH, accessed on November 7, 2025, [https://pmc.ncbi.nlm.nih.gov/articles/PMC4026342/](https://pmc.ncbi.nlm.nih.gov/articles/PMC4026342/)  
3. Bonded interactions \- Harmonic potential \- GROMACS documentation, accessed on November 7, 2025, [https://manual.gromacs.org/documentation/2019-current/reference-manual/functions/bonded-interactions.html](https://manual.gromacs.org/documentation/2019-current/reference-manual/functions/bonded-interactions.html)  
4. 3.1: Potential Energy Surface and Bonding Interactions \- Chemistry LibreTexts, accessed on November 7, 2025, [https://chem.libretexts.org/Courses/Western\_Washington\_University/Biophysical\_Chemistry\_(Smirnov\_and\_McCarty)/03%3A\_Molecular\_Mechanics\_and\_Statistical\_Thermodynamics/3.01%3A\_Potential\_Energy\_Surface\_and\_Bonding\_Interactions](https://chem.libretexts.org/Courses/Western_Washington_University/Biophysical_Chemistry_\(Smirnov_and_McCarty\)/03%3A_Molecular_Mechanics_and_Statistical_Thermodynamics/3.01%3A_Potential_Energy_Surface_and_Bonding_Interactions)  
5. Potential energy functions, accessed on November 7, 2025, [https://www.ks.uiuc.edu/Research/namd/2.9/ug/node22.html](https://www.ks.uiuc.edu/Research/namd/2.9/ug/node22.html)  
6. Optimal Solution to the Torsional Coefficient Fitting Problem in Force Field Parametrization, accessed on November 7, 2025, [https://pmc.ncbi.nlm.nih.gov/articles/PMC8041298/](https://pmc.ncbi.nlm.nih.gov/articles/PMC8041298/)  
7. Building a Torsional Potential between Thiophene Rings to Illustrate the Basics of Molecular Modeling | Journal of Chemical Education \- ACS Publications, accessed on November 7, 2025, [https://pubs.acs.org/doi/10.1021/acs.jchemed.2c00733](https://pubs.acs.org/doi/10.1021/acs.jchemed.2c00733)  
8. Non-bonded interactions \- GROMACS 2025.3 documentation, accessed on November 7, 2025, [https://manual.gromacs.org/current/reference-manual/functions/nonbonded-interactions.html](https://manual.gromacs.org/current/reference-manual/functions/nonbonded-interactions.html)  
9. Efficient Nonbonded Interactions for Molecular Dynamics on a Graphics Processing Unit, accessed on November 7, 2025, [https://pmc.ncbi.nlm.nih.gov/articles/PMC2841709/](https://pmc.ncbi.nlm.nih.gov/articles/PMC2841709/)  
10. A double exponential potential for van der Waals interaction \- PMC \- NIH, accessed on November 7, 2025, [https://pmc.ncbi.nlm.nih.gov/articles/PMC6555761/](https://pmc.ncbi.nlm.nih.gov/articles/PMC6555761/)  
11. accessed on November 7, 2025, [http://lampz.tugraz.at/\~hadley/ss1/molecules/VdW.php\#:\~:text=The%20Van%20der%20Waal%20bond,is%20typically%20a%20small%20number.](http://lampz.tugraz.at/~hadley/ss1/molecules/VdW.php#:~:text=The%20Van%20der%20Waal%20bond,is%20typically%20a%20small%20number.)  
12. Lennard-Jones potential \- Wikipedia, accessed on November 7, 2025, [https://en.wikipedia.org/wiki/Lennard-Jones\_potential](https://en.wikipedia.org/wiki/Lennard-Jones_potential)  
13. Review of force fields and intermolecular potentials used in atomistic computational materials research \- AIP Publishing, accessed on November 7, 2025, [https://pubs.aip.org/aip/apr/article/5/3/031104/123935/Review-of-force-fields-and-intermolecular](https://pubs.aip.org/aip/apr/article/5/3/031104/123935/Review-of-force-fields-and-intermolecular)  
14. Molecular Dynamics Study of Hydration in Ethanol−Water Mixtures Using a Polarizable Force Field | The Journal of Physical Chemistry B \- ACS Publications, accessed on November 7, 2025, [https://pubs.acs.org/doi/10.1021/jp045438q](https://pubs.acs.org/doi/10.1021/jp045438q)  
15. Screened Electrostatic Interactions in Molecular Mechanics | Journal of Chemical Theory and Computation \- ACS Publications, accessed on November 7, 2025, [https://pubs.acs.org/doi/10.1021/ct5005142](https://pubs.acs.org/doi/10.1021/ct5005142)  
16. accessed on November 7, 2025, [https://en.wikipedia.org/wiki/Coulomb%27s\_law\#:\~:text=The%20law%20states%20that%20the,of%20the%20distance%20between%20them.](https://en.wikipedia.org/wiki/Coulomb%27s_law#:~:text=The%20law%20states%20that%20the,of%20the%20distance%20between%20them.)  
17. Coulomb's law \- Wikipedia, accessed on November 7, 2025, [https://en.wikipedia.org/wiki/Coulomb%27s\_law](https://en.wikipedia.org/wiki/Coulomb%27s_law)  
18. 2.1: Coulomb's Law and the Electrostatic Potential \- Chemistry LibreTexts, accessed on November 7, 2025, [https://chem.libretexts.org/Courses/Oregon\_Institute\_of\_Technology/OIT%3A\_CHE\_202\_-\_General\_Chemistry\_II/Unit\_2%3A\_Electrons\_in\_Atoms/2.1%3A\_Coulomb's\_Law\_and\_the\_Electrostatic\_Potential](https://chem.libretexts.org/Courses/Oregon_Institute_of_Technology/OIT%3A_CHE_202_-_General_Chemistry_II/Unit_2%3A_Electrons_in_Atoms/2.1%3A_Coulomb's_Law_and_the_Electrostatic_Potential)  
19. molgeom \- PyPI, accessed on November 7, 2025, [https://pypi.org/project/molgeom/](https://pypi.org/project/molgeom/)  
20. How to parse (read) .XYZ atomic coordinates file using Python? \[TUTORIAL\] \- BragitOff.com, accessed on November 7, 2025, [https://www.bragitoff.com/2023/07/how-to-parse-read-xyz-atomic-coordinates-file-using-python-tutorial/](https://www.bragitoff.com/2023/07/how-to-parse-read-xyz-atomic-coordinates-file-using-python-tutorial/)  
21. How to input 3D coordinates from xyz file and connectivity from SMILES in rdkit?, accessed on November 7, 2025, [https://mattermodeling.stackexchange.com/questions/7234/how-to-input-3d-coordinates-from-xyz-file-and-connectivity-from-smiles-in-rdkit](https://mattermodeling.stackexchange.com/questions/7234/how-to-input-3d-coordinates-from-xyz-file-and-connectivity-from-smiles-in-rdkit)  
22. Extended XYZ cartesian coordinates format (exyz) \- Open Babel \- Read the Docs, accessed on November 7, 2025, [https://open-babel.readthedocs.io/en/latest/FileFormats/Extended\_XYZ\_cartesian\_coordinates\_format.html](https://open-babel.readthedocs.io/en/latest/FileFormats/Extended_XYZ_cartesian_coordinates_format.html)  
23. How to Work with YAML in Python \- Earthly Blog, accessed on November 7, 2025, [https://earthly.dev/blog/yaml-in-python/](https://earthly.dev/blog/yaml-in-python/)  
24. zotko/xyz2graph: Convert an xyz file into a molecular graph and create a 3D visualisation of the graph. \- GitHub, accessed on November 7, 2025, [https://github.com/zotko/xyz2graph](https://github.com/zotko/xyz2graph)  
25. 01 Molecular Geometry Analysis \- Programming Tutorial in Chemistry by Python, accessed on November 7, 2025, [https://pycrawfordprogproj.readthedocs.io/en/latest/Project\_01/Project\_01.html](https://pycrawfordprogproj.readthedocs.io/en/latest/Project_01/Project_01.html)  
26. jensengroup/xyz2mol: Converts an xyz file to an RDKit mol object \- GitHub, accessed on November 7, 2025, [https://github.com/jensengroup/xyz2mol](https://github.com/jensengroup/xyz2mol)  
27. Atom typing behavior — ForceField 2025.1 documentation \- SCM, accessed on November 7, 2025, [https://www.scm.com/doc/ForceField/Atom\_typing\_behaviour.html](https://www.scm.com/doc/ForceField/Atom_typing_behaviour.html)  
28. MATCH: An Atom- Typing Toolset for Molecular Mechanics Force Fields \- PMC \- NIH, accessed on November 7, 2025, [https://pmc.ncbi.nlm.nih.gov/articles/PMC3228871/](https://pmc.ncbi.nlm.nih.gov/articles/PMC3228871/)  
29. 6\. Creating Force Fields — OpenMM User Guide 8.0 documentation, accessed on November 7, 2025, [https://docs.openmm.org/8.0.0/userguide/application/05\_creating\_ffs.html](https://docs.openmm.org/8.0.0/userguide/application/05_creating_ffs.html)  
30. Tweaking and Inspecting Parameters — openff-interchange documentation, accessed on November 7, 2025, [https://docs.openforcefield.org/projects/interchange/en/stable/using/collections.html](https://docs.openforcefield.org/projects/interchange/en/stable/using/collections.html)  
31. cKDTree — SciPy v1.16.2 Manual, accessed on November 7, 2025, [https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.html](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.html)  
32. Cell lists \- Wikipedia, accessed on November 7, 2025, [https://en.wikipedia.org/wiki/Cell\_lists](https://en.wikipedia.org/wiki/Cell_lists)  
33. Running code in parallel — HPC Python \- Jeff Stafford, accessed on November 7, 2025, [https://jstaf.github.io/hpc-python/parallel/](https://jstaf.github.io/hpc-python/parallel/)  
34. multiprocessing — Process-based parallelism — Python 3.14.0 documentation, accessed on November 7, 2025, [https://docs.python.org/3/library/multiprocessing.html](https://docs.python.org/3/library/multiprocessing.html)  
35. Streamlit • A faster way to build and share data apps, accessed on November 7, 2025, [https://streamlit.io/](https://streamlit.io/)  
36. What is the best way to divide work among developers, accessed on November 7, 2025, [https://softwareengineering.stackexchange.com/questions/132622/what-is-the-best-way-to-divide-work-among-developers](https://softwareengineering.stackexchange.com/questions/132622/what-is-the-best-way-to-divide-work-among-developers)  
37. What are the different methods of dividing and organizing work among project team members? \- Software Engineering Stack Exchange, accessed on November 7, 2025, [https://softwareengineering.stackexchange.com/questions/67054/what-are-the-different-methods-of-dividing-and-organizing-work-among-project-tea](https://softwareengineering.stackexchange.com/questions/67054/what-are-the-different-methods-of-dividing-and-organizing-work-among-project-tea)  
38. Getting Started with the RDKit in Python, accessed on November 7, 2025, [https://rdkit.readthedocs.io/en/latest/GettingStartedInPython.html](https://rdkit.readthedocs.io/en/latest/GettingStartedInPython.html)  
39. I need to find a list of all the possible Bonds, angles between atoms in a molecule from a smiles string (or .xyz file) \- Stack Overflow, accessed on November 7, 2025, [https://stackoverflow.com/questions/58390152/i-need-to-find-a-list-of-all-the-possible-bonds-angles-between-atoms-in-a-molec](https://stackoverflow.com/questions/58390152/i-need-to-find-a-list-of-all-the-possible-bonds-angles-between-atoms-in-a-molec)  
40. Systems Header for YAML Files — YANK 0.21.2 documentation, accessed on November 7, 2025, [http://getyank.org/0.21.2/yamlpages/systems.html](http://getyank.org/0.21.2/yamlpages/systems.html)  
41. YAML Tutorial: Everything You Need to Get Started in Minutes \- CloudBees, accessed on November 7, 2025, [https://www.cloudbees.com/blog/yaml-tutorial-everything-you-need-get-started](https://www.cloudbees.com/blog/yaml-tutorial-everything-you-need-get-started)  
42. 4\. SMARTS \- A Language for Describing Molecular Patterns \- Daylight, accessed on November 7, 2025, [https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html](https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html)  
43. Getting Started with the RDKit in Python, accessed on November 7, 2025, [https://www.rdkit.org/docs/GettingStartedInPython.html](https://www.rdkit.org/docs/GettingStartedInPython.html)  
44. Molecular dynamics \- Python in Chemistry, accessed on November 7, 2025, [https://pythoninchemistry.org/sim\_and\_scat/molecular\_dynamics/intro.html](https://pythoninchemistry.org/sim_and_scat/molecular_dynamics/intro.html)  
45. query\_pairs \- cKDTree \- Numpy and Scipy Documentation, accessed on November 7, 2025, [https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.query\_pairs.html](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.query_pairs.html)  
46. Parallel programming in Python: multiprocessing (part 1\) – PDC Blog \- KTH, accessed on November 7, 2025, [https://www.kth.se/blogs/pdc/2019/02/parallel-programming-in-python-multiprocessing-part-1/](https://www.kth.se/blogs/pdc/2019/02/parallel-programming-in-python-multiprocessing-part-1/)  
47. Visualize your multiprocessing calculations in python with parallelbar module \- Medium, accessed on November 7, 2025, [https://medium.com/pythoneers/visualize-your-multiprocessing-calculations-in-python-with-parallelbar-5395651f35aa](https://medium.com/pythoneers/visualize-your-multiprocessing-calculations-in-python-with-parallelbar-5395651f35aa)  
48. Building a dashboard in Python using Streamlit, accessed on November 7, 2025, [https://blog.streamlit.io/crafting-a-dashboard-app-in-python-using-streamlit/](https://blog.streamlit.io/crafting-a-dashboard-app-in-python-using-streamlit/)  
49. Streamlit vs Gradio: The Ultimate Showdown for Python Dashboards | UI Bakery Blog, accessed on November 7, 2025, [https://uibakery.io/blog/streamlit-vs-gradio](https://uibakery.io/blog/streamlit-vs-gradio)  
50. Gradio vs Streamlit vs Dash vs Flask | Towards Data Science, accessed on November 7, 2025, [https://towardsdatascience.com/gradio-vs-streamlit-vs-dash-vs-flask-d3defb1209a2/](https://towardsdatascience.com/gradio-vs-streamlit-vs-dash-vs-flask-d3defb1209a2/)  
51. How to make code that can read an .xyz file and calculate distances between atoms?, accessed on November 7, 2025, [https://stackoverflow.com/questions/67186114/how-to-make-code-that-can-read-an-xyz-file-and-calculate-distances-between-atom](https://stackoverflow.com/questions/67186114/how-to-make-code-that-can-read-an-xyz-file-and-calculate-distances-between-atom)  
52. Force YAML values to be strings \- python \- Stack Overflow, accessed on November 7, 2025, [https://stackoverflow.com/questions/12012774/force-yaml-values-to-be-strings](https://stackoverflow.com/questions/12012774/force-yaml-values-to-be-strings)  
53. yaml.load, force dict keys to strings \- python \- Stack Overflow, accessed on November 7, 2025, [https://stackoverflow.com/questions/50045617/yaml-load-force-dict-keys-to-strings](https://stackoverflow.com/questions/50045617/yaml-load-force-dict-keys-to-strings)  
54. Infer bond orders and formal charges · Issue \#1828 · openforcefield/openff-toolkit \- GitHub, accessed on November 7, 2025, [https://github.com/openforcefield/openff-toolkit/issues/1828](https://github.com/openforcefield/openff-toolkit/issues/1828)  
55. Generating Conformers from MOL/XYZ Files and Aligning Common Atoms in RDKit \#8062, accessed on November 7, 2025, [https://github.com/rdkit/rdkit/discussions/8062](https://github.com/rdkit/rdkit/discussions/8062)  
56. xyz2mol: converting an xyz file to an RDKit mol object \- Proteins and Wave Functions, accessed on November 7, 2025, [http://proteinsandwavefunctions.blogspot.com/2018/01/xyz2mol-converting-xyz-file-to-rdkit.html](http://proteinsandwavefunctions.blogspot.com/2018/01/xyz2mol-converting-xyz-file-to-rdkit.html)  
57. Creating YAML Mechanism Files from Scratch — Cantera 3.1.0 documentation, accessed on November 7, 2025, [https://cantera.org/stable/userguide/creating-mechanisms.html](https://cantera.org/stable/userguide/creating-mechanisms.html)  
58. Accuracy Test of the OPLS-AA Force Field for Calculating Free Energies of Mixing and Comparison with PAC-MAC \- NIH, accessed on November 7, 2025, [https://pmc.ncbi.nlm.nih.gov/articles/PMC5425945/](https://pmc.ncbi.nlm.nih.gov/articles/PMC5425945/)  
59. Influence of Ethanol Parametrization on Diffusion Coefficients Using OPLS-AA Force Field, accessed on November 7, 2025, [https://www.mdpi.com/1422-0067/24/8/7316](https://www.mdpi.com/1422-0067/24/8/7316)  
60. The different bond parameters obtained from LigParGen and from oplsaa.ff in Gromacs, accessed on November 7, 2025, [https://gromacs.bioexcel.eu/t/the-different-bond-parameters-obtained-from-ligpargen-and-from-oplsaa-ff-in-gromacs/4286](https://gromacs.bioexcel.eu/t/the-different-bond-parameters-obtained-from-ligpargen-and-from-oplsaa-ff-in-gromacs/4286)  
61. \[gmx-users\] N-Acetylglucosamine (NAG) in OPLS-AA \- KTH, accessed on November 7, 2025, [https://mailman-1.sys.kth.se/pipermail/gromacs.org\_gmx-users/2011-November/066071.html](https://mailman-1.sys.kth.se/pipermail/gromacs.org_gmx-users/2011-November/066071.html)  
62. Refinement of the OPLS Force Field for Thermodynamics and Dynamics of Liquid Alkanes \- PMC \- NIH, accessed on November 7, 2025, [https://pmc.ncbi.nlm.nih.gov/articles/PMC9939004/](https://pmc.ncbi.nlm.nih.gov/articles/PMC9939004/)  
63. MCCCS Towhee: OPLS-aa \- SourceForge, accessed on November 7, 2025, [https://towhee.sourceforge.net/forcefields/oplsaa.html](https://towhee.sourceforge.net/forcefields/oplsaa.html)  
64. MODELS/oplsaa.ff/ethanol.itp · main · Claire LOISON / MD\_water\_methanol\_Gromacs · GitLab, accessed on November 7, 2025, [https://cameleon.univ-lyon1.fr/cloison/md\_water\_methanol\_gromacs/-/blob/main/MODELS/oplsaa.ff/ethanol.itp](https://cameleon.univ-lyon1.fr/cloison/md_water_methanol_gromacs/-/blob/main/MODELS/oplsaa.ff/ethanol.itp)  
65. Molecular-dynamics simulations of the ethanol liquid–vapor interface \- American Institute of Physics, accessed on November 7, 2025, [https://pubs.aip.org/aip/jcp/article-pdf/119/23/12569/19266700/12569\_1\_online.pdf](https://pubs.aip.org/aip/jcp/article-pdf/119/23/12569/19266700/12569_1_online.pdf)  
66. Implementing OPLS-AA/M Force Field in Gromacs: Lessons Learnt \- LigParGen Server, accessed on November 7, 2025, [https://traken.chem.yale.edu/ligpargen/oplsaam\_gmx\_tutorial.html](https://traken.chem.yale.edu/ligpargen/oplsaam_gmx_tutorial.html)  
67. NonbondedForce — OpenMM Python API 8.4.0.dev-4768436 documentation, accessed on November 7, 2025, [https://docs.openmm.org/development/api-python/generated/openmm.openmm.NonbondedForce.html](https://docs.openmm.org/development/api-python/generated/openmm.openmm.NonbondedForce.html)  
68. Neighbor List Artifacts in Molecular Dynamics Simulations | Journal of Chemical Theory and Computation \- ACS Publications, accessed on November 7, 2025, [https://pubs.acs.org/doi/10.1021/acs.jctc.3c00777](https://pubs.acs.org/doi/10.1021/acs.jctc.3c00777)  
69. Optimizing Python KD Tree Searches \- Stack Overflow, accessed on November 7, 2025, [https://stackoverflow.com/questions/13079010/optimizing-python-kd-tree-searches](https://stackoverflow.com/questions/13079010/optimizing-python-kd-tree-searches)  
70. Solving embarassingly parallel problems using Python multiprocessing \- Stack Overflow, accessed on November 7, 2025, [https://stackoverflow.com/questions/2359253/solving-embarassingly-parallel-problems-using-python-multiprocessing](https://stackoverflow.com/questions/2359253/solving-embarassingly-parallel-problems-using-python-multiprocessing)  
71. Python: Multiprocessing vs Multithreading \- Capital One, accessed on November 7, 2025, [https://www.capitalone.com/tech/software-engineering/python-multiprocessing-multithreading/](https://www.capitalone.com/tech/software-engineering/python-multiprocessing-multithreading/)  
72. Multiprocessing vs Threading Python \[duplicate\] \- Stack Overflow, accessed on November 7, 2025, [https://stackoverflow.com/questions/3044580/multiprocessing-vs-threading-python](https://stackoverflow.com/questions/3044580/multiprocessing-vs-threading-python)  
73. Should I be using multi-threading or multi-processing? : r/learnpython \- Reddit, accessed on November 7, 2025, [https://www.reddit.com/r/learnpython/comments/1grca6v/should\_i\_be\_using\_multithreading\_or/](https://www.reddit.com/r/learnpython/comments/1grca6v/should_i_be_using_multithreading_or/)  
74. How to use multiprocessing pool.map with multiple arguments \- Codemia, accessed on November 7, 2025, [https://codemia.io/knowledge-hub/path/how\_to\_use\_multiprocessing\_poolmap\_with\_multiple\_arguments](https://codemia.io/knowledge-hub/path/how_to_use_multiprocessing_poolmap_with_multiple_arguments)  
75. Parallelization in Python \- Kanishk Varshney \- Medium, accessed on November 7, 2025, [https://kanishkvarshney.medium.com/parallelization-in-python-32d2add653f](https://kanishkvarshney.medium.com/parallelization-in-python-32d2add653f)  
76. 10 Powerful Python GUI Frameworks for 2024: Simplify Your Desktop App Development, accessed on November 7, 2025, [https://fullscale.io/blog/python-gui-frameworks/](https://fullscale.io/blog/python-gui-frameworks/)  
77. Which Python GUI library should you use in 2025?, accessed on November 7, 2025, [https://www.pythonguis.com/faq/which-python-gui-library/](https://www.pythonguis.com/faq/which-python-gui-library/)  
78. which GUI is good : r/learnpython \- Reddit, accessed on November 7, 2025, [https://www.reddit.com/r/learnpython/comments/1dhuwhd/which\_gui\_is\_good/](https://www.reddit.com/r/learnpython/comments/1dhuwhd/which_gui_is_good/)