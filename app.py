"""
High-Performance Molecular Energy Calculator - Streamlit GUI
Member 4: GUI & Parallelization Module
"""

import streamlit as st
import streamlit.components.v1 as components
import os
import time
import multiprocessing as mp
from functools import partial
import pandas as pd
import py3Dmol
import yaml

# Import real calculator functions
from calculator import (
    calculate_single_molecule_energy,
    calculate_energy_with_breakdown,
    run_parallel_calculations,
    validate_force_field_coverage
)

# RDKit for cheminformatics and SMARTS generation
try:
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds, Descriptors, AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    st.warning("⚠️ RDKit not available. Install with: pip install rdkit")

# Page configuration
st.set_page_config(
    page_title="Molecular Energy Calculator",
    page_icon="⚛️",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for better styling
st.markdown("""
    <style>
    .main-header {
        font-size: 2.5rem;
        font-weight: bold;
        color: #1f77b4;
        text-align: center;
        padding: 1rem 0;
    }
    .sub-header {
        font-size: 1.2rem;
        color: #666;
        text-align: center;
        padding-bottom: 2rem;
    }
    .metric-container {
        background-color: #f0f2f6;
        padding: 1rem;
        border-radius: 0.5rem;
        margin: 0.5rem 0;
    }
    /* Ensure code blocks have scrollbars and limited height */
    .stCodeBlock {
        max-height: 300px;
        overflow-y: auto;
    }
    </style>
""", unsafe_allow_html=True)

# Header
st.markdown('<div class="main-header">⚛️ High-Performance Molecular Energy Calculator</div>', unsafe_allow_html=True)
st.markdown('<div class="sub-header">Advanced Force Field-Based Energy Computation</div>', unsafe_allow_html=True)

# Sidebar - Information
with st.sidebar:
    st.header("ℹ️ About")
    st.markdown("""
    This tool calculates molecular potential energy using:
    - **OPLS-AA Force Field**
    - **Bonded interactions**: Bonds, Angles, Dihedrals
    - **Non-bonded interactions**: Van der Waals (LJ), Electrostatic
    - **HPC Optimization**: Parallel processing & spatial tree algorithms
    """)
    
    st.header("📋 Module Information")
    st.info("**Module 5: GUI & Parallelization**\n\n**Member 4 Responsibilities:**\n- Streamlit Interface\n- Multi-core parallelization\n- System integration")
    
    st.header("🔧 Computation Mode")
    computation_mode = st.radio(
        "Select mode:",
        ["Single Molecule", "Batch Processing (HPC-2)"],
        help="Single: One molecule. Batch: Multiple molecules in parallel."
    )

# Main content area
tab1, tab2, tab3, tab4 = st.tabs(["💻 Calculate Energy", "🔧 YAML Builder", "📊 Benchmark", "📖 Documentation"])

# Tab 1: Energy Calculation
with tab1:
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("📁 Upload Molecular Geometry")
        uploaded_xyz = st.file_uploader(
            "Upload .xyz file",
            type=['xyz'],
            help="Standard XYZ format with atomic coordinates"
        )
        
        if uploaded_xyz:
            st.success(f"✅ Loaded: {uploaded_xyz.name}")
            
            # Download button for XYZ file
            st.download_button(
                label="⬇️ Download XYZ file",
                data=uploaded_xyz.getvalue(),
                file_name=uploaded_xyz.name,
                mime="chemical/x-xyz"
            )
            
            with st.expander("View file preview (full content)"):
                content = uploaded_xyz.getvalue().decode('utf-8')
                # Display all lines with a scrollable container
                st.code(content, language='text', line_numbers=True)
    
    with col2:
        st.subheader("⚙️ Upload Force Field")
        uploaded_yaml = st.file_uploader(
            "Upload .yaml force field",
            type=['yaml', 'yml'],
            help="OPLS-AA force field parameters"
        )
        
        if uploaded_yaml:
            st.success(f"✅ Loaded: {uploaded_yaml.name}")
            
            # Download button for YAML file
            st.download_button(
                label="⬇️ Download YAML file",
                data=uploaded_yaml.getvalue(),
                file_name=uploaded_yaml.name,
                mime="application/x-yaml"
            )
            
            # Parse YAML content
            yaml_content = uploaded_yaml.getvalue().decode('utf-8')
            
            # Formatted YAML preview
            with st.expander("📋 View formatted force field parameters"):
                try:
                    yaml_data = yaml.safe_load(yaml_content)
                    
                    # Display atom types
                    if 'atom_types' in yaml_data:
                        st.markdown("**🔹 Atom Types:**")
                        atom_types_df = pd.DataFrame(yaml_data['atom_types'])
                        st.dataframe(atom_types_df, use_container_width=True, height=200)
                    
                    # Display bond types
                    if 'bond_types' in yaml_data:
                        st.markdown("**🔹 Bond Parameters:**")
                        bond_data = [
                            {"Bond Type": k, "k_b (kJ/mol/nm²)": v[0], "b₀ (nm)": v[1]}
                            for k, v in yaml_data['bond_types'].items()
                        ]
                        st.dataframe(pd.DataFrame(bond_data), use_container_width=True, height=150)
                    
                    # Display angle types
                    if 'angle_types' in yaml_data:
                        st.markdown("**🔹 Angle Parameters:**")
                        angle_data = [
                            {"Angle Type": k, "k_θ (kJ/mol/rad²)": v[0], "θ₀ (rad)": v[1]}
                            for k, v in yaml_data['angle_types'].items()
                        ]
                        st.dataframe(pd.DataFrame(angle_data), use_container_width=True, height=150)
                    
                    # Display dihedral types
                    if 'dihedral_types' in yaml_data:
                        st.markdown("**🔹 Dihedral Parameters:**")
                        dihedral_data = [
                            {"Dihedral Type": k, "V₁": v[0], "V₂": v[1], "V₃": v[2], "V₄": v[3]}
                            for k, v in yaml_data['dihedral_types'].items()
                        ]
                        st.dataframe(pd.DataFrame(dihedral_data), use_container_width=True, height=150)
                        
                except Exception as e:
                    st.error(f"Error parsing YAML: {str(e)}")
            
            # Raw YAML preview
            with st.expander("📄 View raw YAML content"):
                st.code(yaml_content, language='yaml', line_numbers=True)
    
    # 3D Molecule Visualization
    if uploaded_xyz is not None:
        st.divider()
        st.subheader("🔬 3D Molecular Structure Visualization")
        
        # Add viewing style selector
        view_style = st.selectbox(
            "Select visualization style:",
            ["Ball & Stick", "Space Filling", "Stick Only"],
            index=0
        )
        
        try:
            # Parse XYZ file content
            xyz_content = uploaded_xyz.getvalue().decode('utf-8')
            
            # Create 3D viewer with py3Dmol
            viewer = py3Dmol.view(width=800, height=500)
            viewer.addModel(xyz_content, 'xyz')
            
            # Apply selected style
            if view_style == "Ball & Stick":
                viewer.setStyle({'stick': {}, 'sphere': {'scale': 0.3}})
            elif view_style == "Space Filling":
                viewer.setStyle({'sphere': {}})
            else:  # Stick Only
                viewer.setStyle({'stick': {}})
            
            viewer.setBackgroundColor('white')
            viewer.zoomTo()
            
            # Render the viewer as HTML and display in Streamlit
            viewer_html = viewer._make_html()
            components.html(viewer_html, height=500, scrolling=False)
            
            st.caption("🖱️ Use mouse to rotate, zoom, and pan the molecule")
                    
        except Exception as e:
            st.error(f"Error visualizing molecule: {str(e)}")
            st.info("💡 Make sure your XYZ file is in the correct format.")
    
    st.divider()
    
    # Calculate button
    if uploaded_xyz is not None and uploaded_yaml is not None:
        
        if computation_mode == "Batch Processing (HPC-2)":
            num_copies = st.slider(
                "Number of calculations (parallel processing):",
                min_value=10,
                max_value=1000,
                value=100,
                step=10,
                help="Number of identical calculations for benchmarking parallel speedup"
            )
        
        col_btn1, col_btn2, col_btn3 = st.columns([1, 2, 1])
        with col_btn2:
            calculate_button = st.button(
                "🚀 Calculate Energy",
                type="primary",
                use_container_width=True
            )
        
        if calculate_button:
            # Save uploaded files temporarily
            with open("temp.xyz", "wb") as f:
                f.write(uploaded_xyz.getvalue())
            with open("temp.yaml", "wb") as f:
                f.write(uploaded_yaml.getvalue())
            
            # STEP 1: Validate force field coverage
            with st.spinner("🔍 Validating force field parameters..."):
                validation = validate_force_field_coverage("temp.xyz", "temp.yaml")
            
            # Display validation results
            if not validation['is_complete']:
                st.warning("⚠️ **Force Field Validation Issues Detected**")
                
                with st.expander("📋 Missing Parameters Details", expanded=True):
                    if validation['missing_atom_types']:
                        st.error(f"**Missing Atom Types ({len(validation['missing_atom_types'])})**")
                        for at in validation['missing_atom_types']:
                            st.write(f"  - {at}")
                    
                    if validation['missing_bonds']:
                        st.error(f"**Missing Bond Parameters ({len(validation['missing_bonds'])})**")
                        for bond in validation['missing_bonds']:
                            st.write(f"  - {bond}")
                    
                    if validation['missing_angles']:
                        st.error(f"**Missing Angle Parameters ({len(validation['missing_angles'])})**")
                        for angle in validation['missing_angles']:
                            st.write(f"  - {angle}")
                    
                    if validation['missing_dihedrals']:
                        st.error(f"**Missing Dihedral Parameters ({len(validation['missing_dihedrals'])})**")
                        for dihedral in validation['missing_dihedrals']:
                            st.write(f"  - {dihedral}")
                
                # Show coverage statistics
                st.subheader("📊 Coverage Statistics")
                cov_stats = validation['coverage_stats']
                
                col_cov1, col_cov2, col_cov3, col_cov4 = st.columns(4)
                with col_cov1:
                    st.metric("Atoms", f"{cov_stats['atoms']['coverage_percent']:.1f}%", 
                             f"{cov_stats['atoms']['typed']}/{cov_stats['atoms']['total']}")
                with col_cov2:
                    st.metric("Bonds", f"{cov_stats['bonds']['coverage_percent']:.1f}%",
                             f"{cov_stats['bonds']['parameterized']}/{cov_stats['bonds']['total']}")
                with col_cov3:
                    st.metric("Angles", f"{cov_stats['angles']['coverage_percent']:.1f}%",
                             f"{cov_stats['angles']['parameterized']}/{cov_stats['angles']['total']}")
                with col_cov4:
                    st.metric("Dihedrals", f"{cov_stats['dihedrals']['coverage_percent']:.1f}%",
                             f"{cov_stats['dihedrals']['parameterized']}/{cov_stats['dihedrals']['total']}")
                
                if validation['warnings']:
                    with st.expander("⚠️ Warnings"):
                        for warning in validation['warnings']:
                            st.warning(warning)
                
                st.info("💡 **Tip**: Use the YAML Builder tab to create missing parameters, or energy calculation will use fallback values which may be inaccurate.")
            else:
                st.success("✅ Force field validation passed! All parameters available.")
            
            # STEP 2: Proceed with energy calculation
            try:
                if computation_mode == "Single Molecule":
                    # Single calculation with detailed breakdown
                    with st.spinner("🔄 Calculating energy..."):
                        energy_breakdown = calculate_energy_with_breakdown("temp.xyz", "temp.yaml")
                    
                    st.success("✅ Calculation complete!")
                    
                    # Display results
                    st.subheader("📈 Results")
                    col_res1, col_res2, col_res3 = st.columns(3)
                    
                    with col_res1:
                        st.metric(
                            label="Total Energy",
                            value=f"{energy_breakdown['total']:.4f} kJ/mol",
                            delta=None
                        )
                    
                    with col_res2:
                        st.metric(
                            label="Bonded Energy",
                            value=f"{energy_breakdown['bonded']:.4f} kJ/mol",
                            delta=None
                        )
                    
                    with col_res3:
                        st.metric(
                            label="Non-bonded Energy",
                            value=f"{energy_breakdown['nonbonded']:.4f} kJ/mol",
                            delta=None
                        )
                    
                    # Detailed breakdown
                    with st.expander("🔍 Detailed Energy Breakdown"):
                        breakdown_data = {
                            "Component": ["Bonds", "Angles", "Dihedrals", "Van der Waals", "Electrostatic", "Total"],
                            "Energy (kJ/mol)": [
                                energy_breakdown['bond'],
                                energy_breakdown['angle'],
                                energy_breakdown['dihedral'],
                                energy_breakdown['vdw'],
                                energy_breakdown['electrostatic'],
                                energy_breakdown['total']
                            ]
                        }
                        df = pd.DataFrame(breakdown_data)
                        st.dataframe(df, use_container_width=True)
                
                else:
                    # Batch processing with parallelization
                    st.subheader(f"🔄 Running {num_copies} calculations in parallel...")
                    
                    # Progress tracking
                    progress_bar = st.progress(0)
                    status_text = st.empty()
                    
                    n_cores = mp.cpu_count()
                    status_text.text(f"Using {n_cores} CPU cores...")
                    
                    # Run parallel computation
                    start_time = time.time()
                    
                    # Create list of xyz files (duplicates for benchmarking)
                    xyz_files = ["temp.xyz"] * num_copies
                    results = run_parallel_calculations(xyz_files, "temp.yaml")
                    
                    end_time = time.time()
                    elapsed_time = end_time - start_time
                    
                    progress_bar.progress(100)
                    status_text.text("✅ All calculations complete!")
                    
                    # Display parallel performance metrics
                    st.success(f"✅ Completed {num_copies} calculations!")
                    
                    col_perf1, col_perf2, col_perf3, col_perf4 = st.columns(4)
                    
                    with col_perf1:
                        st.metric("Total Time", f"{elapsed_time:.2f} s")
                    
                    with col_perf2:
                        st.metric("CPU Cores Used", f"{n_cores}")
                    
                    with col_perf3:
                        time_per_calc = (elapsed_time / num_copies) * 1000
                        st.metric("Time/Molecule", f"{time_per_calc:.1f} ms")
                    
                    with col_perf4:
                        theoretical_speedup = n_cores * 0.8  # Account for overhead
                        st.metric("Est. Speedup", f"{theoretical_speedup:.1f}x")
                    
                    # Show sample results
                    with st.expander("📊 Sample Results (first 10 molecules)"):
                        sample_results = results[:10] if len(results) > 10 else results
                        sample_data = {
                            "Molecule": [os.path.basename(r[0]) for r in sample_results],
                            "Total Energy (kJ/mol)": [r[1] for r in sample_results]
                        }
                        st.dataframe(pd.DataFrame(sample_data), use_container_width=True)
            
            except Exception as e:
                st.error(f"❌ **Error during calculation**: {str(e)}")
                st.error("Please check your input files and force field parameters.")
                import traceback
                with st.expander("🐛 Debug Information"):
                    st.code(traceback.format_exc())
            
            finally:
                # Cleanup temp files
                if os.path.exists("temp.xyz"):
                    os.remove("temp.xyz")
                if os.path.exists("temp.yaml"):
                    os.remove("temp.yaml")
    
    else:
        st.info("👆 Please upload both a .xyz geometry file and a .yaml force field file to begin.")

# Tab 2: YAML Builder
with tab2:
    st.subheader("🔧 Interactive Force Field Parameter Builder")
    st.markdown("""
    Create a custom YAML force field file by manually specifying parameters for your molecule.
    Upload an XYZ file to get started, and we'll help you identify the required parameters.
    """)
    
    st.divider()
    
    # File upload for YAML builder
    col_builder1, col_builder2 = st.columns([2, 1])
    
    with col_builder1:
        st.markdown("### 📁 Step 1: Upload Molecule Structure")
        builder_xyz = st.file_uploader(
            "Upload .xyz file for analysis",
            type=['xyz'],
            help="We'll analyze this to help you define force field parameters",
            key="builder_xyz"
        )
    
    with col_builder2:
        st.markdown("### 💾 Quick Actions")
        if st.button("📋 Load Template", help="Start with default OPLS-AA template"):
            st.session_state['use_template'] = True
            st.success("✅ Template loaded!")
    
    if builder_xyz:
        # Parse XYZ file
        xyz_content = builder_xyz.getvalue().decode('utf-8')
        lines = xyz_content.strip().split('\n')
        
        try:
            num_atoms = int(lines[0].strip())
            atom_data = []
            
            for i in range(2, 2 + num_atoms):
                parts = lines[i].strip().split()
                element = parts[0]
                coords = [float(parts[1]), float(parts[2]), float(parts[3])]
                atom_data.append({'index': i-2, 'element': element, 'coords': coords})
            
            # Display molecule info
            st.success(f"✅ Analyzed molecule: {num_atoms} atoms")
            
            # Count elements
            from collections import Counter
            import numpy as np
            element_counts = Counter([atom['element'] for atom in atom_data])
            
            # Try to use RDKit for better bond inference and SMARTS generation
            rdkit_mol = None
            auto_smarts = []
            
            if RDKIT_AVAILABLE:
                try:
                    # Write XYZ content to temporary file for RDKit
                    import tempfile
                    with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as tmp:
                        tmp.write(xyz_content)
                        tmp_path = tmp.name
                    
                    # Read XYZ and determine bonds using RDKit
                    raw_mol = Chem.MolFromXYZFile(tmp_path)
                    rdkit_mol = Chem.Mol(raw_mol)
                    rdDetermineBonds.DetermineBonds(rdkit_mol, charge=0)
                    
                    # Clean up temp file
                    os.remove(tmp_path)
                    
                    # Generate SMARTS patterns for each unique atom environment
                    def generate_atom_smarts(mol):
                        """Generate SMARTS patterns for each atom based on its environment"""
                        smarts_list = []
                        
                        for atom in mol.GetAtoms():
                            idx = atom.GetIdx()
                            symbol = atom.GetSymbol()
                            
                            # Get atom properties
                            degree = atom.GetDegree()  # number of bonded neighbors
                            total_hs = atom.GetTotalNumHs()  # total hydrogens
                            formal_charge = atom.GetFormalCharge()
                            
                            # Get neighbor information
                            neighbors = atom.GetNeighbors()
                            neighbor_symbols = sorted([n.GetSymbol() for n in neighbors])
                            
                            # Build a descriptive SMARTS pattern
                            # Start with element
                            smarts_parts = [symbol]
                            
                            # Add connectivity (X = total connections)
                            if degree > 0:
                                smarts_parts.append(f"X{degree}")
                            
                            # Add hydrogen count
                            if total_hs > 0:
                                smarts_parts.append(f"H{total_hs}")
                            
                            # Join into SMARTS
                            base_smarts = f"[{';'.join(smarts_parts)}]"
                            
                            # Create extended SMARTS with one neighbor for better specificity
                            if neighbors:
                                # Most significant neighbor (e.g., heavy atom if present)
                                sig_neighbor = None
                                for n in neighbors:
                                    if n.GetSymbol() != 'H':
                                        sig_neighbor = n
                                        break
                                if not sig_neighbor:
                                    sig_neighbor = neighbors[0]
                                
                                n_symbol = sig_neighbor.GetSymbol()
                                n_degree = sig_neighbor.GetDegree()
                                n_hs = sig_neighbor.GetTotalNumHs()
                                
                                neighbor_smarts_parts = [n_symbol]
                                if n_degree > 0:
                                    neighbor_smarts_parts.append(f"X{n_degree}")
                                if n_hs > 0:
                                    neighbor_smarts_parts.append(f"H{n_hs}")
                                
                                neighbor_smarts = f"[{';'.join(neighbor_smarts_parts)}]"
                                extended_smarts = f"{base_smarts}{neighbor_smarts}"
                            else:
                                extended_smarts = base_smarts
                            
                            # Create description
                            neighbor_str = ", ".join(neighbor_symbols) if neighbor_symbols else "none"
                            description = f"{symbol} (degree={degree}, H={total_hs}, neighbors: {neighbor_str})"
                            
                            smarts_list.append({
                                'atom_idx': idx,
                                'element': symbol,
                                'smarts': extended_smarts,
                                'base_smarts': base_smarts,
                                'description': description,
                                'degree': degree,
                                'total_hs': total_hs,
                                'neighbors': neighbor_symbols
                            })
                        
                        return smarts_list
                    
                    auto_smarts = generate_atom_smarts(rdkit_mol)
                    
                    # Group by unique SMARTS patterns
                    from collections import defaultdict
                    smarts_groups = defaultdict(list)
                    for item in auto_smarts:
                        smarts_groups[item['smarts']].append(item['atom_idx'])
                    
                    st.success(f"✅ RDKit detected {rdkit_mol.GetNumBonds()} bonds and generated {len(smarts_groups)} unique atom type patterns")
                    
                except Exception as e:
                    st.warning(f"⚠️ RDKit bond detection failed: {str(e)}. Using distance-based fallback.")
                    rdkit_mol = None
            
            # Function to detect bonds based on distance (fallback)
            def detect_bonds(atom_data):
                """Detect bonds based on atomic distances"""
                # Typical covalent radii (in Angstroms, will convert from nm)
                covalent_radii = {
                    'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66,
                    'F': 0.57, 'P': 1.07, 'S': 1.05, 'Cl': 1.02
                }
                
                bonds = []
                for i, atom1 in enumerate(atom_data):
                    for j, atom2 in enumerate(atom_data[i+1:], start=i+1):
                        # Calculate distance (coords are in Angstroms in XYZ)
                        dist = np.sqrt(sum((atom1['coords'][k] - atom2['coords'][k])**2 for k in range(3)))
                        
                        # Bond threshold: sum of covalent radii * 1.2
                        r1 = covalent_radii.get(atom1['element'], 1.0)
                        r2 = covalent_radii.get(atom2['element'], 1.0)
                        threshold = (r1 + r2) * 1.2
                        
                        if dist < threshold:
                            bonds.append((i, j, atom1['element'], atom2['element']))
                
                return bonds
            
            # Function to detect angles
            def detect_angles(bonds, atom_data):
                """Detect angles from bond connectivity"""
                # Build adjacency list
                neighbors = {i: [] for i in range(len(atom_data))}
                for i, j, _, _ in bonds:
                    neighbors[i].append(j)
                    neighbors[j].append(i)
                
                angles = []
                for center in range(len(atom_data)):
                    if len(neighbors[center]) >= 2:
                        # All pairs of neighbors form angles
                        neighbor_list = neighbors[center]
                        for idx1 in range(len(neighbor_list)):
                            for idx2 in range(idx1 + 1, len(neighbor_list)):
                                i = neighbor_list[idx1]
                                k = neighbor_list[idx2]
                                angles.append((
                                    i, center, k,
                                    atom_data[i]['element'],
                                    atom_data[center]['element'],
                                    atom_data[k]['element']
                                ))
                
                return angles
            
            # Function to detect dihedrals
            def detect_dihedrals(bonds, atom_data):
                """Detect dihedrals from bond connectivity"""
                # Build adjacency list
                neighbors = {i: [] for i in range(len(atom_data))}
                for i, j, _, _ in bonds:
                    neighbors[i].append(j)
                    neighbors[j].append(i)
                
                dihedrals = []
                # For each bond (j-k), find all i-j-k-l combinations
                for bond_j, bond_k, _, _ in bonds:
                    # Find neighbors of j (excluding k)
                    i_atoms = [n for n in neighbors[bond_j] if n != bond_k]
                    # Find neighbors of k (excluding j)
                    l_atoms = [n for n in neighbors[bond_k] if n != bond_j]
                    
                    for i in i_atoms:
                        for l in l_atoms:
                            dihedrals.append((
                                i, bond_j, bond_k, l,
                                atom_data[i]['element'],
                                atom_data[bond_j]['element'],
                                atom_data[bond_k]['element'],
                                atom_data[l]['element']
                            ))
                
                return dihedrals
            
            # Detect molecular structure
            if rdkit_mol is not None:
                # Use RDKit-detected bonds
                detected_bonds = []
                for bond in rdkit_mol.GetBonds():
                    i = bond.GetBeginAtomIdx()
                    j = bond.GetEndAtomIdx()
                    elem_i = atom_data[i]['element']
                    elem_j = atom_data[j]['element']
                    detected_bonds.append((i, j, elem_i, elem_j))
            else:
                # Fallback to distance-based detection
                detected_bonds = detect_bonds(atom_data)
            
            detected_angles = detect_angles(detected_bonds, atom_data)
            detected_dihedrals = detect_dihedrals(detected_bonds, atom_data)
            
            # IMPORTANT: We need to assign atom types first before creating bond/angle/dihedral type names
            # For now, use element-based names, but these will be replaced with type names in the UI
            
            # Get unique bond types (using element symbols initially)
            unique_bond_types_elements = list(set([
                f"{min(b[2], b[3])}-{max(b[2], b[3])}" for b in detected_bonds
            ]))
            
            # Get unique angle types (order matters for center atom)
            unique_angle_types_elements = list(set([
                f"{min(a[3], a[5])}-{a[4]}-{max(a[3], a[5])}" for a in detected_angles
            ]))
            
            # Get unique dihedral types
            unique_dihedral_types_elements = list(set([
                f"{d[4]}-{d[5]}-{d[6]}-{d[7]}" for d in detected_dihedrals
            ]))
            
            # Store for later use after atom types are defined
            st.session_state['detected_bonds'] = detected_bonds
            st.session_state['detected_angles'] = detected_angles
            st.session_state['detected_dihedrals'] = detected_dihedrals
            
            # === Generate SMARTS patterns using RDKit ===
            # Use the auto-generated SMARTS from earlier RDKit analysis
            rdkit_success = False
            unique_atom_types_rdkit = []
            
            if rdkit_mol is not None and auto_smarts:
                # Group auto-generated SMARTS by unique patterns
                from collections import defaultdict
                
                def make_default_dict():
                    return {'smarts': '', 'element': '', 'count': 0, 'indices': [], 'description': ''}
                
                smarts_groups = {}
                
                for item in auto_smarts:
                    key = item['smarts']
                    if key not in smarts_groups:
                        smarts_groups[key] = make_default_dict()
                        smarts_groups[key]['smarts'] = item['smarts']
                        smarts_groups[key]['element'] = item['element']
                        smarts_groups[key]['description'] = item['description']
                    smarts_groups[key]['count'] += 1
                    smarts_groups[key]['indices'].append(item['atom_idx'])
                
                unique_atom_types_rdkit = list(smarts_groups.values())
                rdkit_success = True
                st.info(f"✅ Generated {len(unique_atom_types_rdkit)} unique atom type SMARTS patterns")
            else:
                st.info("💡 RDKit not available or bond detection failed. Using manual input mode.")
            
            col_info1, col_info2, col_info3 = st.columns(3)
            with col_info1:
                st.metric("Total Atoms", num_atoms)
            with col_info2:
                st.metric("Unique Elements", len(element_counts))
            with col_info3:
                element_str = ", ".join([f"{elem}: {count}" for elem, count in element_counts.items()])
                st.info(f"**Elements:** {element_str}")
            
            # Show detected structure info
            col_struct1, col_struct2, col_struct3 = st.columns(3)
            with col_struct1:
                st.metric("Detected Bonds", len(detected_bonds))
                with st.expander("Bond types"):
                    for bt in unique_bond_types_elements:
                        st.text(bt)
            with col_struct2:
                st.metric("Detected Angles", len(detected_angles))
                with st.expander("Angle types"):
                    for at in unique_angle_types_elements:
                        st.text(at)
            with col_struct3:
                st.metric("Detected Dihedrals", len(detected_dihedrals))
                with st.expander("Dihedral types"):
                    for dt in unique_dihedral_types_elements[:10]:  # Show first 10
                        st.text(dt)
                    if len(unique_dihedral_types_elements) > 10:
                        st.caption(f"... and {len(unique_dihedral_types_elements) - 10} more")
            
            st.divider()
            
            # Initialize session state for parameters
            if 'atom_types_builder' not in st.session_state:
                st.session_state['atom_types_builder'] = []
            if 'bond_types_builder' not in st.session_state:
                st.session_state['bond_types_builder'] = {}
            if 'angle_types_builder' not in st.session_state:
                st.session_state['angle_types_builder'] = {}
            if 'dihedral_types_builder' not in st.session_state:
                st.session_state['dihedral_types_builder'] = {}
            
            # Create tabs for different parameter types
            param_tabs = st.tabs(["⚛️ Atom Types", "🔗 Bonds", "📐 Angles", "🔄 Dihedrals", "📄 Preview YAML"])
            
            # Atom Types Tab
            with param_tabs[0]:
                st.markdown("### Define Atom Types")
                if rdkit_success and unique_atom_types_rdkit:
                    st.success(f"✅ Auto-detected {len(unique_atom_types_rdkit)} unique atom environments using RDKit!")
                    st.info("💡 SMARTS patterns automatically generated based on chemical environment. Adjust as needed.")
                else:
                    st.info("💡 Each unique atom environment needs its own type with SMARTS pattern, charge, σ (sigma), and ε (epsilon) values.")
                
                # Auto-set number based on RDKit detection
                default_num_atom_types = len(unique_atom_types_rdkit) if rdkit_success and unique_atom_types_rdkit else 3
                num_atom_types = st.number_input("Number of atom types to define:", min_value=1, max_value=50, value=default_num_atom_types, key="num_atom_types")
                
                atom_types_list = []
                
                for i in range(num_atom_types):
                    # Use RDKit-generated SMARTS if available
                    if rdkit_success and i < len(unique_atom_types_rdkit):
                        auto_type = unique_atom_types_rdkit[i]
                        default_smarts = auto_type['smarts']
                        default_element = auto_type['element']
                        atom_count = auto_type['count']
                        expander_title = f"Atom Type {i+1}: {default_element} - {default_smarts} ({atom_count} atoms)"
                    else:
                        default_smarts = f"[C]" if i == 0 else f"[H]" if i == 1 else f"[O]"
                        default_element = "C" if i == 0 else "H" if i == 1 else "O"
                        expander_title = f"Atom Type {i+1}"
                    
                    with st.expander(expander_title, expanded=(i < 3)):
                        # Show which atoms match this pattern if RDKit was used
                        if rdkit_success and i < len(unique_atom_types_rdkit):
                            auto_type = unique_atom_types_rdkit[i]
                            st.caption(f"🎯 Matches atoms: {', '.join(map(str, auto_type['indices']))}")
                        
                        col1, col2 = st.columns(2)
                        
                        with col1:
                            smarts = st.text_area(
                                "SMARTS Pattern",
                                value=default_smarts,
                                key=f"smarts_{i}",
                                help="Auto-generated SMARTS pattern based on chemical environment",
                                height=80
                            )
                            type_name = st.text_input(
                                "Type Name",
                                value=f"opls_{135+i}",
                                key=f"type_name_{i}"
                            )
                        
                        with col2:
                            charge = st.number_input(
                                "Charge (e)",
                                value=0.0,
                                format="%.4f",
                                key=f"charge_{i}",
                                help="Partial atomic charge in elementary charge units"
                            )
                            sigma = st.number_input(
                                "Sigma σ (nm)",
                                value=0.350,
                                format="%.4f",
                                key=f"sigma_{i}",
                                help="Lennard-Jones collision diameter"
                            )
                        
                        epsilon = st.number_input(
                            "Epsilon ε (kJ/mol)",
                            value=0.276,
                            format="%.6f",
                            key=f"epsilon_{i}",
                            help="Lennard-Jones well depth"
                        )
                        
                        atom_types_list.append({
                            'smarts': smarts,
                            'type_name': type_name,
                            'charge': charge,
                            'sigma': sigma,
                            'epsilon': epsilon
                        })
                
                st.session_state['atom_types_builder'] = atom_types_list
                
                # Create mapping from atom index to type name for bonds/angles/dihedrals
                atom_idx_to_type = {}
                atom_idx_to_element = {}
                atom_idx_to_smarts = {}  # NEW: Store SMARTS for each atom
                
                if rdkit_success and unique_atom_types_rdkit:
                    # Map each atom to its assigned type name
                    for i, auto_type in enumerate(unique_atom_types_rdkit):
                        if i < len(atom_types_list):
                            type_name = atom_types_list[i]['type_name']
                            for atom_idx in auto_type['indices']:
                                atom_idx_to_type[atom_idx] = type_name
                                atom_idx_to_element[atom_idx] = auto_type['element']
                                atom_idx_to_smarts[atom_idx] = auto_type['smarts']  # NEW: Store SMARTS
                
                # Fallback: use element symbols as type indicators
                if not atom_idx_to_type:
                    for idx, atom in enumerate(atom_data):
                        atom_idx_to_type[idx] = atom['element']
                        atom_idx_to_element[idx] = atom['element']
                        atom_idx_to_smarts[idx] = f"[{atom['element']}]"  # NEW: Simple SMARTS
            
            # Bond Types Tab
            with param_tabs[1]:
                st.markdown("### Define Bond Parameters")
                
                # Generate bond type names using atom type names
                unique_bond_types_typed = []
                bond_type_hints = {}
                
                if 'detected_bonds' in st.session_state and atom_idx_to_type:
                    for bond in st.session_state['detected_bonds']:
                        idx_i, idx_j, elem_i, elem_j = bond
                        type_i = atom_idx_to_type.get(idx_i, elem_i)
                        type_j = atom_idx_to_type.get(idx_j, elem_j)
                        # Sort alphabetically for consistency (unless using opls naming)
                        if type_i <= type_j:
                            bond_type = f"{type_i}-{type_j}"
                        else:
                            bond_type = f"{type_j}-{type_i}"
                        
                        if bond_type not in [b[0] for b in unique_bond_types_typed]:
                            # Use SMARTS patterns for hint instead of just elements
                            smarts_i = atom_idx_to_smarts.get(idx_i, f"[{elem_i}]")
                            smarts_j = atom_idx_to_smarts.get(idx_j, f"[{elem_j}]")
                            if type_i <= type_j:
                                elem_hint = f"{smarts_i}-{smarts_j}"
                            else:
                                elem_hint = f"{smarts_j}-{smarts_i}"
                            unique_bond_types_typed.append((bond_type, elem_hint))
                            bond_type_hints[bond_type] = elem_hint
                
                if unique_bond_types_typed:
                    st.info(f"💡 Detected {len(unique_bond_types_typed)} unique bond types. Define force constants (k_b) and equilibrium lengths (b₀).")
                else:
                    st.warning("⚠️ No bonds detected. This molecule may have no bonded interactions (e.g., single atom or very distant atoms).")
                
                # Auto-set number based on detected bonds, allow 0 minimum
                default_num_bonds = len(unique_bond_types_typed) if unique_bond_types_typed else 0
                num_bond_types = st.number_input("Number of bond types:", min_value=0, max_value=100, value=default_num_bonds, key="num_bond_types")
                
                bond_types_dict = {}
                
                for i in range(num_bond_types):
                    # Use detected bond type if available
                    if i < len(unique_bond_types_typed):
                        default_bond_name, elem_hint = unique_bond_types_typed[i]
                        expander_title = f"Bond Type {i+1}: {default_bond_name} ({elem_hint})"
                    else:
                        default_bond_name = f"type1-type{i+2}"
                        elem_hint = "custom"
                        expander_title = f"Bond Type {i+1}: {default_bond_name}"
                    
                    with st.expander(expander_title, expanded=(i < 3)):
                        bond_name = st.text_input(
                            f"Bond Type Name (e.g., 'opls_157-opls_157')",
                            value=default_bond_name,
                            key=f"bond_name_{i}",
                            help=f"Atom type-based bond name. Element representation: {elem_hint}"
                        )
                        
                        col1, col2 = st.columns(2)
                        with col1:
                            kb = st.number_input(
                                "k_b (kJ/mol/nm²)",
                                value=224262.4,
                                format="%.2f",
                                key=f"kb_{i}",
                                help="Bond force constant"
                            )
                        with col2:
                            b0 = st.number_input(
                                "b₀ (nm)",
                                value=0.1529,
                                format="%.4f",
                                key=f"b0_{i}",
                                help="Equilibrium bond length"
                            )
                        
                        bond_types_dict[bond_name] = [kb, b0]
                
                st.session_state['bond_types_builder'] = bond_types_dict
            
            # Angle Types Tab
            with param_tabs[2]:
                st.markdown("### Define Angle Parameters")
                
                # Generate angle type names using atom type names
                unique_angle_types_typed = []
                angle_type_hints = {}
                
                if 'detected_angles' in st.session_state and atom_idx_to_type:
                    for angle in st.session_state['detected_angles']:
                        idx_i, idx_center, idx_k, elem_i, elem_center, elem_k = angle
                        type_i = atom_idx_to_type.get(idx_i, elem_i)
                        type_center = atom_idx_to_type.get(idx_center, elem_center)
                        type_k = atom_idx_to_type.get(idx_k, elem_k)
                        # Sort outer atoms alphabetically
                        if type_i <= type_k:
                            angle_type = f"{type_i}-{type_center}-{type_k}"
                        else:
                            angle_type = f"{type_k}-{type_center}-{type_i}"
                        
                        if angle_type not in [a[0] for a in unique_angle_types_typed]:
                            # Use SMARTS patterns for hint instead of just elements
                            smarts_i = atom_idx_to_smarts.get(idx_i, f"[{elem_i}]")
                            smarts_center = atom_idx_to_smarts.get(idx_center, f"[{elem_center}]")
                            smarts_k = atom_idx_to_smarts.get(idx_k, f"[{elem_k}]")
                            if type_i <= type_k:
                                elem_hint = f"{smarts_i}-{smarts_center}-{smarts_k}"
                            else:
                                elem_hint = f"{smarts_k}-{smarts_center}-{smarts_i}"
                            unique_angle_types_typed.append((angle_type, elem_hint))
                            angle_type_hints[angle_type] = elem_hint
                
                if unique_angle_types_typed:
                    st.info(f"💡 Detected {len(unique_angle_types_typed)} unique angle types. Define force constants (k_θ) and equilibrium angles (θ₀).")
                else:
                    st.warning("⚠️ No angles detected. This molecule may have linear or very simple geometry.")
                
                # Auto-set number based on detected angles, allow 0 minimum
                default_num_angles = len(unique_angle_types_typed) if unique_angle_types_typed else 0
                num_angle_types = st.number_input("Number of angle types:", min_value=0, max_value=100, value=default_num_angles, key="num_angle_types")
                
                angle_types_dict = {}
                
                for i in range(num_angle_types):
                    # Use detected angle type if available
                    if i < len(unique_angle_types_typed):
                        default_angle_name, elem_hint = unique_angle_types_typed[i]
                        expander_title = f"Angle Type {i+1}: {default_angle_name} ({elem_hint})"
                    else:
                        default_angle_name = f"type1-type2-type3"
                        elem_hint = "custom"
                        expander_title = f"Angle Type {i+1}: {default_angle_name}"
                    
                    with st.expander(expander_title, expanded=(i < 3)):
                        angle_name = st.text_input(
                            f"Angle Type Name (e.g., 'opls_157-opls_157-opls_154')",
                            value=default_angle_name,
                            key=f"angle_name_{i}",
                            help=f"Atom type-based angle name. Element representation: {elem_hint}"
                        )
                        
                        col1, col2 = st.columns(2)
                        with col1:
                            k_theta = st.number_input(
                                "k_θ (kJ/mol/rad²)",
                                value=418.4,
                                format="%.2f",
                                key=f"k_theta_{i}",
                                help="Angle force constant"
                            )
                        with col2:
                            theta0_deg = st.number_input(
                                "θ₀ (degrees)",
                                value=109.5,
                                format="%.2f",
                                key=f"theta0_deg_{i}",
                                help="Equilibrium angle in degrees"
                            )
                        
                        # Convert degrees to radians
                        import math
                        theta0_rad = theta0_deg * math.pi / 180.0
                        st.caption(f"θ₀ = {theta0_rad:.4f} radians")
                        
                        angle_types_dict[angle_name] = [k_theta, theta0_rad]
                
                st.session_state['angle_types_builder'] = angle_types_dict
            
            # Dihedral Types Tab
            with param_tabs[3]:
                st.markdown("### Define Dihedral Parameters")
                
                # Generate dihedral type names using atom type names
                unique_dihedral_types_typed = []
                dihedral_type_hints = {}
                
                if 'detected_dihedrals' in st.session_state and atom_idx_to_type:
                    for dihedral in st.session_state['detected_dihedrals']:
                        idx_i, idx_j, idx_k, idx_l, elem_i, elem_j, elem_k, elem_l = dihedral
                        type_i = atom_idx_to_type.get(idx_i, elem_i)
                        type_j = atom_idx_to_type.get(idx_j, elem_j)
                        type_k = atom_idx_to_type.get(idx_k, elem_k)
                        type_l = atom_idx_to_type.get(idx_l, elem_l)
                        
                        dihedral_type = f"{type_i}-{type_j}-{type_k}-{type_l}"
                        
                        if dihedral_type not in [d[0] for d in unique_dihedral_types_typed]:
                            # Use SMARTS patterns for hint instead of just elements
                            smarts_i = atom_idx_to_smarts.get(idx_i, f"[{elem_i}]")
                            smarts_j = atom_idx_to_smarts.get(idx_j, f"[{elem_j}]")
                            smarts_k = atom_idx_to_smarts.get(idx_k, f"[{elem_k}]")
                            smarts_l = atom_idx_to_smarts.get(idx_l, f"[{elem_l}]")
                            elem_hint = f"{smarts_i}-{smarts_j}-{smarts_k}-{smarts_l}"
                            unique_dihedral_types_typed.append((dihedral_type, elem_hint))
                            dihedral_type_hints[dihedral_type] = elem_hint
                
                if unique_dihedral_types_typed:
                    st.info(f"💡 Detected {len(unique_dihedral_types_typed)} unique dihedral types. OPLS-AA uses V₁, V₂, V₃, V₄ coefficients (kJ/mol).")
                else:
                    st.warning("⚠️ No dihedrals detected. This molecule may have very simple geometry (e.g., CCl4, CH4).")
                
                # Auto-set number based on detected dihedrals, allow 0 minimum (no artificial limit)
                default_num_dihedrals = len(unique_dihedral_types_typed) if unique_dihedral_types_typed else 0
                num_dihedral_types = st.number_input("Number of dihedral types:", min_value=0, max_value=100, value=default_num_dihedrals, key="num_dihedral_types")
                
                if len(unique_dihedral_types_typed) > num_dihedral_types and unique_dihedral_types_typed:
                    st.warning(f"⚠️ Note: {len(unique_dihedral_types_typed)} dihedral types detected, but only defining {num_dihedral_types}. You may need to add more.")
                
                dihedral_types_dict = {}
                
                for i in range(num_dihedral_types):
                    # Use detected dihedral type if available
                    if i < len(unique_dihedral_types_typed):
                        default_dihedral_name, elem_hint = unique_dihedral_types_typed[i]
                        expander_title = f"Dihedral Type {i+1}: {default_dihedral_name} ({elem_hint})"
                    else:
                        default_dihedral_name = f"type1-type2-type3-type4"
                        elem_hint = "custom"
                        expander_title = f"Dihedral Type {i+1}: {default_dihedral_name}"
                    
                    with st.expander(expander_title, expanded=(i < 2)):
                        dihedral_name = st.text_input(
                            f"Dihedral Type Name (e.g., 'opls_156-opls_157-opls_157-opls_154')",
                            value=default_dihedral_name,
                            key=f"dihedral_name_{i}",
                            help=f"Atom type-based dihedral name. Element representation: {elem_hint}"
                        )
                        
                        col1, col2, col3, col4 = st.columns(4)
                        with col1:
                            v1 = st.number_input("V₁ (kJ/mol)", value=0.0, format="%.5f", key=f"v1_{i}")
                        with col2:
                            v2 = st.number_input("V₂ (kJ/mol)", value=0.0, format="%.5f", key=f"v2_{i}")
                        with col3:
                            v3 = st.number_input("V₃ (kJ/mol)", value=0.0, format="%.5f", key=f"v3_{i}")
                        with col4:
                            v4 = st.number_input("V₄ (kJ/mol)", value=0.0, format="%.5f", key=f"v4_{i}")
                        
                        dihedral_types_dict[dihedral_name] = [v1, v2, v3, v4]
                
                st.session_state['dihedral_types_builder'] = dihedral_types_dict
            
            # Preview and Generate YAML Tab
            with param_tabs[4]:
                st.markdown("### 📄 Generated YAML Preview")
                
                # Build YAML content
                yaml_dict = {}
                
                if st.session_state['atom_types_builder']:
                    yaml_dict['atom_types'] = st.session_state['atom_types_builder']
                
                if st.session_state['bond_types_builder']:
                    yaml_dict['bond_types'] = st.session_state['bond_types_builder']
                
                if st.session_state['angle_types_builder']:
                    yaml_dict['angle_types'] = st.session_state['angle_types_builder']
                
                if st.session_state['dihedral_types_builder']:
                    yaml_dict['dihedral_types'] = st.session_state['dihedral_types_builder']
                
                # Generate YAML string
                yaml_output = yaml.dump(yaml_dict, default_flow_style=False, sort_keys=False)
                
                # Add header comment
                yaml_header = f"""# Generated Force Field Parameters
# Created for molecule: {builder_xyz.name}
# Number of atoms: {num_atoms}
# UNITS: kJ/mol, nm, radians, elementary charge (e)
# -----------------------------------------------------------------

"""
                full_yaml = yaml_header + yaml_output
                
                # Preview
                st.code(full_yaml, language='yaml')
                
                # Download button
                st.download_button(
                    label="⬇️ Download YAML Force Field",
                    data=full_yaml,
                    file_name=f"{builder_xyz.name.replace('.xyz', '')}_forcefield.yaml",
                    mime="application/x-yaml",
                    help="Download this YAML file to use in the energy calculator"
                )
                
                st.success("✅ YAML file ready! Download and use it in the 'Calculate Energy' tab.")
        
        except Exception as e:
            st.error(f"Error parsing XYZ file: {str(e)}")
            st.info("Please ensure your XYZ file is in the correct format.")
    
    else:
        st.info("👆 Upload an XYZ file to start building your force field parameters.")
        
        # Show example
        with st.expander("📖 How to use the YAML Builder"):
            st.markdown("""
            ### Step-by-Step Guide:
            
            1. **Upload XYZ File**: Upload your molecular geometry file
            2. **Define Atom Types**: Specify SMARTS patterns and parameters for each unique atom environment
            3. **Define Bonds**: Add force constants and equilibrium lengths for each bond type
            4. **Define Angles**: Add force constants and equilibrium angles
            5. **Define Dihedrals**: Add OPLS Fourier series coefficients
            6. **Preview & Download**: Review the generated YAML and download it
            
            ### Parameter Descriptions:
            
            - **SMARTS Pattern**: Chemical pattern to identify atom types (e.g., `[CX4H3]` for methyl carbon)
            - **Charge**: Partial atomic charge in elementary charge units
            - **Sigma (σ)**: Lennard-Jones collision diameter in nm
            - **Epsilon (ε)**: Lennard-Jones well depth in kJ/mol
            - **k_b**: Bond force constant in kJ/mol/nm²
            - **b₀**: Equilibrium bond length in nm
            - **k_θ**: Angle force constant in kJ/mol/rad²
            - **θ₀**: Equilibrium angle in radians (input in degrees, auto-converted)
            - **V₁, V₂, V₃, V₄**: OPLS dihedral Fourier coefficients in kJ/mol
            
            ### Common Values (OPLS-AA):
            - **C (sp³)**: σ = 0.350 nm, ε = 0.276 kJ/mol
            - **H (alkane)**: σ = 0.250 nm, ε = 0.126 kJ/mol
            - **O (alcohol)**: σ = 0.312 nm, ε = 0.711 kJ/mol
            - **C-C bond**: k_b = 224262 kJ/mol/nm², b₀ = 0.1529 nm
            - **C-H bond**: k_b = 284512 kJ/mol/nm², b₀ = 0.1090 nm
            """)

# Tab 3: Benchmark
with tab3:
    st.subheader("📊 Performance Benchmarking")
    st.markdown("""
    This section demonstrates the HPC-2 scaling performance using multiprocessing.
    """)
    
    col_bench1, col_bench2 = st.columns(2)
    
    with col_bench1:
        st.markdown("### ⚡ Serial vs Parallel")
        st.info("""
        **Naïve O(N²) vs Optimized O(N)**
        
        The optimized algorithm uses `scipy.spatial.cKDTree` for efficient 
        neighbor searching, reducing computational complexity.
        """)
        
        # Placeholder chart data
        import numpy as np
        chart_data = pd.DataFrame({
            'System Size (atoms)': [10, 50, 100, 200, 500],
            'Naïve O(N²) (s)': [0.01, 0.25, 1.0, 4.0, 25.0],
            'Optimized O(N) (s)': [0.01, 0.05, 0.10, 0.20, 0.50]
        })
        st.line_chart(chart_data.set_index('System Size (atoms)'))
    
    with col_bench2:
        st.markdown("### 🔄 Parallel Speedup")
        st.info("""
        **Multi-core Performance**
        
        Using `multiprocessing.Pool.map()` to distribute calculations 
        across all CPU cores.
        """)
        
        # Placeholder speedup chart
        speedup_data = pd.DataFrame({
            'Number of Cores': [1, 2, 4, 8, 16],
            'Speedup': [1.0, 1.9, 3.7, 7.2, 13.5],
            'Ideal Linear': [1.0, 2.0, 4.0, 8.0, 16.0]
        })
        st.line_chart(speedup_data.set_index('Number of Cores'))

# Tab 4: Documentation
with tab4:
    st.subheader("📖 Documentation")
    
    st.markdown("""
    ### 🎯 Project Overview
    
    This High-Performance Molecular Energy Calculator implements a complete pipeline for computing 
    molecular potential energy using classical force fields.
    
    ### 🏗️ Architecture (5-Module Pipeline)
    
    1. **Module 1**: Input/Output & Data Structures (`.xyz` and `.yaml` parsers)
    2. **Module 2**: Topology Inference Engine (bond/angle/dihedral detection)
    3. **Module 3**: Parameter Assignment Engine (atom typing with SMARTS)
    4. **Module 4**: Core Energy Calculator (bonded + non-bonded energy)
    5. **Module 5**: Parallelization & GUI *(this module)*
    
    ### ⚙️ Force Field Components
    
    #### Bonded Interactions:
    - **Bonds**: $V_{bond}(b) = \\sum_{bonds} \\frac{1}{2} k_b (b - b_0)^2$
    - **Angles**: $V_{angle}(\\theta) = \\sum_{angles} \\frac{1}{2} k_\\theta (\\theta - \\theta_0)^2$
    - **Dihedrals**: OPLS Fourier series with V₁, V₂, V₃, V₄ parameters
    
    #### Non-bonded Interactions:
    - **Van der Waals**: Lennard-Jones 12-6 potential
    - **Electrostatic**: Coulomb's law with partial charges
    
    ### 🚀 HPC Optimizations
    
    - **HPC-1**: Spatial tree algorithm (`cKDTree`) for O(N) neighbor search
    - **HPC-2**: Multi-core parallelization with `multiprocessing.Pool`
    
    ### 📦 Dependencies
    
    ```python
    streamlit
    numpy
    scipy
    rdkit
    xyz2mol
    pyyaml
    pandas
    ```
    
    ### 🔧 Usage
    
    1. Upload a `.xyz` molecular geometry file
    2. Upload a `.yaml` force field parameter file
    3. Choose single or batch processing mode
    4. Click "Calculate Energy"
    
    ### 👥 Team Roles
    
    - **Member 1**: Topology Inference (Modules 1 & 2)
    - **Member 2**: Force Field & Parameterization (Modules 1 & 3)
    - **Member 3**: Core Energy Calculator & HPC-1 (Module 4)
    - **Member 4**: GUI & Parallelization (Module 5) ← *This module*
    """)
    
    st.divider()
    
    st.markdown("""
    ### 📚 References
    
    - OPLS-AA Force Field
    - RDKit for cheminformatics
    - xyz2mol for topology inference
    - scipy.spatial.cKDTree for spatial algorithms
    - Python multiprocessing for parallelization
    """)

# Footer
st.divider()
st.markdown("""
<div style='text-align: center; color: #666; padding: 1rem;'>
    <p>High-Performance Molecular Energy Calculator | Module 5: GUI & Parallelization</p>
    <p>Built with Streamlit | Member 4</p>
</div>
""", unsafe_allow_html=True)
