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
    initial_sidebar_state="collapsed"
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

# Sidebar - Minimal
# with st.sidebar:
#     st.header("🔧 Computation Mode")
#     computation_mode = st.radio(
#         "Select mode:",
#         ["Single Molecule", "Batch Processing (HPC-2)"],
#         help="Single: One molecule. Batch: Multiple molecules in parallel."
#     )

computation_mode = "Single Molecule"  # Default mode

# Main content area
tab1, tab2, tab3 = st.tabs(["💻 Calculate Energy", "🔧 YAML Builder", "📖 Documentation"])

# Tab 1: Energy Calculation
with tab1:
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("📁 Upload Molecular Geometry")
        uploaded_xyz = st.file_uploader(
            "Upload molecule file (.xyz, .pdb, or .mol)",
            type=['xyz', 'pdb', 'mol'],
            help="Supported formats: XYZ, PDB (Protein Data Bank), MOL (MDL Molfile)"
        )
        
        if uploaded_xyz:
            st.success(f"✅ Loaded: {uploaded_xyz.name}")
            
            # Determine MIME type based on file extension
            file_ext = uploaded_xyz.name.lower().split('.')[-1]
            mime_types = {
                'xyz': 'chemical/x-xyz',
                'pdb': 'chemical/x-pdb',
                'mol': 'chemical/x-mdl-molfile'
            }
            mime_type = mime_types.get(file_ext, 'text/plain')
            
            # Download button for molecule file
            st.download_button(
                label=f"⬇️ Download {file_ext.upper()} file",
                data=uploaded_xyz.getvalue(),
                file_name=uploaded_xyz.name,
                mime=mime_type
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
            # Parse molecule file content based on format
            file_content = uploaded_xyz.getvalue().decode('utf-8')
            file_ext = uploaded_xyz.name.lower().split('.')[-1]
            
            # Create 3D viewer with py3Dmol
            viewer = py3Dmol.view(width=800, height=500)
            
            # Add model based on file format
            viewer.addModel(file_content, file_ext)
            
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
            
            st.caption(f"🖱️ Use mouse to rotate, zoom, and pan the molecule | Format: {file_ext.upper()}")
                    
        except Exception as e:
            st.error(f"Error visualizing molecule: {str(e)}")
            file_ext = uploaded_xyz.name.lower().split('.')[-1]
            st.info(f"💡 Make sure your {file_ext.upper()} file is in the correct format.")
    
    st.divider()
    
    # Calculate button
    if uploaded_xyz is not None and uploaded_yaml is not None:
        
        # Core selection for parallel processing
        st.subheader("⚙️ Performance Settings")
        col_perf1, col_perf2 = st.columns(2)
        
        with col_perf1:
            import multiprocessing as mp
            max_cores = mp.cpu_count()
            n_cores = st.slider(
                "Number of CPU cores to use:",
                min_value=1,
                max_value=max_cores,
                value=max_cores,
                help=f"Your system has {max_cores} CPU cores available"
            )
        
        with col_perf2:
            num_copies = 1  # Default for single molecule mode
            if computation_mode == "Batch Processing (HPC-2)":
                num_copies = st.slider(
                    "Number of calculations (parallel):",
                    min_value=10,
                    max_value=1000,
                    value=100,
                    step=10,
                    help="Number of identical calculations for benchmarking"
                )
        
        col_btn1, col_btn2, col_btn3 = st.columns([1, 2, 1])
        with col_btn2:
            calculate_button = st.button(
                "🚀 Calculate Energy",
                type="primary",
                use_container_width=True
            )
        
        if calculate_button:
            # Save uploaded files temporarily with correct extension
            file_ext = uploaded_xyz.name.lower().split('.')[-1]
            temp_mol_file = f"temp.{file_ext}"
            
            with open(temp_mol_file, "wb") as f:
                f.write(uploaded_xyz.getvalue())
            with open("temp.yaml", "wb") as f:
                f.write(uploaded_yaml.getvalue())
            
            # STEP 1: Validate force field coverage
            with st.spinner("🔍 Validating force field parameters..."):
                validation = validate_force_field_coverage(temp_mol_file, "temp.yaml")
            
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
                    with st.spinner("🔄 Calculating energy and running performance benchmarks..."):
                        energy_breakdown = calculate_energy_with_breakdown(temp_mol_file, "temp.yaml", n_cores=n_cores)
                    
                    st.success("✅ Calculation complete!")
                    
                    # Display results
                    st.subheader("📈 Energy Results")
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
                    
                    # Performance Benchmarks
                    st.subheader("⚡ Performance Benchmarks")
                    timing = energy_breakdown['timing']
                    
                    st.info(f"**Molecule Size:** {timing['n_atoms']} atoms")
                    
                    col_bench1, col_bench2, col_bench3 = st.columns(3)
                    
                    with col_bench1:
                        st.metric(
                            label="Brute Force O(N²)",
                            value=f"{timing['brute_force_time']*1000:.2f} ms",
                            delta=None,
                            help="Time for naive all-pairs calculation"
                        )
                    
                    with col_bench2:
                        st.metric(
                            label="Single-Core Optimized",
                            value=f"{timing['single_core_time']*1000:.2f} ms",
                            delta=f"{timing['speedup_optimized']:.1f}x faster",
                            delta_color="normal",
                            help="Time with k-d tree optimization"
                        )
                    
                    with col_bench3:
                        st.metric(
                            label=f"Multi-Core ({timing['n_cores_used']} cores)",
                            value=f"{timing['multi_core_time']*1000:.2f} ms",
                            delta=f"{timing['speedup_parallel']:.1f}x faster",
                            delta_color="normal",
                            help="Time with parallel processing"
                        )
                    
                    # Speedup summary
                    with st.expander("📊 Performance Analysis"):
                        st.write("**Optimization Speedup:**")
                        st.write(f"- Brute force → Optimized: **{timing['speedup_optimized']:.2f}x** speedup")
                        st.write(f"- Single-core → Multi-core: **{timing['speedup_parallel']:.2f}x** speedup")
                        st.write(f"- Overall speedup: **{timing['speedup_optimized'] * timing['speedup_parallel']:.2f}x**")
                        st.write("")
                        st.write("**Complexity:**")
                        st.write(f"- Brute force: O(N²) = O({timing['n_atoms']}²) = {timing['n_atoms']**2:,} pair checks")
                        st.write(f"- Optimized: O(N log N) with spatial indexing (k-d tree)")
                        st.write(f"- Parallel: {timing['n_cores_used']} cores for distributed workload")
                    
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
                    
                    # Create list of molecule files (duplicates for benchmarking)
                    mol_files = [temp_mol_file] * num_copies
                    results = run_parallel_calculations(mol_files, "temp.yaml")
                    
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
    Upload a molecule file (.xyz, .pdb, or .mol) to get started, and we'll help you identify the required parameters.
    """)
    
    st.divider()
    
    # File upload for YAML builder
    col_builder1, col_builder2 = st.columns([2, 1])
    
    with col_builder1:
        st.markdown("### 📁 Step 1: Upload Molecule Structure")
        builder_xyz = st.file_uploader(
            "Upload molecule file (.xyz, .pdb, or .mol)",
            type=['xyz', 'pdb', 'mol'],
            help="We'll analyze this to help you define force field parameters",
            key="builder_xyz"
        )
    
    with col_builder2:
        st.markdown("### 💾 Quick Actions")
        if st.button("📋 Load Template", help="Start with default OPLS-AA template"):
            st.session_state['use_template'] = True
            st.success("✅ Template loaded!")
    
    if builder_xyz:
        # Determine file format and load appropriately
        file_ext = builder_xyz.name.lower().split('.')[-1]
        st.info(f"📄 Processing {file_ext.upper()} file: {builder_xyz.name}")
        
        # Parse molecule file
        file_content = builder_xyz.getvalue().decode('utf-8')
        
        # Save temporarily to load with RDKit
        temp_builder_file = f"temp_builder.{file_ext}"
        with open(temp_builder_file, "wb") as f:
            f.write(builder_xyz.getvalue())
        
        try:
            from calculator import load_molecule
            
            # Load molecule using universal loader
            temp_mol = load_molecule(temp_builder_file)
            num_atoms = len(temp_mol.atoms)
            
            atom_data = []
            for atom in temp_mol.atoms:
                atom_data.append({
                    'index': atom.index,
                    'element': atom.element,
                    'coords': atom.coords
                })

            
            # Display molecule info
            st.success(f"✅ Analyzed molecule: {num_atoms} atoms | Format: {file_ext.upper()}")
            
            # Count elements
            from collections import Counter
            import numpy as np
            element_counts = Counter([atom['element'] for atom in atom_data])
            
            # Try to use RDKit for better bond inference and SMARTS generation
            rdkit_mol = temp_mol.rdkit_mol
            auto_smarts = []
            
            if rdkit_mol is not None:
                try:
                    # Generate SMARTS patterns for each unique atom environment
                    def generate_atom_smarts(mol):
                        """Generate SMARTS patterns for each atom based on its environment"""
                        smarts_list = []
                        
                        for atom in mol.GetAtoms():
                            idx = atom.GetIdx()
                            symbol = atom.GetSymbol()
                            atomic_num = atom.GetAtomicNum()
                            
                            # Get atom properties
                            degree = atom.GetDegree()  # number of bonded neighbors
                            formal_charge = atom.GetFormalCharge()
                            
                            # Get neighbor information
                            neighbors = atom.GetNeighbors()
                            neighbor_symbols = sorted([n.GetSymbol() for n in neighbors])
                            
                            # For explicit H molecules (XYZ files), count explicit H neighbors
                            # GetTotalNumHs() returns 0 for explicit H, so we count manually
                            explicit_h_count = sum(1 for n in neighbors if n.GetSymbol() == 'H')
                            total_hs = explicit_h_count  # Use explicit H count for XYZ molecules
                            
                            # Build a descriptive SMARTS pattern
                            # For hydrogen, use #1 notation to avoid SMARTS ambiguity
                            # For other elements, use element symbol
                            if symbol == 'H':
                                smarts_parts = ['#1']  # Use atomic number for hydrogen
                            else:
                                smarts_parts = [symbol]
                            
                            # Add connectivity (X = total connections, D = degree)
                            if degree > 0:
                                smarts_parts.append(f"D{degree}")  # Use D instead of X for clarity
                            
                            # Add hydrogen count (only for non-hydrogen atoms)
                            if symbol != 'H' and total_hs > 0:
                                smarts_parts.append(f"H{total_hs}")
                            
                            # Join into SMARTS
                            base_smarts = f"[{';'.join(smarts_parts)}]"
                            
                            # Create extended SMARTS with one neighbor for better specificity
                            if neighbors:
                                # Choose most significant neighbor by priority:
                                # 1. Most electronegative heavy atom (O > N > Cl > S > C > ...)
                                # 2. This ensures CH2-O-H gets typed differently than CH2-C-H
                                
                                # Electronegativity order (simplified for common atoms)
                                ELECTRONEGATIVITY = {
                                    'F': 4.0, 'O': 3.5, 'N': 3.0, 'Cl': 3.0, 'Br': 2.8,
                                    'S': 2.5, 'C': 2.5, 'P': 2.1, 'H': 2.1, 'Si': 1.8
                                }
                                
                                heavy_neighbors = [n for n in neighbors if n.GetSymbol() != 'H']
                                
                                if heavy_neighbors:
                                    # Pick neighbor with highest electronegativity
                                    sig_neighbor = max(heavy_neighbors, 
                                                      key=lambda n: ELECTRONEGATIVITY.get(n.GetSymbol(), 2.0))
                                else:
                                    # Only H neighbors
                                    sig_neighbor = neighbors[0]
                                
                                n_symbol = sig_neighbor.GetSymbol()
                                n_atomic_num = sig_neighbor.GetAtomicNum()
                                n_degree = sig_neighbor.GetDegree()
                                
                                # Count explicit H neighbors for this neighbor atom too
                                n_neighbors = sig_neighbor.GetNeighbors()
                                n_explicit_h_count = sum(1 for nn in n_neighbors if nn.GetSymbol() == 'H')
                                n_hs = n_explicit_h_count
                                
                                # Use #1 for hydrogen neighbors too
                                if n_symbol == 'H':
                                    neighbor_smarts_parts = ['#1']
                                else:
                                    neighbor_smarts_parts = [n_symbol]
                                    
                                if n_degree > 0:
                                    neighbor_smarts_parts.append(f"D{n_degree}")  # Use D
                                if n_symbol != 'H' and n_hs > 0:
                                    neighbor_smarts_parts.append(f"H{n_hs}")
                                
                                neighbor_smarts = f"[{';'.join(neighbor_smarts_parts)}]"
                                extended_smarts = f"{base_smarts}{neighbor_smarts}"
                            else:
                                extended_smarts = base_smarts
                            
                            # Create description
                            neighbor_str = ", ".join(neighbor_symbols) if neighbor_symbols else "none"
                            description = f"{symbol} (degree={degree}, H={total_hs}, neighbors: {neighbor_str})"
                            
                            # Create unique key based on element + sorted neighbor elements
                            # This groups atoms by their complete chemical environment
                            neighbor_key = "_".join(sorted(neighbor_symbols)) if neighbor_symbols else "isolated"
                            environment_key = f"{symbol}_{neighbor_key}"
                            
                            smarts_list.append({
                                'atom_idx': idx,
                                'element': symbol,
                                'smarts': extended_smarts,
                                'base_smarts': base_smarts,
                                'description': description,
                                'degree': degree,
                                'total_hs': total_hs,
                                'neighbors': neighbor_symbols,
                                'environment_key': environment_key  # NEW: unique key for grouping
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
                    st.error(f"❌ RDKit SMARTS generation failed: {str(e)}")
                    st.exception(e)  # Show full traceback
                    rdkit_mol = None
                    auto_smarts = []  # Ensure it's empty on failure
            
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
                    # Group by BASE SMARTS pattern
                    # Base SMARTS like [C;D4;H2] or [#1;D1] encode the local chemistry
                    # This groups all atoms with identical local environment
                    key = item['base_smarts']
                    
                    if key not in smarts_groups:
                        smarts_groups[key] = make_default_dict()
                        smarts_groups[key]['smarts'] = item['base_smarts']
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
                            # Generate intuitive default type names based on element
                            # Users should edit these to match their specific force field
                            if rdkit_success and i < len(unique_atom_types_rdkit):
                                elem = unique_atom_types_rdkit[i]['element']
                                # Use element-based naming for clarity
                                default_type = f"{elem}_type_{i+1}"
                            else:
                                default_type = f"atom_type_{i+1}"
                            
                            type_name = st.text_input(
                                "Type Name",
                                value=default_type,
                                key=f"type_name_{i}",
                                help="Edit to match your force field (e.g., opls_154, opls_155, etc.)"
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
                atom_idx_to_smarts = {}
                
                # APPROACH 1: Use RDKit-detected groupings if available
                if rdkit_success and unique_atom_types_rdkit and len(atom_types_list) == len(unique_atom_types_rdkit):
                    # Map each atom to its assigned type name using RDKit groupings
                    for i, auto_type in enumerate(unique_atom_types_rdkit):
                        type_name = atom_types_list[i]['type_name']
                        for atom_idx in auto_type['indices']:
                            atom_idx_to_type[atom_idx] = type_name
                            atom_idx_to_element[atom_idx] = auto_type['element']
                            atom_idx_to_smarts[atom_idx] = auto_type['smarts']
                else:
                    # APPROACH 2: Map atoms sequentially to user-defined types
                    # This handles when RDKit fails OR user changed number of types
                    for idx, atom in enumerate(atom_data):
                        if idx < len(atom_types_list):
                            # Use the type name from atom_types_list
                            atom_idx_to_type[idx] = atom_types_list[idx]['type_name']
                            atom_idx_to_element[idx] = atom['element']
                            atom_idx_to_smarts[idx] = atom_types_list[idx].get('smarts', f"[{atom['element']}]")
                        else:
                            # Fallback: more types needed than defined
                            atom_idx_to_type[idx] = f"atom_type_{idx+1}"
                            atom_idx_to_element[idx] = atom['element']
                            atom_idx_to_smarts[idx] = f"[{atom['element']}]"
            
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
                num_bond_types = st.number_input("Number of bond types:", min_value=0, value=default_num_bonds, key="num_bond_types")
                
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
                        # AUTO-GENERATE bond name from current atom type names - don't allow editing
                        # This ensures consistency with atom_types section
                        bond_name = default_bond_name
                        st.text(f"Bond Type Name: {bond_name}")
                        st.caption(f"💡 Auto-generated from atom types. Element representation: {elem_hint}")
                        st.caption("⚠️ To change this, modify the atom type names in the 'Define Atom Types' tab.")
                        
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
                num_angle_types = st.number_input("Number of angle types:", min_value=0, value=default_num_angles, key="num_angle_types")
                
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
                        # AUTO-GENERATE angle name from current atom type names - don't allow editing
                        angle_name = default_angle_name
                        st.text(f"Angle Type Name: {angle_name}")
                        st.caption(f"💡 Auto-generated from atom types. Element representation: {elem_hint}")
                        st.caption("⚠️ To change this, modify the atom type names in the 'Define Atom Types' tab.")
                        
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
                num_dihedral_types = st.number_input("Number of dihedral types:", min_value=0, value=default_num_dihedrals, key="num_dihedral_types")
                
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
                        # AUTO-GENERATE dihedral name from current atom type names - don't allow editing
                        dihedral_name = default_dihedral_name
                        st.text(f"Dihedral Type Name: {dihedral_name}")
                        st.caption(f"💡 Auto-generated from atom types. Element representation: {elem_hint}")
                        st.caption("⚠️ To change this, modify the atom type names in the 'Define Atom Types' tab.")
                        
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

# Tab 3: Documentation
with tab3:
    st.markdown("# 📖 Comprehensive Documentation")
    st.markdown("### High-Performance Molecular Energy Calculator")
    
    # Table of Contents
    st.markdown("""
    ## Table of Contents
    1. [Input Formats](#input-formats)
    2. [Force Field Format](#force-field-format)
    3. [YAML Builder](#yaml-builder)
    4. [Molecule Visualizer](#molecule-visualizer)
    5. [YAML Verification](#yaml-verification)
    6. [Energy Calculation](#energy-calculation)
    7. [Performance Optimizations](#performance-optimizations)
    """)
    
    st.markdown("---")
    
    # Section 1: Input Formats
    st.markdown("## 1. Input Formats for Molecules")
    st.markdown("""
    The calculator supports three molecular structure file formats:
    
    ### 1.1 XYZ Format (.xyz)
    
    The XYZ format is a simple text-based format containing:
    - **Line 1**: Number of atoms (integer)
    - **Line 2**: Comment line (usually molecule name or description)
    - **Lines 3+**: Atom symbol followed by x, y, z coordinates in Ångströms
    
    **Example (Ethanol):**
    ```
    9
    Ethanol molecule
    C    0.000000    0.000000    0.000000
    C    1.522900    0.000000    0.000000
    O    2.018100   -1.324300    0.000000
    H   -0.386600    1.027100    0.000000
    H   -0.386600   -0.513500   -0.889200
    H   -0.386600   -0.513500    0.889200
    H    1.909500    0.513500    0.889200
    H    1.909500    0.513500   -0.889200
    H    1.631900   -1.769100   -0.000000
    ```
    
    **Coordinate Conversion**: XYZ files use Ångströms (Å), which are automatically converted to nanometers (nm) internally:
    $$1 \\text{ Å} = 0.1 \\text{ nm}$$
    
    ### 1.2 PDB Format (.pdb)
    
    Protein Data Bank format with explicit connectivity information. PDB files contain:
    - ATOM/HETATM records with 3D coordinates
    - CONECT records defining bonds
    - Residue and chain information
    
    **Advantages**: 
    - Explicit bond information (no inference needed)
    - Standard format for biomolecules
    - Contains secondary structure information
    
    ### 1.3 MOL Format (.mol)
    
    MDL Molfile format (also called SDF) containing:
    - Atom block with coordinates and types
    - Bond block with explicit connectivity
    - Properties block
    
    **Advantages**:
    - Explicit bond orders (single, double, triple, aromatic)
    - Stereochemistry information
    - Standard in computational chemistry
    
    ### 1.4 Bond Inference for XYZ Files
    
    Since XYZ files lack bond information, bonds are inferred using **covalent radii**:
    
    $$d_{ij} \\leq f \\times (r_i + r_j)$$
    
    Where:
    - $d_{ij}$ = distance between atoms $i$ and $j$
    - $r_i, r_j$ = covalent radii of atoms
    - $f$ = scaling factor (typically 1.2)
    
    **Covalent Radii** (in Ångströms):
    - H: 0.31, C: 0.76, N: 0.71, O: 0.66, S: 1.05, P: 1.07, Cl: 1.02
    """)
    
    st.markdown("---")
    
    # Section 2: Force Field Format
    st.markdown("## 2. Force Field Format (YAML)")
    st.markdown("""
    The force field parameters are stored in YAML format with four main sections:
    
    ### 2.1 YAML Structure
    
    ```yaml
    atom_types:
      - smarts: "[C;D4;H3]"
        type_name: "C_type_1"
        charge: -0.18
        sigma: 0.350
        epsilon: 0.276
    
    bond_types:
      C_type_1-C_type_2: [224262.4, 0.1529]  # [kb, b0]
    
    angle_types:
      H_type_4-C_type_1-C_type_2: [313.8, 1.91114]  # [k_theta, theta0]
    
    dihedral_types:
      H_type_4-C_type_1-C_type_2-O_type_3: [0.0, 0.0, 1.2552, 0.0]  # [V1, V2, V3, V4]
    ```
    
    ### 2.2 Atom Types
    
    Each atom type is defined using **SMARTS patterns** for automatic atom typing:
    
    - **smarts**: SMARTS pattern matching chemical environment
    - **type_name**: Unique identifier for this atom type
    - **charge**: Partial charge in elementary charge units (e)
    - **sigma** (σ): Lennard-Jones collision diameter (nm)
    - **epsilon** (ε): Lennard-Jones well depth (kJ/mol)
    
    **SMARTS Pattern Examples**:
    - `[C;D4;H3]`: Carbon with 4 bonds, 3 hydrogens (methyl -CH₃)
    - `[C;D4;H2]`: Carbon with 4 bonds, 2 hydrogens (methylene -CH₂-)
    - `[O;D2;H1]`: Oxygen with 2 bonds, 1 hydrogen (hydroxyl -OH)
    - `[#1;D1]`: Hydrogen with 1 bond
    
    ### 2.3 Bond Types
    
    Bond parameters use harmonic potential:
    
    $$V_{\\text{bond}}(b) = \\frac{1}{2} k_b (b - b_0)^2$$
    
    Format: `type1-type2: [kb, b0]`
    - **kb**: Force constant (kJ/mol/nm²)
    - **b0**: Equilibrium bond length (nm)
    
    ### 2.4 Angle Types
    
    Angle parameters use harmonic potential:
    
    $$V_{\\text{angle}}(\\theta) = \\frac{1}{2} k_\\theta (\\theta - \\theta_0)^2$$
    
    Format: `type1-type2-type3: [k_theta, theta0]`
    - **k_theta**: Force constant (kJ/mol/rad²)
    - **theta0**: Equilibrium angle (radians)
    
    ### 2.5 Dihedral Types
    
    Dihedral parameters use OPLS Fourier series:
    
    $$V_{\\text{dihedral}}(\\phi) = \\frac{V_1}{2}(1 + \\cos\\phi) + \\frac{V_2}{2}(1 - \\cos 2\\phi) + \\frac{V_3}{2}(1 + \\cos 3\\phi) + \\frac{V_4}{2}(1 - \\cos 4\\phi)$$
    
    Format: `type1-type2-type3-type4: [V1, V2, V3, V4]`
    - All parameters in kJ/mol
    - $\\phi$ is the dihedral angle
    
    ### 2.6 Units Summary
    
    | Parameter | Unit |
    |-----------|------|
    | Energy | kJ/mol |
    | Length | nm (nanometers) |
    | Angle | radians |
    | Charge | e (elementary charge) |
    | Force constant (bond) | kJ/mol/nm² |
    | Force constant (angle) | kJ/mol/rad² |
    """)
    
    st.markdown("---")
    
    # Section 3: YAML Builder
    st.markdown("## 3. YAML Builder - Automatic Force Field Generation")
    st.markdown("""
    The YAML Builder automates the creation of force field parameter files.
    
    ### 3.1 Workflow
    
    1. **Upload Molecule**: Upload XYZ/PDB/MOL file
    2. **RDKit Analysis**: Automatic detection of molecular topology
    3. **SMARTS Generation**: Chemical environment patterns for each atom
    4. **Atom Type Grouping**: Group atoms by identical environments
    5. **Topology Detection**: Identify all bonds, angles, dihedrals
    6. **Parameter Input**: User provides force field parameters
    7. **YAML Export**: Download complete parameter file
    
    ### 3.2 RDKit Integration
    
    The builder uses RDKit to:
    - **Detect bonds** from 3D coordinates (for XYZ files)
    - **Analyze chemical environments** around each atom
    - **Generate SMARTS patterns** encoding local chemistry
    
    **SMARTS Pattern Construction**:
    
    For each atom, the pattern includes:
    - Element symbol (or #N for atomic number)
    - Degree (D): Number of bonded neighbors
    - Hydrogen count (H): Explicit hydrogen neighbors
    - Connectivity to most electronegative neighbor
    
    Example: `[C;D4;H2][O;D2;H1]` = Carbon bonded to 4 atoms (2 H) connected to oxygen
    
    ### 3.3 Atom Type Grouping
    
    Atoms are grouped by their **base SMARTS pattern**:
    
    - All atoms matching `[C;D4;H3]` → Same type (methyl carbons)
    - All atoms matching `[#1;D1]` → Same type (hydrogens)
    
    This reduces redundancy: a molecule with 100 H atoms needs only 1 H atom type definition.
    
    ### 3.4 Type Name Consistency
    
    **Critical Design**: Bond/angle/dihedral type names are **auto-generated** from atom type names.
    
    If atom types are:
    - `C_type_1`, `C_type_2`, `O_type_3`, `H_type_4`
    
    Then bonds are automatically named:
    - `C_type_1-C_type_2` (C-C bond)
    - `C_type_2-O_type_3` (C-O bond)
    - `C_type_1-H_type_4` (C-H bond)
    
    This ensures **100% consistency** between atom types and topology parameters.
    
    ### 3.5 Default Parameters
    
    The builder provides typical OPLS-AA default values:
    - **Bonds**: kb ≈ 224000-284000 kJ/mol/nm², b0 ≈ 0.109-0.153 nm
    - **Angles**: k_θ ≈ 313-418 kJ/mol/rad², θ0 ≈ 109.5° (1.911 rad)
    - **LJ Parameters**: σ ≈ 0.25-0.37 nm, ε ≈ 0.066-0.711 kJ/mol
    
    Users should replace these with actual force field values for production calculations.
    """)
    
    st.markdown("---")
    
    # Section 4: Molecule Visualizer  
    st.markdown("## 4. Molecule Visualizer (3Dmol.js)")
    st.markdown("""
    Interactive 3D visualization powered by **3Dmol.js**.
    
    ### 4.1 Features
    
    - **Rendering Styles**:
      - Stick: Shows bonds as sticks
      - Sphere: Van der Waals spheres
      - Ball and Stick: Atoms as spheres, bonds as sticks
      - Line: Simple wireframe
    
    - **Atom Coloring**:
      - CPK: Standard element colors (C=gray, O=red, N=blue, H=white)
      - Custom schemes supported
    
    - **Interactive Controls**:
      - Rotate: Click and drag
      - Zoom: Mouse wheel
      - Pan: Right-click drag
    
    ### 4.2 Implementation
    
    The visualizer converts molecule coordinates to XYZ format and renders using py3Dmol:
    
    ```python
    view = py3Dmol.view(width=800, height=600)
    view.addModel(xyz_string, "xyz")
    view.setStyle({'stick': {}, 'sphere': {'scale': 0.3}})
    view.setBackgroundColor('white')
    view.zoomTo()
    ```
    
    ### 4.3 Coordinate System
    
    - Input coordinates in Ångströms (standard XYZ)
    - Internally converted to nm for calculations
    - Visualizer displays in Ångströms (standard molecular graphics)
    """)
    
    st.markdown("---")
    
    # Section 5: YAML Verification
    st.markdown("## 5. YAML Verification - Force Field Coverage Analysis")
    st.markdown("""
    Before energy calculation, the system validates that the force field covers all molecular features.
    
    ### 5.1 Coverage Metrics
    
    The validator checks coverage for:
    
    1. **Atoms**: % of atoms successfully typed by SMARTS patterns
    2. **Bonds**: % of bonds with defined parameters
    3. **Angles**: % of angles with defined parameters  
    4. **Dihedrals**: % of dihedrals with defined parameters
    
    ### 5.2 SMARTS Matching
    
    For each atom in the molecule:
    
    1. Try to match each SMARTS pattern in atom_types
    2. If match found → assign type_name to atom
    3. If no match → atom remains untyped (generic fallback)
    
    **RDKit Matching**:
    ```python
    pattern = Chem.MolFromSmarts(smarts_string)
    matches = mol.GetSubstructMatches(pattern)
    ```
    
    ### 5.3 Parameter Lookup
    
    For each bond (i,j):
    - Get atom types: `type_i`, `type_j`
    - Create sorted key: `f"{min(type_i, type_j)}-{max(type_i, type_j)}"`
    - Look up in `bond_types` dictionary
    
    Same process for angles (3 atoms) and dihedrals (4 atoms).
    
    ### 5.4 Missing Parameter Handling
    
    If parameters are missing:
    - **Warning displayed** with list of missing types
    - **Generic fallback** used (basic element-based estimates)
    - **Coverage % reported** to user
    
    **Example Output**:
    ```
    ✅ Atom Coverage: 100.0% (9/9)
    ✅ Bond Coverage: 100.0% (8/8)
    ✅ Angle Coverage: 100.0% (13/13)
    ✅ Dihedral Coverage: 100.0% (12/12)
    ```
    
    ### 5.5 Validation Report
    
    The validator generates a detailed report showing:
    - Total count of each topology type
    - Number successfully parameterized
    - Coverage percentage
    - List of missing parameters (if any)
    - Warnings for low coverage (<90%)
    """)
    
    st.markdown("---")
    
    # Section 6: Energy Calculation - Detailed Theory
    st.markdown("## 6. Molecular Energy Calculation")
    st.markdown(r"""
    The total potential energy is calculated using classical molecular mechanics:
    
    $$E_{\text{total}} = E_{\text{bonded}} + E_{\text{non-bonded}}$$
    
    ### 6.1 Bonded Interactions
    
    #### 6.1.1 Bond Stretching Energy
    
    Harmonic potential modeling covalent bond deformation:
    
    $$E_{\text{bonds}} = \sum_{\text{bonds}} \frac{1}{2} k_b (b - b_0)^2$$
    
    Where:
    - $b$ = current bond length (nm)
    - $b_0$ = equilibrium bond length (nm)
    - $k_b$ = force constant (kJ/mol/nm²)
    
    **Calculation Steps**:
    1. For each bond (i, j):
       - Calculate distance: $b = \sqrt{(x_i - x_j)^2 + (y_i - y_j)^2 + (z_i - z_j)^2}$
       - Get parameters from force field: $(k_b, b_0)$
       - Compute energy: $E = \frac{1}{2} k_b (b - b_0)^2$
    2. Sum over all bonds
    
    **Example** (C-C bond):
    - $k_b = 224262$ kJ/mol/nm²
    - $b_0 = 0.1529$ nm
    - If $b = 0.1535$ nm: $E = \frac{1}{2} \times 224262 \times (0.1535 - 0.1529)^2 = 0.404$ kJ/mol
    
    #### 6.1.2 Angle Bending Energy
    
    Harmonic potential for bond angle deformation:
    
    $$E_{\text{angles}} = \sum_{\text{angles}} \frac{1}{2} k_\theta (\theta - \theta_0)^2$$
    
    Where:
    - $\theta$ = current angle (radians)
    - $\theta_0$ = equilibrium angle (radians)
    - $k_\theta$ = force constant (kJ/mol/rad²)
    
    **Calculation Steps**:
    1. For each angle (i, j, k) where j is the central atom:
       - Calculate vectors: $\vec{v}_1 = \vec{r}_i - \vec{r}_j$, $\vec{v}_2 = \vec{r}_k - \vec{r}_j$
       - Compute angle: $\theta = \arccos\left(\frac{\vec{v}_1 \cdot \vec{v}_2}{|\vec{v}_1| |\vec{v}_2|}\right)$
       - Get parameters: $(k_\theta, \theta_0)$
       - Compute energy: $E = \frac{1}{2} k_\theta (\theta - \theta_0)^2$
    2. Sum over all angles
    
    **Example** (H-C-C angle):
    - $k_\theta = 313.8$ kJ/mol/rad²
    - $\theta_0 = 1.911$ rad (109.5°)
    - If $\theta = 1.92$ rad: $E = \frac{1}{2} \times 313.8 \times (1.92 - 1.911)^2 = 0.127$ kJ/mol
    
    #### 6.1.3 Dihedral Torsion Energy
    
    OPLS Fourier series for rotation around bonds:
    
    $$E_{\text{dihedrals}} = \sum_{\text{dihedrals}} \left[\frac{V_1}{2}(1 + \cos\phi) + \frac{V_2}{2}(1 - \cos 2\phi) + \frac{V_3}{2}(1 + \cos 3\phi) + \frac{V_4}{2}(1 - \cos 4\phi)\right]$$
    
    Where:
    - $\phi$ = dihedral angle (radians)
    - $V_1, V_2, V_3, V_4$ = Fourier coefficients (kJ/mol)
    
    **Calculation Steps**:
    1. For each dihedral (i, j, k, l):
       - Calculate normal vectors to planes (i,j,k) and (j,k,l)
       - Compute dihedral angle $\phi$ using normal vectors
       - Get Fourier coefficients: $(V_1, V_2, V_3, V_4)$
       - Compute energy using OPLS formula
    2. Sum over all dihedrals
    
    **Dihedral Angle Calculation**:
    ```
    b1 = r_j - r_i
    b2 = r_k - r_j  
    b3 = r_l - r_k
    n1 = b1 × b2  (normal to plane 1)
    n2 = b2 × b3  (normal to plane 2)
    φ = atan2((n1 × n2)·b2/|b2|, n1·n2)
    ```
    
    **Example** (H-C-C-O dihedral):
    - $V_1 = 0.0$, $V_2 = 0.0$, $V_3 = 1.2552$, $V_4 = 0.0$ kJ/mol
    - If $\phi = 60°$ (1.047 rad): $E = \frac{1.2552}{2}(1 + \cos(3 \times 1.047)) = 0.628$ kJ/mol
    
    ### 6.2 Non-Bonded Interactions
    
    #### 6.2.1 Van der Waals (Lennard-Jones)
    
    12-6 Lennard-Jones potential between non-bonded atom pairs:
    
    $$E_{\text{LJ}} = \sum_{i<j} 4\epsilon_{ij} \left[\left(\frac{\sigma_{ij}}{r_{ij}}\right)^{12} - \left(\frac{\sigma_{ij}}{r_{ij}}\right)^6\right]$$
    
    Where:
    - $r_{ij}$ = distance between atoms i and j (nm)
    - $\epsilon_{ij}$ = well depth (kJ/mol)
    - $\sigma_{ij}$ = collision diameter (nm)
    
    **Combining Rules** (Lorentz-Berthelot):
    
    $$\sigma_{ij} = \frac{\sigma_i + \sigma_j}{2}$$
    
    $$\epsilon_{ij} = \sqrt{\epsilon_i \times \epsilon_j}$$
    
    **Calculation Steps**:
    1. Build list of non-bonded pairs (exclude 1-2, 1-3, 1-4 neighbors)
    2. For each pair (i, j):
       - Calculate distance: $r_{ij}$
       - Apply combining rules to get $\sigma_{ij}$, $\epsilon_{ij}$
       - Compute LJ energy
    3. Sum over all pairs
    
    **Example** (C···H interaction at 0.3 nm):
    - $\sigma_C = 0.355$ nm, $\sigma_H = 0.242$ nm → $\sigma_{CH} = 0.2985$ nm
    - $\epsilon_C = 0.276$ kJ/mol, $\epsilon_H = 0.126$ kJ/mol → $\epsilon_{CH} = 0.187$ kJ/mol
    - $E = 4 \times 0.187 \times [(0.2985/0.3)^{12} - (0.2985/0.3)^6] = -0.183$ kJ/mol
    
    #### 6.2.2 Electrostatic (Coulomb)
    
    Coulomb potential between partial charges:
    
    $$E_{\text{elec}} = \sum_{i<j} \frac{q_i q_j}{4\pi\epsilon_0 r_{ij}} = \sum_{i<j} 138.935 \frac{q_i q_j}{r_{ij}}$$
    
    Where:
    - $q_i, q_j$ = partial charges (elementary charge e)
    - $r_{ij}$ = distance (nm)
    - Constant: $\frac{1}{4\pi\epsilon_0} = 138.935$ kJ·nm/(mol·e²)
    
    **Calculation Steps**:
    1. For each non-bonded pair (i, j):
       - Get partial charges from force field
       - Calculate distance
       - Compute electrostatic energy
    2. Sum over all pairs
    
    **Example** (C⁻···H⁺ at 0.25 nm):
    - $q_C = -0.18$ e, $q_H = +0.06$ e
    - $E = 138.935 \times \frac{(-0.18) \times 0.06}{0.25} = -6.00$ kJ/mol
    
    ### 6.3 Exclusions
    
    **1-2 exclusions**: Bonded atoms (bond between i-j)
    **1-3 exclusions**: Atoms separated by 2 bonds (i-j-k)
    **1-4 exclusions**: Atoms separated by 3 bonds (i-j-k-l)
    
    These pairs are **excluded** from non-bonded calculations to avoid double-counting with bonded terms.
    """)
    
    st.markdown("---")
    
    # Section 7: Performance Optimizations
    st.markdown("## 7. Performance Optimizations")
    
    st.markdown(r"""
    ### 7.1 Naive Approach - O(N²) Complexity
    
    **Algorithm**: Nested loop over all atom pairs
    
    ```python
    def calculate_nonbonded_naive(atoms):
        energy = 0.0
        for i in range(len(atoms)):
            for j in range(i+1, len(atoms)):
                if (i,j) not in exclusions:
                    r = distance(atoms[i], atoms[j])
                    energy += lennard_jones(r, sigma, epsilon)
                    energy += coulomb(r, q_i, q_j)
        return energy
    ```
    
    **Complexity**: $O(N^2)$ where N = number of atoms
    
    **Performance**:
    - 10 atoms: ~45 pairs → Fast
    - 100 atoms: ~4,950 pairs → Slow
    - 1,000 atoms: ~499,500 pairs → Very slow
    - 10,000 atoms: ~49,995,000 pairs → Prohibitive
    
    **Problem**: Computational cost scales quadratically. Large molecules become intractable.
    
    ### 7.2 Single-Core Optimization - Spatial Trees O(N log N)
    
    **Key Insight**: Van der Waals interactions decay rapidly with distance ($\sim r^{-6}$). 
    Atoms beyond a cutoff distance contribute negligibly to energy.
    
    **Algorithm**: Use **spatial tree** (k-d tree) for efficient neighbor search
    
    ```python
    from scipy.spatial import cKDTree
    
    def calculate_nonbonded_optimized(atoms, cutoff=1.0):
        # Build spatial tree
        coords = np.array([atom.coords for atom in atoms])
        tree = cKDTree(coords)
        
        # Query pairs within cutoff
        pairs = tree.query_pairs(r=cutoff)
        
        energy = 0.0
        for (i, j) in pairs:
            if (i,j) not in exclusions:
                r = distance(atoms[i], atoms[j])
                energy += lennard_jones(r, sigma, epsilon)
                energy += coulomb(r, q_i, q_j)
        return energy
    ```
    
    **Complexity**: $O(N \log N)$ for tree construction + $O(M)$ for pair iteration
    
    **Cutoff Distance**: Typically 1.0-1.2 nm reduces pairs from $O(N^2)$ to $O(N)$ on average
    
    **Speedup**: 10x for 100 atoms, 50x for 1,000 atoms, 285x for 10,000 atoms!
    
    ### 7.3 Multi-Core Optimization - Parallel Processing
    
    **Key Insight**: Different molecules are independent → embarrassingly parallel
    
    **Algorithm**: Distribute molecules across CPU cores using `multiprocessing`
    
    ```python
    from multiprocessing import Pool, cpu_count
    
    def calculate_batch_parallel(molecules, force_field):
        n_cores = cpu_count()
        with Pool(processes=n_cores) as pool:
            results = pool.map(
                partial(calculate_energy, ff=force_field),
                molecules
            )
        return results
    ```
    
    **Expected Speedup** (Amdahl's Law):
    
    $$S(P) \approx P \times \eta$$
    
    Where $P$ = cores, $\eta$ = efficiency (0.85-0.95)
    
    **Speedup Table**:
    
    | Cores | Ideal | Actual | Efficiency |
    |-------|-------|--------|------------|
    | 1 | 1.0x | 1.0x | 100% |
    | 2 | 2.0x | 1.9x | 95% |
    | 4 | 4.0x | 3.7x | 92% |
    | 8 | 8.0x | 7.2x | 90% |
    | 16 | 16.0x | 13.5x | 84% |
    
    ### 7.4 Combined Performance
    
    **Maximum Performance**: Spatial tree + multi-core
    
    $$\text{Speedup}_{\text{total}} = S_{\text{tree}} \times S_{\text{parallel}}$$
    
    **Example** (4 cores, 1000-atom molecules):
    - Tree speedup: 50x
    - Parallel speedup: 3.7x  
    - **Total: 185x faster than naive!**
    """)
    
    st.markdown("---")
    
    st.markdown("""
    ## 8. Summary
    
    This molecular energy calculator implements:
    
    ✅ **Multi-format input**: XYZ, PDB, MOL files  
    ✅ **SMARTS-based typing**: Automatic atom type assignment  
    ✅ **Complete force field**: Bonds, angles, dihedrals, LJ, Coulomb  
    ✅ **YAML builder**: Automated parameter file generation  
    ✅ **Validation**: Force field coverage verification  
    ✅ **Visualization**: Interactive 3D molecule viewer  
    ✅ **HPC-1 optimization**: Spatial trees for O(N log N) scaling  
    ✅ **HPC-2 optimization**: Multi-core parallelization  
    ✅ **User-friendly GUI**: Streamlit-based interface  
    
    ### Applications
    
    - Drug design and screening
    - Molecular dynamics pre-processing
    - Force field development and testing
    - Educational demonstrations
    - High-throughput virtual screening
    """)

# Footer
st.divider()
st.markdown("""
<div style='text-align: center; color: #666; padding: 1rem;'>
    <p>High-Performance Molecular Energy Calculator</p>
    <p>Made by Arunangshu, Ankana, Rohit and Prashanthi.</p>
</div>
""", unsafe_allow_html=True)
