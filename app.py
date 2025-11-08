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
tab1, tab2, tab3 = st.tabs(["💻 Calculate Energy", "📊 Benchmark", "📖 Documentation"])

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
            
            try:
                if computation_mode == "Single Molecule":
                    # Single calculation
                    with st.spinner("🔄 Calculating energy..."):
                        # TODO: Replace with actual calculation function
                        # from calculator import calculate_single_molecule_energy
                        # result = calculate_single_molecule_energy("temp.xyz", "temp.yaml")
                        
                        # Placeholder simulation
                        time.sleep(1.5)
                        result = ("temp.xyz", -123.4567)  # Dummy result
                    
                    st.success("✅ Calculation complete!")
                    
                    # Display results
                    st.subheader("📈 Results")
                    col_res1, col_res2, col_res3 = st.columns(3)
                    
                    with col_res1:
                        st.metric(
                            label="Total Energy",
                            value=f"{result[1]:.4f} kJ/mol",
                            delta=None
                        )
                    
                    with col_res2:
                        st.metric(
                            label="Bonded Energy",
                            value="-45.23 kJ/mol",
                            delta=None
                        )
                    
                    with col_res3:
                        st.metric(
                            label="Non-bonded Energy",
                            value="-78.23 kJ/mol",
                            delta=None
                        )
                    
                    # Detailed breakdown
                    with st.expander("🔍 Detailed Energy Breakdown"):
                        breakdown_data = {
                            "Component": ["Bonds", "Angles", "Dihedrals", "Van der Waals", "Electrostatic", "Total"],
                            "Energy (kJ/mol)": [-12.45, -18.32, -14.46, -35.67, -42.56, -123.46]
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
                    
                    # Simulate parallel computation
                    start_time = time.time()
                    
                    # TODO: Replace with actual parallel function
                    # from calculator import run_parallel_calculations
                    # results = run_parallel_calculations(["temp.xyz"] * num_copies, "temp.yaml")
                    
                    # Placeholder simulation
                    for i in range(100):
                        time.sleep(0.02)
                        progress_bar.progress(i + 1)
                    
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
                        sample_data = {
                            "Molecule": [f"Molecule {i+1}" for i in range(10)],
                            "Total Energy (kJ/mol)": [-123.46 + (i * 0.1) for i in range(10)]
                        }
                        st.dataframe(pd.DataFrame(sample_data), use_container_width=True)
            
            finally:
                # Cleanup temp files
                if os.path.exists("temp.xyz"):
                    os.remove("temp.xyz")
                if os.path.exists("temp.yaml"):
                    os.remove("temp.yaml")
    
    else:
        st.info("👆 Please upload both a .xyz geometry file and a .yaml force field file to begin.")

# Tab 2: Benchmark
with tab2:
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

# Tab 3: Documentation
with tab3:
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
