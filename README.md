# ⚛️ High-Performance Molecular Energy Calculator

A web-based application for calculating molecular potential energy using classical force fields (OPLS-AA). This project implements a complete 5-module pipeline with HPC optimizations and parallel processing capabilities.

## 🎯 Project Overview

This calculator computes the total potential energy of molecular systems using:
- **Bonded interactions**: Bonds, Angles, Dihedrals
- **Non-bonded interactions**: Van der Waals (Lennard-Jones), Electrostatic (Coulomb)
- **HPC-1 Optimization**: Spatial tree algorithm (O(N) neighbor search)
- **HPC-2 Parallelization**: Multi-core processing for batch calculations

## 🏗️ Architecture

### 5-Module Pipeline

1. **Module 1**: Input/Output & Data Structures
   - `.xyz` molecular geometry parser
   - `.yaml` force field parameter parser
   
2. **Module 2**: Topology Inference Engine
   - Bond detection from coordinates
   - Angle and dihedral identification
   
3. **Module 3**: Parameter Assignment Engine
   - Atom typing using SMARTS patterns
   - Force field parameter lookup
   
4. **Module 4**: Core Energy Calculator (HPC-1)
   - Bonded energy calculations
   - Optimized non-bonded calculations
   
5. **Module 5**: Parallelization & GUI (HPC-2) ⭐ **This Module**
   - Streamlit web interface
   - Multi-core parallel processing

## 👥 Team Responsibilities

- **Member 1**: Topology Inference (Modules 1 & 2)
- **Member 2**: Force Field & Parameterization (Modules 1 & 3)
- **Member 3**: Core Energy Calculator & HPC-1 (Module 4)
- **Member 4**: GUI & Parallelization (Module 5) ← **Current Focus**

## 🚀 Quick Start

### Local Installation

1. **Clone the repository**
   ```bash
   git clone https://github.com/arunangshu/hpmec.git
   cd hpmec
   ```

2. **Install dependencies**
   ```bash
   pip install -r requirements.txt
   ```

3. **Run the Streamlit app**
   ```bash
   streamlit run app.py
   ```

4. **Open your browser**
   - The app will automatically open at `http://localhost:8501`

## ☁️ Deploy to Streamlit Cloud

### Step 1: Push to GitHub
Your code is already on GitHub at: https://github.com/arunangshu/hpmec

### Step 2: Deploy on Streamlit Cloud

1. **Go to Streamlit Cloud**
   - Visit: https://share.streamlit.io/

2. **Sign in with GitHub**
   - Click "Sign in with GitHub"
   - Authorize Streamlit Cloud

3. **Deploy your app**
   - Click "New app"
   - Repository: `arunangshu/hpmec`
   - Branch: `main`
   - Main file path: `app.py`
   - Click "Deploy!"

4. **Wait for deployment**
   - Streamlit will install dependencies from `requirements.txt`
   - Your app will be live at: `https://share.streamlit.io/arunangshu/hpmec/main/app.py`

### Step 3: Share Your App
Once deployed, you'll get a permanent URL like:
```
https://arunangshu-hpmec-app-[random-id].streamlit.app
```

## 📁 Project Structure

```
hpmec/
├── app.py                 # Main Streamlit GUI application
├── calculator.py          # Placeholder calculator module
├── requirements.txt       # Python dependencies
├── ethanol.yaml          # Sample OPLS-AA force field for ethanol
├── ethanol_sample.xyz    # Sample ethanol molecular geometry
├── idea.md               # Project blueprint and documentation
├── .streamlit/
│   └── config.toml       # Streamlit configuration
└── README.md             # This file
```

## 🔧 Usage

### Testing with Sample Files
The repository includes sample files for testing:
- **ethanol_sample.xyz**: Sample ethanol molecular geometry
- **ethanol.yaml**: OPLS-AA force field parameters for ethanol

Simply upload these files in the web interface to test the application!

### Single Molecule Calculation
1. Upload a `.xyz` molecular geometry file (or use `ethanol_sample.xyz`)
2. Upload a `.yaml` force field parameter file (or use `ethanol.yaml`)
3. View the interactive 3D molecule visualization
4. Download uploaded files if needed
5. Select "Single Molecule" mode
6. Click "Calculate Energy"

### Batch Processing (HPC-2)
1. Upload input files
2. Select "Batch Processing" mode
3. Choose number of calculations
4. Click "Calculate Energy"
5. View parallel performance metrics

## 📊 Features

### GUI Features
- ✅ File upload for `.xyz` and `.yaml` files
- ✅ **3D Interactive Molecule Visualization** with py3Dmol
- ✅ **Download buttons** for uploaded files
- ✅ Multiple visualization styles (Ball & Stick, Space Filling, Stick Only)
- ✅ Single and batch processing modes
- ✅ Real-time progress tracking
- ✅ Detailed energy breakdown
- ✅ Performance benchmarking
- ✅ Interactive documentation

### HPC Features
- ⚡ Multi-core parallelization using `multiprocessing.Pool`
- 📈 Linear speedup demonstration
- 🔍 Spatial tree optimization (cKDTree)
- 📊 Performance metrics and benchmarking

## 🧪 Example Files

Sample files are included in the repository:

**ethanol_sample.xyz** - Ethanol molecular geometry:
```
9
Ethanol molecule - Sample test file
C    0.000   0.000   0.000
C    1.520   0.000   0.000
O    2.020   1.350   0.000
H   -0.400   1.000   0.000
H   -0.400  -0.500   0.890
H   -0.400  -0.500  -0.890
H    1.920  -0.500   0.890
H    1.920  -0.500  -0.890
H    2.800   1.350   0.000
```

**ethanol.yaml** - OPLS-AA force field parameters (included in repo)

## 📦 Dependencies

Core requirements:
- `streamlit` - Web application framework
- `pandas` - Data manipulation
- `numpy` - Numerical computing
- `scipy` - Scientific computing (cKDTree)
- `pyyaml` - YAML parsing
- `py3Dmol` - 3D molecule visualization

Future integration (Modules 1-3):
- `rdkit` - Cheminformatics
- `xyz2mol` - Topology inference

## 🎓 Scientific Background

### Energy Function

$$E_{Total} = E_{Bonded} + E_{Non-Bonded}$$

**Bonded:**
- Bonds: $V_{bond}(b) = \sum \frac{1}{2} k_b (b - b_0)^2$
- Angles: $V_{angle}(\theta) = \sum \frac{1}{2} k_\theta (\theta - \theta_0)^2$
- Dihedrals: OPLS Fourier series

**Non-bonded:**
- Van der Waals: Lennard-Jones 12-6 potential
- Electrostatic: Coulomb's law

## 🔬 Performance Optimizations

### HPC-1: Spatial Tree Algorithm
- Replaces O(N²) pairwise loop with O(N log N) or O(N) tree search
- Uses `scipy.spatial.cKDTree` for neighbor finding
- Applies distance cutoff (typically 1.0 nm)

### HPC-2: Parallel Processing
- Uses `multiprocessing.Pool.map()` for distribution
- Bypasses Python's GIL for true parallelism
- Scales linearly with number of CPU cores

## 📝 Notes for Integration

The current `app.py` uses a placeholder `calculator.py` module. To integrate with real implementations:

1. **Import actual modules** (uncomment in `app.py`):
   ```python
   from calculator import calculate_single_molecule_energy
   from calculator import run_parallel_calculations
   ```

2. **Replace placeholder calls** with real function calls

3. **Update energy breakdown** with actual component values from Module 4

## 🐛 Troubleshooting

### Local Issues
- **Port already in use**: Use `streamlit run app.py --server.port 8502`
- **Module not found**: Run `pip install -r requirements.txt`

### Deployment Issues
- **Build fails**: Check `requirements.txt` for incompatible versions
- **App crashes**: Check Streamlit Cloud logs for errors
- **Slow startup**: Large dependencies may take time to install

## 📚 References

- OPLS-AA Force Field
- RDKit Cheminformatics Toolkit
- Streamlit Documentation: https://docs.streamlit.io/
- Python Multiprocessing: https://docs.python.org/3/library/multiprocessing.html

## 🤝 Contributing

This is a 4-member team project. For Module 5 (GUI & Parallelization):
- Maintain the Streamlit interface
- Optimize parallel processing performance
- Integrate with other modules
- Document deployment procedures

## 📄 License

This project is for educational purposes.

## 👤 Author

**Member 4** - GUI & Parallelization Lead
- Module 5 Implementation
- Streamlit Interface Design
- Multi-core Optimization
- System Integration

---

**Repository**: https://github.com/arunangshu/hpmec  
**Streamlit Cloud**: Deploy at https://share.streamlit.io/
