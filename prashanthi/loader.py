
import re
import numpy as np
from typing import Tuple, List

# -------- Public API --------
def load_structure(path: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Detect .xyz / .pdb / .mol (V2000) and return (symbols[N], positions[N,3]) in Å.
    Falls back to sniffing if extension is missing/misleading.
    """
    lower = path.lower()
    if lower.endswith(".xyz"):
        return read_xyz(path)
    if lower.endswith(".pdb"):
        return read_pdb(path)
    if lower.endswith(".mol") or lower.endswith(".sdf"):
        return read_mol_v2000(path)

    # Fallback: sniff by content
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        head = f.read(2048)
    if re.search(r"\bV2000\b", head):
        return read_mol_v2000(path)
    if ("ATOM" in head) or ("HETATM" in head):
        return read_pdb(path)
    return read_xyz(path)

# -------- Parsers --------
def read_xyz(path: str) -> Tuple[np.ndarray, np.ndarray]:
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        first = f.readline().strip()
        symbols, pos = [], []
        try:
            n = int(first)
            _ = f.readline()  # comment
            for _ in range(n):
                parts = f.readline().split()
                if len(parts) < 4:  # skip blank/short lines
                    continue
                symbols.append(parts[0])
                pos.append([float(parts[1]), float(parts[2]), float(parts[3])])
        except ValueError:
            # No count header: treat whole file as coordinates
            rest = [first] + f.readlines()
            for ln in rest:
                parts = ln.split()
                if len(parts) >= 4:
                    symbols.append(parts[0])
                    pos.append([float(parts[1]), float(parts[2]), float(parts[3])])
    if not symbols:
        raise ValueError("XYZ parse error: no atoms found.")
    return np.array(symbols), np.asarray(pos, float)

def read_pdb(path: str) -> Tuple[np.ndarray, np.ndarray]:
    symbols: List[str] = []
    coords: List[List[float]] = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            rec = line[:6].strip()
            if rec in ("ATOM", "HETATM"):
                # Try fixed columns first, then fallback to split.
                try:
                    x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                except ValueError:
                    parts = line.split()
                    # Last 6 fields often include occupancy/temp + element; be lenient:
                    floats = [p for p in parts if _is_float(p)]
                    if len(floats) < 3:
                        raise
                    x, y, z = map(float, floats[-3:])
                elem = line[76:78].strip()
                if not elem:
                    name = line[12:16].strip()
                    elem = re.sub(r"[^A-Za-z]", "", name)[:2].title() or "X"
                symbols.append(elem)
                coords.append([x, y, z])
    if not symbols:
        raise ValueError("PDB parse error: no ATOM/HETATM records found.")
    return np.array(symbols), np.asarray(coords, float)

def read_mol_v2000(path: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Robust MOL V2000 parser:
    - Locates the counts line with 'V2000'
    - Parses atom block by tokens (split), not fixed columns
    - Accepts squished spacing like '0.0000   0.0000   0.0000 C'
    """
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        lines = [ln.rstrip("\n\r") for ln in f]

    # Find counts line (usually line index 3), but scan a bit to be safe.
    counts_i = None
    for i in range(0, min(12, len(lines))):
        if re.search(r"\bV2000\b", lines[i]):
            counts_i = i
            break
    if counts_i is None:
        raise ValueError("MOL V2000 parse error: counts line with 'V2000' not found.")

    # Parse counts line by tokens (handles variable spacing)
    parts = lines[counts_i].split()
    # First two tokens must be n_atoms and n_bonds (per V2000)
    try:
        n_atoms = int(parts[0]); n_bonds = int(parts[1])
    except Exception as e:
        # Some generators pad fields; try fixed-field fallback (0:3, 3:6)
        try:
            n_atoms = int(lines[counts_i][0:3])
            n_bonds = int(lines[counts_i][3:6])
        except Exception:
            raise ValueError(f"MOL V2000 counts parse error: {lines[counts_i]!r}") from e

    # Atom block
    atom_start = counts_i + 1
    atom_end = atom_start + n_atoms
    if atom_end > len(lines):
        raise ValueError("MOL V2000 atom block exceeds file length.")

    symbols: List[str] = []
    coords: List[List[float]] = []
    for ln in lines[atom_start:atom_end]:
        if not ln.strip():
            continue
        toks = ln.split()
        # Expect: x y z sym ...
        if len(toks) < 4:
            # fallback: try fixed columns if super-cramped
            x = float(ln[0:10]); y = float(ln[10:20]); z = float(ln[20:30])
            sym = ln[31:34].strip().title() or "X"
        else:
            # Be tolerant: element could be toks[3]; if it's numeric (weird), try toks[4]
            x, y, z = float(toks[0]), float(toks[1]), float(toks[2])
            sym = toks[3].title()
            if _is_float(sym) and len(toks) >= 5:
                sym = toks[4].title()
            # Normalize common multi-letter elements (Cl, Br) capitalization
            if len(sym) >= 2:
                sym = sym[0].upper() + sym[1:].lower()
        symbols.append(sym)
        coords.append([x, y, z])

    if not symbols or len(symbols) != n_atoms:
        raise ValueError("MOL V2000 parse error: atom count mismatch.")

    return np.array(symbols), np.asarray(coords, float)

# -------- Helpers --------
def _is_float(s: str) -> bool:
    try:
        float(s)
        return True
    except Exception:
        return False

# -------- CLI quick test --------
if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python loader.py file.xyz|file.pdb|file.mol")
        raise SystemExit(1)
    syms, xyz = load_structure(sys.argv[1])
    print(f"Loaded N={len(syms)} atoms")
    print("Symbols:", syms.tolist())
    print("First 3 coords:\n", xyz[:3])
