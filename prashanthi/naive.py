
import sys, time
import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple
from loader import load_structure

RC = 2.5  

def minimum_image(dr: np.ndarray, box: np.ndarray) -> np.ndarray:
    return dr - np.round(dr / box) * box

def make_box(positions: np.ndarray, pad: float = 1.5) -> Tuple[np.ndarray, np.ndarray]:
    lo = positions.min(axis=0); hi = positions.max(axis=0)
    L = float(max(hi - lo) + 2 * pad)
    center = (hi + lo) / 2.0
    return positions - center, np.array([L, L, L])

def compute_forces_naive(pos: np.ndarray, box: np.ndarray, rc: float = RC) -> np.ndarray:
    n = len(pos); forces = np.zeros_like(pos); rc2 = rc * rc
    for i in range(n - 1):
        rij = pos[i+1:] - pos[i]
        rij = minimum_image(rij, box)
        r2 = np.einsum("ij,ij->i", rij, rij)
        mask = r2 <= rc2
        idxs = np.where(mask)[0]
        for k in idxs:
            j = i + 1 + k
            r2_ij = r2[k]
            if r2_ij == 0: continue
            inv_r2 = 1.0 / r2_ij
            sr2 = inv_r2
            sr6 = sr2 * sr2 * sr2
            sr12 = sr6 * sr6
            f_over_r = 24.0 * (2.0 * sr12 - sr6) * inv_r2
            fij = f_over_r * rij[k]
            forces[i] += fij; forces[j] -= fij
    return forces

def visualize_static(symbols: np.ndarray, pos: np.ndarray, box: np.ndarray, bond_cutoff: float = 1.8):
    from mpl_toolkits.mplot3d import Axes3D  
    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection="3d")
    L = box[0]
    ax.set_xlim(-L/2, L/2); ax.set_ylim(-L/2, L/2); ax.set_zlim(-L/2, L/2)
    ax.set_box_aspect((1,1,1))
    ax.set_title("3D Visualization (Naive Algorithm)")
    colors = {"H":"#FFFFFF","C":"#444444","N":"#3050F8","O":"#FF0D0D","S":"#FFFF30"}
    c = [colors.get(s, "#AAAAAA") for s in symbols]
    sizes = np.array([80 if s=="H" else 160 for s in symbols])
    ax.scatter(pos[:,0], pos[:,1], pos[:,2], s=sizes, c=c, depthshade=False, edgecolor="k", linewidths=0.3)
    # bonds (no PBC for display)
    n = len(pos); cutoff2 = bond_cutoff**2
    for i in range(n-1):
        rij = pos[i+1:] - pos[i]
        d2 = np.einsum("ij,ij->i", rij, rij)
        js = np.where((d2>1e-6) & (d2<=cutoff2))[0] + (i+1)
        for j in js:
            ax.plot([pos[i,0], pos[j,0]],[pos[i,1], pos[j,1]],[pos[i,2], pos[j,2]], linewidth=1.2, alpha=0.8)
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python naive.py file.xyz|file.pdb|file.mol [bond_cutoff=1.8]")
        sys.exit(1)
    path = sys.argv[1]
    bond_cutoff = float(sys.argv[2]) if len(sys.argv) >= 3 else 1.8
    symbols, pos0 = load_structure(path)
    pos, box = make_box(pos0, pad=2.0)

    build_time = 0.0
    t0 = time.perf_counter()
    _ = compute_forces_naive(pos, box, rc=RC)
    force_time = time.perf_counter() - t0

    print(f"[Naive] N={len(pos)} | build=0.000s | force={force_time:.6f}s")
    visualize_static(symbols, pos, box, bond_cutoff=bond_cutoff)
