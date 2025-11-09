
import sys, time
import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple, Dict, List
from itertools import product
from loader import load_structure

RC = 2.5  

def minimum_image(dr: np.ndarray, box: np.ndarray) -> np.ndarray:
    return dr - np.round(dr / box) * box

def make_box(positions: np.ndarray, pad: float = 1.5) -> Tuple[np.ndarray, np.ndarray]:
    lo = positions.min(axis=0); hi = positions.max(axis=0)
    L = float(max(hi - lo) + 2 * pad)
    center = (hi + lo) / 2.0
    return positions - center, np.array([L, L, L])

def build_cells(pos: np.ndarray, box: np.ndarray, rc: float = RC):
    ncell = np.floor(box / rc).astype(int)
    ncell = np.maximum(ncell, 1)
    cell_size = box / ncell
    shifted = (pos + 0.5 * box) % box
    idxs = np.floor(shifted / cell_size).astype(int)
    cells: Dict[tuple, List[int]] = {}
    for p, (ix, iy, iz) in enumerate(idxs):
        cells.setdefault((ix, iy, iz), []).append(p)
    neighbor_offsets = list(product([-1, 0, 1], repeat=3))
    return cells, ncell, neighbor_offsets

def compute_forces_celllist(pos: np.ndarray, box: np.ndarray, rc: float = RC) -> np.ndarray:
    forces = np.zeros_like(pos); rc2 = rc * rc
    cells, ncell, neighbor_offsets = build_cells(pos, box, rc)

    def wrap(t, n): return (t + n) % n

    for (ix, iy, iz), plist in cells.items():
        for dx, dy, dz in neighbor_offsets:
            jx, jy, jz = wrap(ix+dx, ncell[0]), wrap(iy+dy, ncell[1]), wrap(iz+dz, ncell[2])
            qlist = cells.get((jx, jy, jz), [])
            for a in plist:
                for b in qlist:
                    if (ix, iy, iz) == (jx, jy, jz) and b <= a:
                        continue
                    rij = pos[b] - pos[a]
                    rij = minimum_image(rij, box)
                    r2 = float(rij @ rij)
                    if 0.0 < r2 <= rc2:
                        inv_r2 = 1.0 / r2
                        sr2 = inv_r2
                        sr6 = sr2 * sr2 * sr2
                        sr12 = sr6 * sr6
                        f_over_r = 24.0 * (2.0 * sr12 - sr6) * inv_r2
                        fij = f_over_r * rij
                        forces[a] += fij; forces[b] -= fij
    return forces

def visualize_static(symbols, pos, box, bond_cutoff=1.8):
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection="3d")
    L = box[0]
    ax.set_xlim(-L/2, L/2); ax.set_ylim(-L/2, L/2); ax.set_zlim(-L/2, L/2)
    ax.set_box_aspect((1,1,1))
    ax.set_title("3D Visualization (Cell-list)")
    colors = {"H":"#FFFFFF","C":"#444444","N":"#3050F8","O":"#FF0D0D","S":"#FFFF30"}
    c = [colors.get(s, "#AAAAAA") for s in symbols]
    sizes = np.array([80 if s=="H" else 160 for s in symbols])
    ax.scatter(pos[:,0], pos[:,1], pos[:,2], s=sizes, c=c, depthshade=False, edgecolor="k", linewidths=0.3)
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
        print("Usage: python celllist.py file.xyz|file.pdb|file.mol [bond_cutoff=1.8]")
        sys.exit(1)
    path = sys.argv[1]
    bond_cutoff = float(sys.argv[2]) if len(sys.argv) >= 3 else 1.8
    symbols, pos0 = load_structure(path)
    pos, box = make_box(pos0, pad=2.0)

    t0 = time.perf_counter(); cells, ncell, neigh = build_cells(pos, box, rc=RC); build_time = time.perf_counter() - t0
    t1 = time.perf_counter(); _ = compute_forces_celllist(pos, box, rc=RC); force_time = time.perf_counter() - t1

    print(f"[CellList] N={len(pos)} | build={build_time:.6f}s | force={force_time:.6f}s")
    visualize_static(symbols, pos, box, bond_cutoff=bond_cutoff)
