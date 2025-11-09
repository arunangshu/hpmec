
import sys, time
import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple, List
from loader import load_structure

RC, SKIN = 2.5, 0.3  

def minimum_image(dr: np.ndarray, box: np.ndarray) -> np.ndarray:
    return dr - np.round(dr / box) * box

def make_box(positions: np.ndarray, pad: float = 1.5) -> Tuple[np.ndarray, np.ndarray]:
    lo = positions.min(axis=0); hi = positions.max(axis=0)
    L = float(max(hi - lo) + 2 * pad)
    center = (hi + lo) / 2.0
    return positions - center, np.array([L, L, L])

class VerletList:
    def __init__(self, rc: float, skin: float, box: np.ndarray):
        self.rc, self.skin = rc, skin
        self.rc_skin2 = (rc + skin) ** 2
        self.box = box.copy()
        self.list: List[np.ndarray] = []

    def build(self, pos: np.ndarray):
        n = len(pos)
        nbrs = [[] for _ in range(n)]
        for i in range(n-1):
            dr = pos[i+1:] - pos[i]
            dr = minimum_image(dr, self.box)
            d2 = np.einsum("ij,ij->i", dr, dr)
            js = np.where(d2 <= self.rc_skin2)[0] + (i + 1)
            for j in js:
                nbrs[i].append(j); nbrs[j].append(i)
        self.list = [np.array(lst, dtype=int) for lst in nbrs]

    def pairs_within_rc(self, pos: np.ndarray) -> List[tuple]:
        rc2 = self.rc * self.rc
        pairs = []
        for i, neigh in enumerate(self.list):
            if len(neigh) == 0: continue
            rij = pos[neigh] - pos[i]
            rij = minimum_image(rij, self.box)
            d2 = np.einsum("ij,ij->i", rij, rij)
            good = np.where(d2 <= rc2)[0]
            for idx in good:
                j = neigh[idx]
                if j > i: pairs.append((i, j))
        return pairs

def compute_forces_verlet(pos: np.ndarray, box: np.ndarray, vlist: VerletList) -> np.ndarray:
    forces = np.zeros_like(pos)
    pairs = vlist.pairs_within_rc(pos)
    for i, j in pairs:
        rij = pos[j] - pos[i]
        rij = minimum_image(rij, box)
        r2 = float(rij @ rij)
        if r2 == 0: continue
        inv_r2 = 1.0 / r2
        sr2 = inv_r2
        sr6 = sr2 * sr2 * sr2
        sr12 = sr6 * sr6
        f_over_r = 24.0 * (2.0 * sr12 - sr6) * inv_r2
        fij = f_over_r * rij
        forces[i] += fij; forces[j] -= fij
    return forces

def visualize_static(symbols, pos, box, bond_cutoff=1.8):
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection="3d")
    L = box[0]
    ax.set_xlim(-L/2, L/2); ax.set_ylim(-L/2, L/2); ax.set_zlim(-L/2, L/2)
    ax.set_box_aspect((1,1,1))
    ax.set_title("3D Visualization (Verlet list)")
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
        print("Usage: python verlet.py file.xyz|file.pdb|file.mol [bond_cutoff=1.8]")
        sys.exit(1)
    path = sys.argv[1]
    bond_cutoff = float(sys.argv[2]) if len(sys.argv) >= 3 else 1.8

    symbols, pos0 = load_structure(path)
    pos, box = make_box(pos0, pad=2.0)

    v = VerletList(RC, SKIN, box)
    t0 = time.perf_counter(); v.build(pos); build_time = time.perf_counter() - t0
    t1 = time.perf_counter(); _ = compute_forces_verlet(pos, box, v); force_time = time.perf_counter() - t1

    print(f"[Verlet] N={len(pos)} | build={build_time:.6f}s | force={force_time:.6f}s")
    visualize_static(symbols, pos, box, bond_cutoff=bond_cutoff)
