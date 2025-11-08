#!/usr/bin/env python3
import sys, math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

# ----------------------------
# Config (tweak as you like)
# ----------------------------
DT = 0.002            # time step
STEPS_PER_FRAME = 2   # sim steps per animation frame
RC = 2.5              # neighbor cutoff (sigma units)
SKIN = 0.3            # Verlet skin distance
EPSILON = 1.0         # LJ epsilon
SIGMA = 1.0           # LJ sigma
MASS = 1.0            # particle mass
DAMPING = 0.0         # simple vel damping (0 = none)
MAX_V0 = 0.2          # initial random velocity scale
# ----------------------------

def read_xyz(path):
    with open(path, "r") as f:
        n = int(f.readline().strip())
        comment = f.readline()
        symbols, pos = [], []
        for _ in range(n):
            parts = f.readline().split()
            symbols.append(parts[0])
            pos.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return np.array(symbols), np.asarray(pos, dtype=float)

def make_box(positions, pad=1.5):
    # Build a cubic periodic box slightly larger than the data span.
    lo = positions.min(axis=0)
    hi = positions.max(axis=0)
    span = hi - lo
    L = float(max(span) + 2 * pad)
    center = (hi + lo) / 2.0
    # shift positions to center the box at origin
    shifted = positions - center
    return shifted, np.array([L, L, L])

def pbc_wrap(r, box):
    # Wrap positions into [-L/2, L/2)
    return (r + 0.5*box) % box - 0.5*box

def minimum_image(dr, box):
    # Apply minimum-image convention
    return dr - np.round(dr / box) * box

class VerletList:
    def __init__(self, rc, skin, box):
        self.rc = rc
        self.skin = skin
        self.rc_skin = rc + skin
        self.rc_skin2 = (rc + skin) ** 2
        self.box = box.copy()
        self.ref_pos = None
        self.list = None

    def needs_rebuild(self, pos):
        if self.ref_pos is None:
            return True
        # If any particle moved more than skin/2 since last build -> rebuild
        disp = pos - self.ref_pos
        disp = minimum_image(disp, self.box)
        max_disp = np.max(np.linalg.norm(disp, axis=1))
        return max_disp > self.skin * 0.5

    def build(self, pos):
        n = len(pos)
        nbrs = [[] for _ in range(n)]
        # O(N^2) build is okay for small N; for big N swap to cell-lists
        for i in range(n-1):
            dr = pos[i+1:] - pos[i]
            dr = minimum_image(dr, self.box)
            dist2 = np.einsum("ij,ij->i", dr, dr)
            mask = dist2 <= self.rc_skin2
            js = np.where(mask)[0] + (i+1)
            for j in js:
                nbrs[i].append(j)
                nbrs[j].append(i)
        self.list = [np.array(lst, dtype=int) for lst in nbrs]
        self.ref_pos = pos.copy()

    def pairs_within_rc(self, pos):
        # Iterate pairs from cached list but filter by real rc
        rc2 = self.rc * self.rc
        pairs = []
        for i, neigh in enumerate(self.list):
            if len(neigh) == 0:
                continue
            rij = pos[neigh] - pos[i]
            rij = minimum_image(rij, self.box)
            d2 = np.einsum("ij,ij->i", rij, rij)
            good = np.where(d2 <= rc2)[0]
            for idx in good:
                j = neigh[idx]
                if j > i:
                    pairs.append((i, j))
        return pairs

def lj_force(r2, epsilon=EPSILON, sigma=SIGMA):
    # Returns scalar f_over_r = |F|/r and potential (optional)
    inv_r2 = 1.0 / r2
    sr2 = (sigma * sigma) * inv_r2
    sr6 = sr2 * sr2 * sr2
    sr12 = sr6 * sr6
    # F = 24*epsilon*(2*sr12 - sr6) * (1/r^2) * r_vec
    f_over_r = 24.0 * epsilon * (2.0 * sr12 - sr6) * inv_r2
    return f_over_r

def compute_forces(pos, box, vlist):
    n = len(pos)
    forces = np.zeros_like(pos)
    pairs = vlist.pairs_within_rc(pos)
    for i, j in pairs:
        rij = pos[j] - pos[i]
        rij = minimum_image(rij, box)
        r2 = np.dot(rij, rij)
        if r2 == 0:
            continue
        f_over_r = lj_force(r2)
        fij = f_over_r * rij
        forces[i] += fij
        forces[j] -= fij
    return forces

def simulate(path):
    symbols, pos0 = read_xyz(path)
    pos, box = make_box(pos0)
    n = len(pos)
    vel = (np.random.rand(n, 3) - 0.5) * 2 * MAX_V0
    acc = np.zeros_like(pos)

    vlist = VerletList(RC, SKIN, box)
    vlist.build(pos)

    # Matplotlib figure
    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection='3d')
    scat = ax.scatter(pos[:, 0], pos[:, 1], pos[:, 2], s=12)
    L = box[0]
    ax.set_xlim(-L/2, L/2)
    ax.set_ylim(-L/2, L/2)
    ax.set_zlim(-L/2, L/2)
    ax.set_box_aspect((1,1,1))
    ax.set_title("XYZ Interactive LJ Simulation (Verlet list)")
    txt = ax.text2D(0.02, 0.95, "", transform=ax.transAxes)

    paused = {"v": False}
    def on_key(event):
        if event.key == " ":
            paused["v"] = not paused["v"]
        elif event.key == "r":
            # rebuild neighbor list manually
            vlist.build(pos)
        elif event.key == "q":
            plt.close(fig)
    fig.canvas.mpl_connect("key_press_event", on_key)

    # Initial force (velocity-Verlet)
    if vlist.needs_rebuild(pos):
        vlist.build(pos)
    acc[:] = compute_forces(pos, box, vlist) / MASS

    def step():
        nonlocal pos, vel, acc
        # velocity-Verlet
        vel *= (1.0 - DAMPING)
        pos = pos + vel * DT + 0.5 * acc * DT * DT
        pos = pbc_wrap(pos, box)
        if vlist.needs_rebuild(pos):
            vlist.build(pos)
        new_acc = compute_forces(pos, box, vlist) / MASS
        vel = vel + 0.5 * (acc + new_acc) * DT
        acc = new_acc

    def update(_frame):
        if not paused["v"]:
            for _ in range(STEPS_PER_FRAME):
                step()
        scat._offsets3d = (pos[:,0], pos[:,1], pos[:,2])
        txt.set_text(f"N={len(pos)}  rc={RC}  skin={SKIN}  (space=pause, r=rebuild, q=quit)")
        return scat, txt

    anim = FuncAnimation(fig, update, interval=20, blit=False)
    plt.show()

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("xyz")
    p.add_argument("--static", action="store_true", help="Show a static molecule view (no dynamics)")
    p.add_argument("--bond-cutoff", type=float, default=1.8, help="Å cutoff for drawing bonds in static mode")
    args = p.parse_args()

    if not args.static:
        simulate(args.xyz)
        sys.exit(0)

    # --- Static viewer (no simulation) ---
    symbols, pos0 = read_xyz(args.xyz)
    pos, box = make_box(pos0, pad=2.0)

    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection='3d')
    L = box[0]
    ax.set_xlim(-L/2, L/2); ax.set_ylim(-L/2, L/2); ax.set_zlim(-L/2, L/2)
    ax.set_box_aspect((1,1,1))
    ax.set_title("Molecule Visualization using Verlet List")

    # atom colors by element (very small set; fallback gray)
    color_map = {
        "H":"#FFFFFF","C":"#444444","N":"#3050F8","O":"#FF0D0D","S":"#FFFF30","P":"#FF8000",
        "F":"#90E050","Cl":"#1FF01F","Br":"#A62929","I":"#940094"
    }
    colors = [color_map.get(sym, "#AAAAAA") for sym in symbols]
    sizes = np.array([80 if s=="H" else 160 for s in symbols])  # crude size hint

    ax.scatter(pos[:,0], pos[:,1], pos[:,2], s=sizes, c=colors, depthshade=False, edgecolor="k", linewidths=0.3)

    # optional bonds by simple distance cutoff (no PBC; for molecules this is fine)
    cutoff2 = args.bond_cutoff**2
    n = len(pos)
    for i in range(n-1):
        rij = pos[i+1:] - pos[i]
        d2 = np.einsum("ij,ij->i", rij, rij)
        js = np.where((d2>1e-6) & (d2<=cutoff2))[0] + (i+1)
        for j in js:
            xs = [pos[i,0], pos[j,0]]
            ys = [pos[i,1], pos[j,1]]
            zs = [pos[i,2], pos[j,2]]
            ax.plot(xs, ys, zs, linewidth=1.2, alpha=0.8)

    plt.show()

