"""Full independent checks of a design in the flexible curved family (bloch_curved2.py), and a picture of it.

Usage: python bloch_verify_curved2.py <params as printed by the search>
"""
import sys

import numpy as np
from scipy.optimize import minimize
from scipy.spatial import cKDTree

import bloch_curved2 as b
import bloch_residual as br
from bloch_solver import conformal_radius

PUBLISHED = 0.656319277272  # Ortega-Cerda and Pridhnani (numerical)
PROVED = 0.6563937  # Carroll and Ortega-Cerda


def main():
    p = np.array([float(v) for v in sys.argv[1:]])
    big_r, classes = b.unpack(p)
    ang, tip, b1, b2 = b.expand(classes)
    a = np.log(big_r / tip)
    print(f"R = {big_r}; {len(ang)} slits; crossing test: {b.slits_cross(big_r, ang, tip, b1, b2)}", flush=True)
    for al, tp, u, v in classes:
        s = np.log(tp / big_r)
        print(f"   leaves the circle at +-{np.degrees(al):.4f} deg, tip at radius {tp:.6f} and angle {np.degrees(al + u * s**2 + v * s**4):.4f} deg, bends {u:+.5f}, {v:+.5f}")

    dense = b.curve(-np.linspace(0, 1, 4000)[None, :], big_r, ang[:, None], a[:, None], b1[:, None], b2[:, None])
    gap = np.inf
    for i in range(len(ang)):
        tree = cKDTree(np.column_stack([dense[i].real, dense[i].imag]))
        for j in range(i + 1, len(ang)):
            gap = min(gap, tree.query(np.column_stack([dense[j].real, dense[j].imag]))[0].min())
    print(f"smallest distance between two different slits: {gap:.5f}", flush=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.2, 6.2))
    t = np.linspace(0, 2 * np.pi, 720)
    ax.plot(big_r * np.cos(t), big_r * np.sin(t), color="black", lw=1.2)
    for k in range(len(ang)):
        ax.plot(dense[k].real, dense[k].imag, color="black", lw=1.2)
    ax.add_patch(matplotlib.patches.Circle((0, 0), 1.0, fill=False, color="#c44e52", lw=1.0, ls="--"))
    ax.set_aspect("equal")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig("domain_G.png", dpi=150)

    r_fast = b.inradius(big_r, classes)
    r_fine = b.inradius(big_r, classes, grid=(500, 251), per_slit=1000)
    print(f"inradius: routine {r_fast:.10f} | fine {r_fine:.10f}", flush=True)
    pts = b.curve(-np.linspace(0, 1, 12000)[None, :], big_r, ang[:, None], a[:, None], b1[:, None], b2[:, None]).ravel()
    tree = cKDTree(np.column_stack([pts.real, pts.imag]))
    rng = np.random.default_rng(2026)
    best, where = 0.0, None
    for k in range(60):  # 30 million points over the whole disc, no symmetry assumed
        z = np.sqrt(rng.random(500000)) * big_r * np.exp(2j * np.pi * rng.random(500000))
        m = np.minimum(tree.query(np.column_stack([z.real, z.imag]), workers=-1)[0], big_r - np.abs(z))
        i = m.argmax()
        if m[i] > best:
            best, where = float(m[i]), z[i]
        if k % 15 == 14:
            print(f"   brute force: {(k + 1) * 0.5:.0f} million points, largest clearance so far {best:.7f}", flush=True)
    pol = minimize(lambda q: -min(tree.query([q[0], q[1]])[0], big_r - np.hypot(q[0], q[1])), [where.real, where.imag],
                   method="Nelder-Mead", options={"xatol": 1e-10, "fatol": 1e-12, "maxiter": 3000})
    r_brute = max(best, -pol.fun)
    print(f"inradius by brute force + polish: {r_brute:.9f} (at {where:.4f})", flush=True)
    inr = max(r_fast, r_fine, r_brute)

    slits = b.slit_objects(big_r, classes)
    vals = {n: conformal_radius(big_r, slits, n) for n in (40, 64, 96, 128)}
    print("conformal radius by resolution:", {n: f"{v:.12f}" for n, v in vals.items()}, flush=True)
    rotated = [b.Slit2(s.alpha + 0.7391, s.a, s.b1, s.b2) for s in slits]
    print(f"conformal radius of the rotated design: {conformal_radius(big_r, rotated, 96):.12f}", flush=True)
    xe = np.cos(np.linspace(0.003, np.pi - 0.003, 1501))
    psi, rows, _ = br.solve(big_r, slits, 96)
    worst = max(float(np.abs(-s.a * xe + rows(i, xe) @ psi).max()) for i, s in enumerate(slits))
    print(f"largest |Green function| on the slits at test points (n=96): {worst:.2e}", flush=True)
    conf = min(vals[96], vals[128])
    ratio = inr / conf
    print(f"\nratio = {inr:.9f} / {conf:.9f} = {ratio:.9f}")
    print(f"below the best published numerical value {PUBLISHED} by {PUBLISHED - ratio:+.3e}; below the best proved bound {PROVED} by {PROVED - ratio:+.3e}")


if __name__ == "__main__":
    main()
