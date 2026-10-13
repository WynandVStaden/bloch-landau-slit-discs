"""Independent checks of a design with curved slits.

Usage: python bloch_verify_curved.py R t2 [angle tip bend]...   (angles in radians)
"""
import sys

import numpy as np
from scipy.optimize import minimize
from scipy.spatial import cKDTree

import bloch_curved as bc
import bloch_residual as br
from bloch_solver import Slit, conformal_radius

RECORD = 0.6563937


def main():
    p = [float(v) for v in sys.argv[1:]]
    big_r = p[0]
    classes = [(0.0, 1.0, 0.0), (np.pi / 3, p[1], 0.0)] + [tuple(p[2 + 3 * k:5 + 3 * k]) for k in range((len(p) - 2) // 3)]
    ang, tip, bend = bc.expand(classes)
    a = np.log(big_r / tip)
    print(f"R = {big_r}; {len(ang)} slits")
    for al, tp, bd in classes:
        aa = np.log(big_r / tp)
        print(f"   leaves the circle at +-{np.degrees(al):.4f} deg, bend {bd:+.5f}, tip at radius {tp} and angle {np.degrees(al + bd * aa * aa):.4f} deg")

    # do any two slits cross or touch? (they must not: the model assumes disjoint slits)
    dense = bc.curve(-np.linspace(0, 1, 4000)[None, :], big_r, ang[:, None], a[:, None], bend[:, None])
    gap = np.inf
    for i in range(len(ang)):
        tree = cKDTree(np.column_stack([dense[i].real, dense[i].imag]))
        for j in range(i + 1, len(ang)):
            gap = min(gap, tree.query(np.column_stack([dense[j].real, dense[j].imag]))[0].min())
    print(f"smallest distance between two different slits: {gap:.5f}")

    r_fast = bc.inradius(big_r, classes)
    r_fine = bc.inradius(big_r, classes, grid=(500, 251), per_slit=1000)
    pts = bc.curve(-np.linspace(0, 1, 30000)[None, :], big_r, ang[:, None], a[:, None], bend[:, None]).ravel()
    tree = cKDTree(np.column_stack([pts.real, pts.imag]))
    rng = np.random.default_rng(2026)
    best, where = 0.0, None
    for _ in range(60):  # 30 million points over the whole disc, no symmetry assumed
        z = np.sqrt(rng.random(500000)) * big_r * np.exp(2j * np.pi * rng.random(500000))
        m = np.minimum(tree.query(np.column_stack([z.real, z.imag]))[0], big_r - np.abs(z))
        i = m.argmax()
        if m[i] > best:
            best, where = float(m[i]), z[i]

    def clear(q):
        return min(tree.query([q[0], q[1]])[0], big_r - np.hypot(q[0], q[1]))

    pol = minimize(lambda q: -clear(q), [where.real, where.imag], method="Nelder-Mead", options={"xatol": 1e-10, "fatol": 1e-12, "maxiter": 3000})
    r_brute = max(best, -pol.fun)
    print(f"inradius: routine {r_fast:.9f} | fine {r_fine:.9f} | brute force on 30000-point slits + polish {r_brute:.9f} (at {where:.4f})")
    inr = max(r_fast, r_fine, r_brute)

    slits = bc.slit_objects(big_r, classes)
    vals = {n: conformal_radius(big_r, slits, n) for n in (40, 64, 96, 128)}
    print("conformal radius by resolution:", {n: f"{v:.12f}" for n, v in vals.items()})
    rotated = [Slit(s.alpha + 0.7391, s.a, s.c) for s in slits]
    print(f"conformal radius of the rotated design: {conformal_radius(big_r, rotated, 96):.12f}")
    xe = np.cos(np.linspace(0.003, np.pi - 0.003, 1501))
    psi, rows, log_rho = br.solve(big_r, slits, 96)
    worst = max(float(np.abs(-s.a * xe + rows(i, xe) @ psi).max()) for i, s in enumerate(slits))
    print(f"largest |Green function| on the slits at test points (n=96): {worst:.2e}")
    conf = min(vals[96], vals[128])
    print(f"\nratio = {inr:.9f} / {conf:.9f} = {inr / conf:.9f};  record {RECORD};  difference {inr / conf - RECORD:+.3e}")


if __name__ == "__main__":
    main()
