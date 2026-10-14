"""Independent checks of a straight-slit design.

Usage: python bloch_verify.py R t2 angle3 t3 [angle4 t4 ...]   (angles in radians)

Checks
  * inradius by three routes: the fast routine, the curved-capable routine,
    and brute force (exact distance to every slit from tens of millions of
    random points, followed by a local polish of the best point);
  * conformal radius at several resolutions of the spectral solver, and by an
    independent finite-difference-free route: the same solver on the design
    rotated by a random angle (the value must not change).
"""
import os as _os
import sys as _sys

_sys.path.append(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # shared modules in the parent folder

import sys

import numpy as np
from scipy.optimize import minimize

import bloch_curved as bc
import bloch_design as d
from bloch_solver import Slit, conformal_radius

RECORD = 0.6563937


def main():
    p = [float(v) for v in sys.argv[1:]]
    big_r = p[0]
    classes = [(0.0, 1.0), (np.pi / 3, p[1])] + [(p[2 + 2 * k], p[3 + 2 * k]) for k in range((len(p) - 2) // 2)]
    print(f"R = {big_r}")
    for alpha, tip in classes:
        print(f"   slits at +-{np.degrees(alpha):.5f} deg (and rotations by 120 deg), from radius {tip} to the circle")
    segs = d.segments(big_r, classes)
    print(f"number of slits: {len(segs)}")
    u = np.array([s[0] for s in segs])
    w = np.array([s[1] - s[0] for s in segs])

    r_fast = d.inradius(big_r, classes)
    r_fine = d.inradius(big_r, classes, grid=(900, 451))
    r_curv = bc.inradius(big_r, [(a, t, 0.0) for a, t in classes], grid=(400, 201), per_slit=800)
    rng = np.random.default_rng(12345)
    best, where = 0.0, None
    for _ in range(100):  # 30 million points over the whole disc, no symmetry assumed
        z = np.sqrt(rng.random(300000)) * big_r * np.exp(2j * np.pi * rng.random(300000))
        m = d.feature_distances(z, big_r, u, w)[0].min(axis=1)
        i = m.argmax()
        if m[i] > best:
            best, where = float(m[i]), z[i]
    pol = minimize(lambda q: -d.feature_distances(np.array([complex(q[0], q[1])]), big_r, u, w)[0].min(),
                   [where.real, where.imag], method="Nelder-Mead", options={"xatol": 1e-12, "fatol": 1e-14, "maxiter": 4000})
    r_brute = max(best, -pol.fun)
    print(f"inradius: fast {r_fast:.10f} | fine grid {r_fine:.10f} | curved routine {r_curv:.10f} | brute force + polish {r_brute:.10f}")
    inr = max(r_fast, r_fine, r_curv, r_brute)

    slits = d.slit_list(big_r, classes)
    vals = {n: conformal_radius(big_r, slits, n) for n in (24, 40, 64, 96, 128)}
    print("conformal radius by resolution:", {n: f"{v:.12f}" for n, v in vals.items()})
    rot = 0.7391
    rotated = [Slit(s.alpha + rot, s.a, s.c) for s in slits]
    print(f"conformal radius of the rotated design: {conformal_radius(big_r, rotated, 96):.12f}")
    conf = min(vals[96], vals[128])
    print(f"\nratio = largest inradius / smallest conformal radius = {inr:.10f} / {conf:.10f} = {inr / conf:.9f}")
    print(f"record 0.6563937; difference {inr / conf - RECORD:+.3e}  ->  {'BELOW the record' if inr / conf < RECORD else 'not below the record'}")


if __name__ == "__main__":
    main()
