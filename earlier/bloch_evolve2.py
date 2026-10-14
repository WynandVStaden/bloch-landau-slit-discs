"""Differential evolution over designs with curved slits (model in bloch_curved.py).

Parameters: R, tip radius of the straight 60-degree slits, then (angle, tip
radius, bend) for each extra class. The slits at 0 and 60 degrees lie on
symmetry axes and stay straight.

Usage: python bloch_evolve2.py n_extra_classes [seed] [maxiter] [start values...]
"""
import os as _os
import sys as _sys

_sys.path.append(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # shared modules in the parent folder

import os
import sys

import numpy as np
from scipy.optimize import differential_evolution, minimize

import bloch_curved as bc
from bloch_solver import conformal_radius


def unpack(p):
    classes = [(0.0, 1.0, 0.0), (np.pi / 3, p[1], 0.0)]
    for k in range((len(p) - 2) // 3):
        classes.append((p[2 + 3 * k], p[3 + 3 * k], p[4 + 3 * k]))
    return p[0], classes


def objective(p):
    big_r, classes = unpack(p)
    if not 2.0 < big_r < 9.0:
        return 9.0
    for alpha, tip, bend in classes:
        if not (0 <= alpha <= np.pi / 3 + 1e-9) or not (0.6 < tip < big_r - 1e-3) or abs(bend) > 1.5:
            return 9.0
    try:
        return bc.ratio(big_r, classes)
    except (np.linalg.LinAlgError, FloatingPointError, ValueError):
        return 9.0


def main():
    extra = int(sys.argv[1])
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    maxiter = int(sys.argv[3]) if len(sys.argv) > 3 else 300
    start = [float(v) for v in sys.argv[4:]]
    bounds = [(2.6, 7.5), (1.6, 2.4)] + [(0.03, np.pi / 3 - 0.03), (1.6, 7.0), (-0.8, 0.8)] * extra
    rng = np.random.default_rng(seed)
    pop = 20 * len(bounds)
    init = np.array([[rng.uniform(lo, hi) for lo, hi in bounds] for _ in range(pop)])
    if start:  # seed part of the population around a known good design
        for i in range(pop // 3):
            init[i] = np.clip(np.array(start) + rng.normal(0, 0.03, len(start)) * (i > 0), [b[0] for b in bounds], [b[1] for b in bounds])
    best = [9.0]

    def report(x, convergence=None):
        f = objective(x)
        if f < best[0] - 1e-9:
            best[0] = f
            print(f"  ratio {f:.9f}  params {np.round(x, 5).tolist()}", flush=True)

    print(f"differential evolution with curved slits, {extra} extra classes, seed {seed}, population {pop}", flush=True)
    res = differential_evolution(objective, bounds, seed=seed, maxiter=maxiter, tol=1e-13, mutation=(0.4, 1.0),
                                 recombination=0.8, polish=False, callback=report, init=init, updating="deferred",
                                 workers=int(os.environ.get("DE_WORKERS", "1")))  # one shape per core
    pol = minimize(objective, res.x, method="Nelder-Mead", options={"xatol": 1e-9, "fatol": 1e-13, "maxiter": 3000, "adaptive": True})
    x = pol.x if pol.fun < res.fun else res.x
    big_r, classes = unpack(x)
    print(f"FINAL: ratio {min(pol.fun, res.fun):.10f}")
    print(f"  params {[float(v) for v in x]}")
    slits = bc.slit_objects(big_r, classes)
    print(f"  inradius {bc.inradius(big_r, classes):.10f}; finer {bc.inradius(big_r, classes, grid=(400, 201), per_slit=600):.10f}; "
          f"conformal radius n=40 {conformal_radius(big_r, slits, 40):.12f}, n=96 {conformal_radius(big_r, slits, 96):.12f}")


if __name__ == "__main__":
    main()
