"""Differential evolution over slit designs (see bloch_design.py for the model).

Usage: python bloch_evolve.py n_classes [seed] [maxiter]
"""
import sys

import numpy as np
from scipy.optimize import differential_evolution, minimize

import bloch_design as d


def main():
    n_classes = int(sys.argv[1])
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    maxiter = int(sys.argv[3]) if len(sys.argv) > 3 else 400
    # R, tip radius of the 60-degree slits, then (angle, tip radius) per extra class
    bounds = [(2.6, 7.0), (1.5, 2.6)] + [(0.02, np.pi / 3 - 0.02), (1.5, 6.5)] * (n_classes - 2)
    best = [9.0]

    def report(x, convergence=None):
        f = d.objective(x, n_classes)
        if f < best[0] - 1e-9:
            best[0] = f
            print(f"  ratio {f:.9f}  params {np.round(x, 5).tolist()}", flush=True)

    print(f"differential evolution, {n_classes} classes, seed {seed}", flush=True)
    res = differential_evolution(d.objective, bounds, args=(n_classes,), seed=seed, popsize=24, maxiter=maxiter,
                                 tol=1e-12, mutation=(0.4, 1.0), recombination=0.8, polish=False, callback=report,
                                 init="sobol", updating="deferred", workers=1)
    pol = minimize(d.objective, res.x, args=(n_classes,), method="Nelder-Mead",
                   options={"xatol": 1e-9, "fatol": 1e-13, "maxiter": 4000, "adaptive": True})
    x = pol.x if pol.fun < res.fun else res.x
    big_r, classes = d.unpack(x, n_classes)
    print(f"FINAL {n_classes} classes: ratio {min(pol.fun, res.fun):.10f}")
    print(f"  R = {big_r:.8f}; classes (angle in degrees, tip radius): {[(round(float(np.degrees(a)), 5), round(float(t), 7)) for a, t in classes]}")
    print(f"  inradius {d.inradius(big_r, classes):.10f}, finer grid {d.inradius(big_r, classes, grid=(600, 301)):.10f}; "
          f"conformal radius n=40: {d.conformal_radius(big_r, d.slit_list(big_r, classes), 40):.12f}, n=80: {d.conformal_radius(big_r, d.slit_list(big_r, classes), 80):.12f}")


if __name__ == "__main__":
    main()
