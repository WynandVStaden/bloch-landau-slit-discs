"""Slit designs with more flexible curved slits, and a differential-evolution search over them.

As in bloch_curved.py, but a slit is now, in the coordinate w = log z,

    w(x) = log R + s + i (alpha + b1 s^2 + b2 s^4),   s = a x,  -1 <= x <= 0,  a = log(R / tip).

Only even powers of s appear, so the slit together with its reflection in the
circle is still one analytic arc and bloch_solver.conformal_radius applies
unchanged.

Usage: python bloch_curved2.py n_extra_classes [seed] [maxiter]
"""
import os
import sys

import numpy as np
from scipy.optimize import differential_evolution, minimize
from scipy.spatial import cKDTree

from bloch_solver import conformal_radius


class Slit2:
    def __init__(self, alpha, depth, b1=0.0, b2=0.0):
        self.alpha, self.a, self.b1, self.b2 = alpha, depth, b1, b2

    def w(self, x, big_l):
        s = self.a * x
        return big_l + s + 1j * (self.alpha + self.b1 * s**2 + self.b2 * s**4)

    def dw(self, x):
        s = self.a * x
        return self.a * (1 + 1j * (2 * self.b1 * s + 4 * self.b2 * s**3))


def expand(classes):
    """All slits as arrays (angle, tip, b1, b2), using rotation by 120 degrees and reflection."""
    out = {}
    for alpha, tip, b1, b2 in classes:
        for r in range(3):
            for sgn in (1, -1):
                ang = (2 * np.pi * r / 3 + sgn * alpha) % (2 * np.pi)
                key = round(ang, 10)
                if key not in out:
                    out[key] = (ang, tip, sgn * b1, sgn * b2)
    arr = np.array(sorted(out.values()))
    return arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]


def curve(x, big_r, ang, a, b1, b2):
    s = a * x
    return big_r * np.exp(s + 1j * (ang + b1 * s**2 + b2 * s**4))


def slits_cross(big_r, ang, tip, b1, b2, margin=0.03):
    """True if two slits cross or come closer than margin. |z| is monotone along a slit, so compare angles radius by radius."""
    a = np.log(big_r / tip)
    for rad in np.linspace(tip.min(), big_r, 80):
        s = np.log(rad / big_r)
        present = tip <= rad + 1e-12
        # slits are listed in order of their angle at the circle; that angular order must persist at every radius
        th = (ang + b1 * s**2 + b2 * s**4)[present]
        if len(th) > 1:
            gaps = np.diff(np.concatenate([th, [th[0] + 2 * np.pi]]))
            if gaps.min() * rad < margin:
                return True
        # a slit must not wander by more than a full sector
        if np.abs((b1 * s**2 + b2 * s**4)[present]).max() > 1.0:
            return True
    return False


def nearest_on_slits(p, big_r, ang, a, b1, b2, x0):
    x = x0.copy()
    for _ in range(12):
        s = a * x
        z = big_r * np.exp(s + 1j * (ang + b1 * s**2 + b2 * s**4))
        k = a * (1 + 1j * (2 * b1 * s + 4 * b2 * s**3))
        dz = z * k
        d2z = z * (k * k + 1j * a * a * (2 * b1 + 12 * b2 * s**2))
        g = ((z - p) * np.conj(dz)).real
        h = np.abs(dz) ** 2 + ((z - p) * np.conj(d2z)).real
        step = np.where(h > 1e-14, g / np.where(h > 1e-14, h, 1.0), 0.0)
        x = np.clip(x - step, -1.0, 0.0)
    return curve(x, big_r, ang, a, b1, b2), x


def inradius(big_r, classes, grid=(150, 76), per_slit=160):
    ang, tip, b1, b2 = expand(classes)
    a = np.log(big_r / tip)
    xs = -np.linspace(0, 1, per_slit)
    slit_pts = curve(xs[None, :], big_r, ang[:, None], a[:, None], b1[:, None], b2[:, None])
    circle = big_r * np.exp(1j * np.linspace(0, 2 * np.pi, 720, endpoint=False))
    boundary = np.concatenate([slit_pts.ravel(), circle])
    tree = cKDTree(np.column_stack([boundary.real, boundary.imag]))
    rr, aa = np.meshgrid(np.linspace(0, big_r, grid[0]), np.linspace(0, np.pi / 3, grid[1]), indexing="ij")
    z = rr * np.exp(1j * aa)
    clear = np.minimum(tree.query(np.column_stack([z.real.ravel(), z.imag.ravel()]))[0].reshape(grid), big_r - rr)
    h = max(big_r / (grid[0] - 1), big_r * (np.pi / 3) / (grid[1] - 1)) + 0.02
    pad = np.pad(clear, 1, constant_values=-1)
    is_max = np.ones_like(clear, bool)
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di or dj:
                is_max &= clear >= pad[1 + di:1 + di + grid[0], 1 + dj:1 + dj + grid[1]]
    is_max &= clear >= clear.max() - 2 * h
    idx = np.argwhere(is_max)
    idx = idx[np.argsort(-clear[is_max])]
    centres = []
    for i, j in idx:
        c = z[i, j]
        if all(abs(c - o) > 0.08 for o in centres):
            centres.append(c)
        if len(centres) >= 60:
            break
    best, n_slits = 0.0, len(ang)
    for p in centres:
        x = xs[np.abs(slit_pts - p).argmin(axis=1)]
        r = None
        for _ in range(40):
            near, x = nearest_on_slits(p, big_r, ang, a, b1, b2, x)
            dist = np.append(np.abs(p - near), big_r - abs(p))
            if r is None:
                r = dist.min()
            act = np.argsort(dist)[:3]
            grad = np.empty((3, 2))
            for row, k in enumerate(act):
                g = (-p / abs(p) if abs(p) > 0 else 0j) if k == n_slits else (p - near[k]) / dist[k]
                grad[row] = (g.real, g.imag)
            try:
                step = np.linalg.solve(np.hstack([grad, -np.ones((3, 1))]), -(dist[act] - r))
            except np.linalg.LinAlgError:
                break
            norm = np.hypot(step[0], step[1])
            if norm > 0.05:
                step *= 0.05 / norm
            p = p + complex(step[0], step[1])
            r += step[2]
            if abs(step).sum() < 1e-14:
                break
        if abs(p) < big_r:
            near, x = nearest_on_slits(p, big_r, ang, a, b1, b2, xs[np.abs(slit_pts - p).argmin(axis=1)])
            best = max(best, float(min(np.abs(p - near).min(), big_r - abs(p))))
    return best


def slit_objects(big_r, classes):
    ang, tip, b1, b2 = expand(classes)
    return [Slit2(float(al), float(np.log(big_r / t)), float(u), float(v)) for al, t, u, v in zip(ang, tip, b1, b2)]


def unpack(p):
    classes = [(0.0, 1.0, 0.0, 0.0), (np.pi / 3, p[1], 0.0, 0.0)]
    for k in range((len(p) - 2) // 4):
        classes.append(tuple(p[2 + 4 * k:6 + 4 * k]))
    return p[0], classes


def objective(p):
    big_r, classes = unpack(p)
    if not 2.0 < big_r < 9.0:
        return 9.0
    for alpha, tip, b1, b2 in classes:
        if not (0 <= alpha <= np.pi / 3 + 1e-9) or not (0.6 < tip < big_r - 1e-3):
            return 9.0
    ang, tip, b1, b2 = expand(classes)
    if slits_cross(big_r, ang, tip, b1, b2):
        return 9.0
    try:
        return inradius(big_r, classes) / conformal_radius(big_r, slit_objects(big_r, classes), 40)
    except (np.linalg.LinAlgError, FloatingPointError, ValueError):
        return 9.0


BEST_SO_FAR = [4.747749358688023, 1.9999999999999036,
               0.12393306687411461, 2.8051388295325297, 0.6738078607105484, 0.0,
               0.7818321423345288, 3.710751140554841, -0.6978713930533904, 0.0]


def main():
    extra = int(sys.argv[1])
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    maxiter = int(sys.argv[3]) if len(sys.argv) > 3 else 500
    bounds = [(2.6, 8.0), (1.6, 2.4)] + [(0.03, np.pi / 3 - 0.03), (1.6, 7.5), (-1.5, 1.5), (-3.0, 3.0)] * extra
    rng = np.random.default_rng(seed)
    pop = 15 * len(bounds)
    lo, hi = np.array([b[0] for b in bounds]), np.array([b[1] for b in bounds])
    init = rng.uniform(lo, hi, (pop, len(bounds)))
    # a third of the population starts from the best known design; an extra class starts as a very short slit
    start = list(BEST_SO_FAR)
    for _ in range(extra - 2):
        start += [0.5, start[0] - 0.02, 0.0, 0.0]
    start = np.array(start[:len(bounds)])
    for i in range(pop // 3):
        cand = start + rng.normal(0, 0.02, len(bounds)) * (i > 0)
        if extra > 2 and i > 0:  # try the new class at assorted places
            cand[10] = rng.uniform(0.05, 1.0)
            cand[11] = rng.uniform(3.2, cand[0] - 0.05)
        init[i] = np.clip(cand, lo, hi)
    best = [9.0]

    def report(x, convergence=None):
        f = objective(x)
        if f < best[0] - 2e-9:
            best[0] = f
            print(f"  ratio {f:.9f}  params {np.round(x, 6).tolist()}", flush=True)

    print(f"flexible curved slits: {extra} extra classes, {len(bounds)} parameters, population {pop}, seed {seed}", flush=True)
    print(f"  starting design scores {objective(start):.9f}", flush=True)
    res = differential_evolution(objective, bounds, seed=seed, maxiter=maxiter, tol=1e-13, mutation=(0.4, 1.0),
                                 recombination=0.8, polish=False, callback=report, init=init, updating="deferred",
                                 workers=int(os.environ.get("DE_WORKERS", "1")))
    pol = minimize(objective, res.x, method="Nelder-Mead", options={"xatol": 1e-9, "fatol": 1e-13, "maxiter": 4000, "adaptive": True})
    x = pol.x if pol.fun < res.fun else res.x
    big_r, classes = unpack(x)
    slits = slit_objects(big_r, classes)
    print(f"FINAL: ratio {min(pol.fun, res.fun):.10f}", flush=True)
    print(f"  params {[float(v) for v in x]}")
    print(f"  inradius {inradius(big_r, classes):.10f}; finer {inradius(big_r, classes, grid=(400, 201), per_slit=600):.10f}; "
          f"conformal radius n=40 {conformal_radius(big_r, slits, 40):.12f}, n=96 {conformal_radius(big_r, slits, 96):.12f}")


if __name__ == "__main__":
    main()
