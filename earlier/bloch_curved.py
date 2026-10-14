"""Slit designs with curved slits.

A design is a disc of radius R with slits, invariant under rotation by 120
degrees and reflection in the real axis. Each class is (angle, tip radius,
bend): the slit leaves the circle at the given polar angle, at a right angle,
and in the coordinate w = log z it is the curve

    w(x) = log R + a x + i (angle + bend * (a x)^2),   -1 <= x <= 0,   a = log(R / tip radius).

bend = 0 is a straight radial slit. The conformal radius comes from
bloch_solver (which already handles these curves exactly); this module
supplies the inradius.
"""
import os as _os
import sys as _sys

_sys.path.append(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))  # shared modules in the parent folder

import numpy as np
from scipy.spatial import cKDTree

from bloch_solver import Slit, conformal_radius


def expand(classes):
    """All slits of a design as arrays (angle, tip radius, bend), using the symmetry."""
    out = {}
    for alpha, tip, bend in classes:
        for r in range(3):
            for sgn in (1, -1):
                ang = (2 * np.pi * r / 3 + sgn * alpha) % (2 * np.pi)
                key = round(ang, 10)
                if key not in out:
                    out[key] = (ang, tip, sgn * bend)
    arr = np.array(list(out.values()))
    return arr[:, 0], arr[:, 1], arr[:, 2]


def curve(x, big_r, ang, a, bend):
    s = a * x
    return big_r * np.exp(s + 1j * (ang + bend * s * s))


def nearest_on_slits(p, big_r, ang, a, bend, x0):
    """Nearest point of every slit to the point p, by Newton's method on the curve parameter."""
    x = x0.copy()
    for _ in range(12):
        s = a * x
        z = big_r * np.exp(s + 1j * (ang + bend * s * s))
        k = a * (1 + 2j * bend * s)
        dz = z * k
        d2z = z * (k * k + 2j * bend * a * a)
        g = ((z - p) * np.conj(dz)).real
        h = (np.abs(dz) ** 2 + ((z - p) * np.conj(d2z)).real)
        step = np.where(h > 1e-14, g / np.where(h > 1e-14, h, 1.0), 0.0)
        x = np.clip(x - step, -1.0, 0.0)
    z = curve(x, big_r, ang, a, bend)
    return z, x


def inradius(big_r, classes, grid=(150, 76), per_slit=160):
    """Radius of the largest disc inside the design."""
    ang, tip, bend = expand(classes)
    a = np.log(big_r / tip)
    xs = -np.linspace(0, 1, per_slit)
    slit_pts = curve(xs[None, :], big_r, ang[:, None], a[:, None], bend[:, None])  # (slits, per_slit)
    circle = big_r * np.exp(1j * np.linspace(0, 2 * np.pi, 720, endpoint=False))
    boundary = np.concatenate([slit_pts.ravel(), circle])
    tree = cKDTree(np.column_stack([boundary.real, boundary.imag]))
    rr, aa = np.meshgrid(np.linspace(0, big_r, grid[0]), np.linspace(0, np.pi / 3, grid[1]), indexing="ij")
    z = rr * np.exp(1j * aa)
    clear = np.minimum(tree.query(np.column_stack([z.real.ravel(), z.imag.ravel()]))[0].reshape(grid), big_r - rr)
    # the distance to the boundary is 1-Lipschitz, so every maximum within the grid spacing of the best one is kept
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
    best = 0.0
    n_slits = len(ang)
    for p in centres:
        # start each slit's parameter at its sampled point nearest to p
        x = xs[np.abs(slit_pts - p).argmin(axis=1)]
        r = None
        for _ in range(40):
            near, x = nearest_on_slits(p, big_r, ang, a, bend, x)
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
            if norm > 0.05:  # damped far from the solution
                step *= 0.05 / norm
            p = p + complex(step[0], step[1])
            r += step[2]
            if abs(step).sum() < 1e-14:
                break
        if abs(p) < big_r:
            near, x = nearest_on_slits(p, big_r, ang, a, bend, xs[np.abs(slit_pts - p).argmin(axis=1)])
            best = max(best, float(min(np.abs(p - near).min(), big_r - abs(p))))
    return max(best, 0.0)


def slit_objects(big_r, classes):
    ang, tip, bend = expand(classes)
    return [Slit(float(al), float(np.log(big_r / t)), float(b)) for al, t, b in zip(ang, tip, bend)]


def ratio(big_r, classes, n=40):
    return inradius(big_r, classes) / conformal_radius(big_r, slit_objects(big_r, classes), n)


if __name__ == "__main__":
    import time

    import bloch_design as d

    print("straight designs: curved-capable routine vs the straight-slit routine")
    tests = [
        (4.11257, [(0.0, 1.0), (np.pi / 3, 2.0), (0.39304, 2.8961)]),
        (2.95617, [(0.0, 1.0), (np.pi / 3, 2.0)]),
        (4.98464, [(0.0, 1.0), (np.pi / 3, 2.00002), (0.60612, 3.92526), (0.39572, 2.88656)]),
    ]
    for big_r, cl in tests:
        t = time.time()
        r1 = inradius(big_r, [(al, tp, 0.0) for al, tp in cl])
        dt = time.time() - t
        r2 = d.inradius(big_r, cl, grid=(480, 241))
        c = conformal_radius(big_r, slit_objects(big_r, [(al, tp, 0.0) for al, tp in cl]), 40)
        print(f"  inradius {r1:.10f} vs {r2:.10f}  difference {r1 - r2:+.1e}  time {dt:.3f}s  ratio {r1 / c:.9f}")
    print("curved design: routine vs brute force over random points (lower estimate)")
    rng = np.random.default_rng(0)
    for big_r, cl in [(4.1, [(0.0, 1.0, 0.0), (np.pi / 3, 2.0, 0.0), (0.4, 2.9, 0.15)]), (4.4, [(0.0, 1.0, 0.0), (np.pi / 3, 2.0, 0.0), (0.45, 2.8, -0.2), (0.8, 3.6, 0.1)])]:
        ang, tip, bend = expand(cl)
        a = np.log(big_r / tip)
        pts = curve(-np.linspace(0, 1, 4000)[None, :], big_r, ang[:, None], a[:, None], bend[:, None]).ravel()
        tree = cKDTree(np.column_stack([pts.real, pts.imag]))
        brute = 0.0
        for _ in range(10):
            q = np.sqrt(rng.random(400000)) * big_r * np.exp(1j * np.pi / 3 * rng.random(400000))
            brute = max(brute, float(np.minimum(tree.query(np.column_stack([q.real, q.imag]))[0], big_r - np.abs(q)).max()))
        print(f"  routine {inradius(big_r, cl):.8f}  brute force {brute:.6f}")
