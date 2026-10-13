"""Search for slit domains with a small ratio inradius / conformal radius.

A design is a disc of radius R with radial slits, invariant under rotation by
120 degrees and reflection in the real axis. It is described by R and a list
of classes (angle in [0, 60 degrees], tip radius); each class stands for 3 or
6 slits. The first class is fixed as (0, 1): slits [1, R] at 0/120/240 degrees.

Usage: python bloch_design.py
"""
import numpy as np
from scipy.optimize import minimize
from scipy.spatial import Voronoi

from bloch_solver import Slit, conformal_radius, symmetric_slits


def segments(big_r, classes):
    """End points (tip, base) of every slit, as complex numbers."""
    out = []
    for alpha, tip in classes:
        angles = set()
        for r in range(3):
            for sgn in (1, -1):
                angles.add(round((2 * np.pi * r / 3 + sgn * alpha) % (2 * np.pi), 12))
        for a in sorted(angles):
            out.append((tip * np.exp(1j * a), big_r * np.exp(1j * a)))
    return out


def clearance(p, big_r, segs):
    """Distance from the point p to the boundary (circle and slits)."""
    z = complex(p[0], p[1])
    u = np.fromiter((s[0] for s in segs), complex, len(segs))
    w = np.fromiter((s[1] - s[0] for s in segs), complex, len(segs))
    lam = np.clip(((z - u) * np.conj(w)).real / np.abs(w) ** 2, 0.0, 1.0)
    return min(big_r - abs(z), float(np.abs(z - (u + lam * w)).min()))


def feature_distances(z, big_r, u, w):
    """Distances from points z (array) to every slit and to the circle, with the nearest boundary points."""
    z = np.atleast_1d(z)[:, None]
    lam = np.clip(((z - u) * np.conj(w)).real / np.abs(w) ** 2, 0.0, 1.0)
    near = u + lam * w
    d_seg = np.abs(z - near)
    absz = np.abs(z)
    safe = np.where(absz == 0, 1.0, absz)
    near_c = np.where(absz == 0, big_r, big_r * z / safe)
    return np.hstack([d_seg, big_r - absz]), np.hstack([near, near_c])


def inradius(big_r, classes, grid=(240, 121)):
    """Radius of the largest disc inside the slit disc.

    A polar grid over the sector 0..60 degrees (enough by symmetry) locates the
    local maxima of the distance to the boundary; each is then refined by
    Newton's method on the three equations 'distance to feature k = r'.
    """
    segs = segments(big_r, classes)
    u = np.array([s[0] for s in segs])
    w = np.array([s[1] - s[0] for s in segs])
    rr, aa = np.meshgrid(np.linspace(0, big_r, grid[0]), np.linspace(0, np.pi / 3, grid[1]), indexing="ij")
    z = (rr * np.exp(1j * aa)).ravel()
    d, _ = feature_distances(z, big_r, u, w)
    clear = d.min(axis=1).reshape(grid)
    # local maxima of the grid
    pad = np.pad(clear, 1, constant_values=-1)
    is_max = np.ones_like(clear, bool)
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di or dj:
                is_max &= clear >= pad[1 + di:1 + di + grid[0], 1 + dj:1 + dj + grid[1]]
    idx = np.argwhere(is_max)
    idx = idx[np.argsort(-clear[is_max])]
    # keep well-separated maxima only, so that a ridge of near-ties cannot crowd out a distinct maximum
    chosen, centres = [], []
    for i, j in idx:
        c = rr[i, j] * np.exp(1j * aa[i, j])
        if all(abs(c - o) > 0.12 for o in centres):
            chosen.append((i, j))
            centres.append(c)
        if len(chosen) >= 40:
            break
    best = float(clear.max())
    for i, j in chosen:
        p = complex(rr[i, j] * np.exp(1j * aa[i, j]))
        r = clear[i, j]
        for _ in range(30):
            dist, near = feature_distances(np.array([p]), big_r, u, w)
            dist, near = dist[0], near[0]
            act = np.argsort(dist)[:3]
            grad = np.empty((3, 2))
            for row, k in enumerate(act):
                if k == len(segs):  # the circle: distance R - |p|
                    g = -p / abs(p) if abs(p) > 0 else 0j
                else:
                    g = (p - near[k]) / dist[k] if dist[k] > 0 else 0j
                grad[row] = (g.real, g.imag)
            jac = np.hstack([grad, -np.ones((3, 1))])
            try:
                step = np.linalg.solve(jac, -(dist[act] - r))
            except np.linalg.LinAlgError:
                break
            p += complex(step[0], step[1])
            r += step[2]
            if abs(step[0]) + abs(step[1]) + abs(step[2]) < 1e-14:
                break
        true = feature_distances(np.array([p]), big_r, u, w)[0][0].min() if abs(p) < big_r else 0.0
        best = max(best, float(true))
    return best


def slit_list(big_r, classes):
    out = []
    for alpha, tip in classes:
        out += symmetric_slits(3, alpha, np.log(big_r / tip))
    return out


def ratio(big_r, classes, n=40):
    return inradius(big_r, classes) / conformal_radius(big_r, slit_list(big_r, classes), n)


def unpack(p, n_classes):
    big_r = p[0]
    classes = [(0.0, 1.0), (np.pi / 3, p[1])]
    for k in range(n_classes - 2):
        classes.append((p[2 + 2 * k], p[3 + 2 * k]))
    return big_r, classes


def objective(p, n_classes):
    big_r, classes = unpack(p, n_classes)
    if big_r < 1.5 or big_r > 12:
        return 9.0
    for alpha, tip in classes:
        if not (0 <= alpha <= np.pi / 3 + 1e-9) or not (0.5 < tip < big_r - 1e-3):
            return 9.0
    return ratio(big_r, classes)


def optimise(start, n_classes, rounds=3):
    p = np.array(start, float)
    best = objective(p, n_classes)
    for _ in range(rounds):
        res = minimize(objective, p, args=(n_classes,), method="Nelder-Mead",
                       options={"xatol": 1e-7, "fatol": 1e-11, "maxiter": 1500, "adaptive": True})
        if res.fun < best:
            best, p = res.fun, res.x
    return best, p


if __name__ == "__main__":
    print("known: Robinson 0.658, Goodman 0.65647, Beller-Hummel 0.6564155, Carroll-Ortega-Cerda 0.6563937\n")
    b2, p2 = optimise([2.9562, 2.0], 2)
    print(f"2 classes (6 slits): ratio {b2:.8f}  R={p2[0]:.5f} t2={p2[1]:.5f}", flush=True)
    rng = np.random.default_rng(0)
    best3 = (9, None)
    for trial in range(12):
        start = [rng.uniform(3.3, 4.6), rng.uniform(1.8, 2.2), rng.uniform(0.3, 0.75), rng.uniform(2.6, 3.6)]
        b3, p3 = optimise(start, 3, rounds=2)
        if b3 < best3[0]:
            best3 = (b3, p3)
            print(f"3 classes (12 slits): ratio {b3:.9f}  R={p3[0]:.6f} t2={p3[1]:.6f} angle3={np.degrees(p3[2]):.4f}deg t3={p3[3]:.6f}", flush=True)
    best4 = (9, None)
    for trial in range(16):
        b, p = best3[1]
        start = list(best3[1]) + [rng.uniform(0.1, 0.9), rng.uniform(3.2, 5.5)]
        start[0] = rng.uniform(best3[1][0], best3[1][0] + 1.5)
        b4, p4 = optimise(start, 4, rounds=2)
        if b4 < best4[0]:
            best4 = (b4, p4)
            print(f"4 classes: ratio {b4:.9f}  params {np.round(p4, 6)}", flush=True)
