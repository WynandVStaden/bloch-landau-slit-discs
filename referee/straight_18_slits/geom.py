"""Geometry of the claimed domain (written from scratch by the referee)."""
import numpy as np

R = 5.99717
DEG = np.pi / 180.0


def claimed_slits():
    """Return list of (angle_rad, tip_radius) for the 18 slits."""
    s = []
    for k in range(3):
        rot = 120.0 * k
        s.append((rot + 0.0, 1.0))
        s.append((rot + 60.0, 1.99998))
        for sg in (+1, -1):
            s.append((rot + sg * 41.47240, 3.87035))
            s.append((rot + sg * 19.61578, 2.84221))
    return [((a % 360.0) * DEG, t) for a, t in s]


def dist_to_boundary(P, slits, R):
    """Exact distance from points P (complex array) to the boundary
    = min(R-|p|, min_j dist(p, segment_j)).  Returns (d, idx) where idx = -1
    for the circle, else the slit index.  Points with |p|>=R get d<=0."""
    P = np.asarray(P, dtype=complex)
    d = R - np.abs(P)
    idx = -np.ones(P.shape, dtype=int)
    for j, (a, t) in enumerate(slits):
        u = np.exp(1j * a)
        q = P * np.conj(u)            # rotate slit onto positive real axis
        x = np.clip(q.real, t, R)     # nearest point on [t, R]
        dj = np.hypot(q.real - x, q.imag)
        m = dj < d
        d = np.where(m, dj, d)
        idx = np.where(m, j, idx)
    return d, idx
