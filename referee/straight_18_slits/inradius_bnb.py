"""Branch-and-bound certificate (floating point, not interval arithmetic) that
the inradius is < 1 + EPS.  The distance-to-boundary function d is 1-Lipschitz,
so on a square with centre c and half-diagonal rho: sup d <= d(c) + rho.
Squares with d(c)+rho <= 1+EPS are discarded, the rest are subdivided.
No symmetry is used; the whole bounding square of the disc is covered."""
import numpy as np
from geom import R, claimed_slits, dist_to_boundary

slits = claimed_slits()
EPS = 1e-7

h = 0.05
xs = np.arange(-R - h, R + 2 * h, 2 * h)
X, Y = np.meshgrid(xs, xs)
C = (X + 1j * Y).ravel()
half = h
level = 0
best = -1.0
while C.size:
    d, _ = dist_to_boundary(C, slits, R)
    best = max(best, d.max())
    rho = half * np.sqrt(2)
    keep = d + rho > 1 + EPS
    n_keep = int(keep.sum())
    print("level %2d half-side %.3e cells %7d  surviving %6d  max d(centre)=%.12f" %
          (level, half, C.size, n_keep, d.max()))
    if n_keep == 0:
        break
    C = C[keep]
    if half < 1e-9:
        print("surviving cells cluster at:", np.unique(np.round(C, 6)))
        break
    half /= 2
    C = (C[:, None] + half * np.array([1 + 1j, 1 - 1j, -1 + 1j, -1 - 1j])[None, :]).ravel()
    level += 1
print("largest d found at any centre: %.12f" % best)
print("CERTIFIED (float arithmetic): inradius <= 1 + %.0e" % EPS if C.size == 0 or n_keep == 0
      else "not fully certified at this EPS")
