"""Convergence study of the SC parameter problem for the claimed domain.
Crowding is strong (gaps between prevertices down to 2e-3 degrees), so the
quadrature panels must be graded geometrically towards the singular
prevertices; we increase the number of graded panels until the answer stalls."""
import time
import numpy as np
from sc_slitmap import conformal_radius
from geom import R, claimed_slits

sl = sorted(claimed_slits())
th = np.array([a for a, t in sl])
tp = np.array([t for a, t in sl])
x0 = None
prev = None
CLAIM = 1.523513988544
# NOTE: npan >= 24 fails in double precision (innermost panel collapses onto the
# prevertex: 0/0).  Results obtained: npan=6 -> 1.523496857, npan=12 ->
# 1.523513988945, npan=18 -> 1.523513988919; (n=128,npan=8) -> 1.523513988947.
for n, npan in ((48, 6), (48, 12), (48, 18)):
    t0 = time.time()
    c, rmax, x0, (al, be, ga) = conformal_radius(R, th, tp, n=n, npan=npan, x0=x0, steps=8)
    print("n=%3d npan=%2d crad=%.15f  res=%.1e  change=%s  minus claim=%+.3e  (%.0fs)" %
          (n, npan, c, rmax, "%.2e" % (c - prev) if prev else "-", c - CLAIM, time.time() - t0),
          flush=True)
    prev = c
print("ratio 1/crad = %.13f" % (1 / c))
gaps = np.concatenate([be[1:], [be[0] + 2 * np.pi]]) - ga
print("min/max outer-arc prevertex gap (deg):", np.degrees(gaps.min()), np.degrees(gaps.max()))
print("total harmonic measure of outer circle from 0:", gaps.sum() / (2 * np.pi))
np.save("prevertices.npy", np.array([be, al, ga]))
