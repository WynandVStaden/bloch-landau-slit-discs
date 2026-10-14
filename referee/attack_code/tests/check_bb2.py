"""Sampling attack on bound_abs: for the boxes that plan() uses, compare the certified (lo, hi) with |kernel| at random
points computed by an independent float implementation (direct formula, no divided difference)."""
import sys, os, json, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "certify"))
import numpy as np
import conformal_certificate as cc
from flint import arb, ctx
ctx.prec = cc.PREC
geo = cc.Geometry("G")
rng = np.random.default_rng(7)
Lf = math.log(float(geo.big_r))
def fw(j, x):
    s = geo.fa[j] * x
    return Lf + s + 1j * (geo.fth[j] + geo.fb1[j] * s**2 + geo.fb2[j] * s**4)
def kern(i, j, x, y):
    v = 2 * np.sinh((fw(i, x) - fw(j, y)) / 2)
    return v / (x - y) if i == j else v
def sample(box, n):
    cols = []
    for lo, hi in box:
        u = rng.random(n)
        u[: n // 4] = rng.integers(0, 2, n // 4)  # a quarter of the samples on faces/corners
        rng.shuffle(u)
        cols.append(lo + (hi - lo) * u)
    return cols[0] + 1j * cols[1], cols[2] + 1j * cols[3]
worst_lo, worst_hi, nb = 9, 9, 0
for i in geo.reps[:1]:
    others = [j for j in range(geo.n_slits) if j != i]
    dist = min([cc.float_distance(geo, i, j, (-1, 1)) for j in others] + [9.0])
    b = min(0.3, cc.SAFETY * dist / geo.speed[i])
    for step in (1, 2):
        for j in [i] + others:
            if step == 1:
                bb = b
                while True:
                    a_half, r = cc.ellipse(bb)
                    box = [(-a_half, a_half), (-bb, bb), (-1.0, 1.0), (0.0, 0.0)]
                    scale = [geo.speed[i], geo.speed[i], geo.speed[j], 0.0]
                    func = (lambda x, y: geo.self_kernel(i, x, y)) if j == i else (lambda x, y: geo.cross(i, j, x, y))
                    bd = cc.bound_abs(func, box, scale)
                    if bd is not None: break
                    bb /= 2
            else:
                bq = 0.3 if j == i else min(0.3, cc.SAFETY * cc.float_distance(geo, i, j, (-1, 0)) / geo.speed[j])
                while True:
                    a_half, r = cc.ellipse(bq)
                    box = [(-1.0, 0.0), (0.0, 0.0), (-a_half, a_half), (-bq, bq)]
                    func = (lambda x, y: geo.self_kernel(i, x, y)) if j == i else (lambda x, y: geo.cross(i, j, x, y))
                    bd = cc.bound_abs(func, box, [geo.speed[i], 0.0, geo.speed[j], geo.speed[j]])
                    if bd is not None: break
                    bq /= 2
            lo, hi = float(bd[0]), float(bd[1])
            assert bd[0].rad() == 0 and bd[1].rad() == 0, "abs_lower/abs_upper not exact"
            x, y = sample(box, 40000)
            v = np.abs(kern(i, j, x, y))
            v = v[np.isfinite(v) & (np.abs(x - y) > 1e-6)]
            mlo, mhi = v.min() / lo, hi / v.max()
            worst_lo, worst_hi, nb = min(worst_lo, mlo), min(worst_hi, mhi), nb + 1
            if i == 0 and j in (0, 1, 3, 6, 7, 12, 18): print(i, j, step, f"lo={lo:.6g} samplemin={v.min():.6g} hi={hi:.6g} samplemax={v.max():.6g} ratios {mlo:.9f} {mhi:.9f}")
            if mlo < 1 - 1e-9 or mhi < 1 - 1e-9:
                print("VIOLATION", i, j, step, lo, hi, v.min(), v.max())
    print(f"rep {i}: done; so far min(sample min / lo) = {worst_lo:.4f}, min(hi / sample max) = {worst_hi:.4f}, boxes {nb}", flush=True)
