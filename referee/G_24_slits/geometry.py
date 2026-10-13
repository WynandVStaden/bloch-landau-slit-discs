"""Step 1: geometry.  In w=log z each slit is a graph theta(s) over s=log(|z|/R) in [-a,0].
Two slits meet iff at a common s their thetas agree mod 2pi  -> 1D check, done with exact
polynomial root finding on the difference (degree 4 in s, i.e. quadratic in s^2 plus shifts 2pi m)."""
import numpy as np
from scipy.spatial import cKDTree
from scipy.optimize import minimize
from dom import R_G as R, slits_G, zslit, dzslit, theta

S = slits_G()
print("R =", R)
print("idx cls  start_deg     tip_mod          tip_arg_deg    a=log(R/tip) arclen    max|angle to radial| deg")
for i, (c, al, tip, b1, b2) in enumerate(S):
    a = np.log(R / tip)
    zt = zslit(-1.0, R, al, tip, b1, b2)
    x = np.linspace(-1, 0, 400001)
    d = dzslit(x, R, al, tip, b1, b2)
    L = np.trapezoid(np.abs(d), x)
    s = a * x
    ang = np.degrees(np.arctan(np.abs(2 * b1 * s + 4 * b2 * s ** 3))).max()
    targ = np.degrees(al + b1 * a * a + b2 * a ** 4)
    print(f"{i:2d} {c} {np.degrees(al)%360:11.6f} {abs(zt):.13f} {targ%360:12.6f} (num {np.degrees(np.angle(zt))%360:10.6f}) "
          f"{a:.6f} {L:.6f} {ang:7.3f}")

# inside disc: |z| = R exp(a x) < R for x<0, = R only at x=0: by construction (a>0)
print("\nall a>0:", all(np.log(R / t[2]) > 0 for t in S), "-> |z| strictly increasing, in disc for x<0, on circle only at x=0")
print("total angular excursion per class (deg):")
for c in "CDE":
    cls, al, tip, b1, b2 = [t for t in S if t[0] == c][0]
    a = np.log(R / tip); s = np.linspace(-a, 0, 100001)
    th = theta(s, 0, b1, b2)
    print(f"  {c}: theta-alpha ranges [{np.degrees(th.min()):.4f}, {np.degrees(th.max()):.4f}] deg; d theta/ds sign changes: "
          f"{np.sum(np.diff(np.sign(2*b1*s[:-1]+4*b2*s[:-1]**3))!=0)}")

# exact intersection test
print("\nIntersection test (roots of theta_i(s)-theta_j(s) = 2 pi m on common s-range):")
nint = 0
gmin = (9, None)
for i in range(len(S)):
    for j in range(i + 1, len(S)):
        ci, ali, ti, b1i, b2i = S[i]
        cj, alj, tj, b1j, b2j = S[j]
        smin = -min(np.log(R / ti), np.log(R / tj))
        for m in range(-2, 3):
            # (b2i-b2j) u^2 + (b1i-b1j) u + (ali-alj-2 pi m) = 0, u=s^2 in [0, smin^2]
            co = [b2i - b2j, b1i - b1j, ali - alj - 2 * np.pi * m]
            if abs(co[0]) < 1e-300 and abs(co[1]) < 1e-300:
                if abs(co[2]) < 1e-12:
                    print("  COINCIDENT", i, j); nint += 1
                continue
            for u in np.roots(co):
                if abs(u.imag) < 1e-12 and -1e-12 <= u.real <= smin ** 2 + 1e-12:
                    print(f"  INTERSECTION slits {i}{ci}-{j}{cj} at s={-np.sqrt(max(u.real,0)):.6f}"); nint += 1
        s = np.linspace(smin, 0, 200001)
        gap = theta(s, ali, b1i, b2i) - theta(s, alj, b1j, b2j)
        gap = np.abs((gap + np.pi) % (2 * np.pi) - np.pi)
        k = gap.argmin()
        if gap[k] < gmin[0]:
            gmin = (gap[k], (i, j, R * np.exp(s[k])))
print("number of intersections found:", nint)
print("smallest angular gap at common modulus: %.5f deg, slits %s at |z|=%.4f" % (np.degrees(gmin[0]), gmin[1][:2], gmin[1][2]))

# Euclidean distances
N = 40001
X = np.linspace(-1, 0, N)
pts = [zslit(X, R, *t[1:]) for t in S]
trees = [cKDTree(np.c_[p.real, p.imag]) for p in pts]
rows = []
for i in range(len(S)):
    for j in range(i + 1, len(S)):
        d, idx = trees[j].query(np.c_[pts[i].real, pts[i].imag])
        k = d.argmin()
        # polish with exact curves
        f = lambda v: abs(zslit(np.clip(v[0], -1, 0), R, *S[i][1:]) - zslit(np.clip(v[1], -1, 0), R, *S[j][1:]))
        r = minimize(f, [X[k], X[idx[k]]], method="Nelder-Mead", options=dict(xatol=1e-13, fatol=1e-15))
        xi, xj = np.clip(r.x, -1, 0)
        rows.append((r.fun, i, j, xi, xj))
rows.sort()
print("\nclosest pairs of distinct slits (dist, pair, param on each [x=-1 tip, 0 circle], |z| of the two points):")
for dd, i, j, xi, xj in rows[:30]:
    zi = zslit(xi, R, *S[i][1:]); zj = zslit(xj, R, *S[j][1:])
    print(f"  {dd:.9f}  {i:2d}{S[i][0]}-{j:2d}{S[j][0]}  x=({xi:+.6f},{xj:+.6f})  |z|=({abs(zi):.5f},{abs(zj):.5f}) "
          f"arg=({np.degrees(np.angle(zi))%360:.3f},{np.degrees(np.angle(zj))%360:.3f})")
print("\nper class-pair minimum distance:")
best = {}
for dd, i, j, xi, xj in rows:
    key = "".join(sorted(S[i][0] + S[j][0]))
    if key not in best: best[key] = (dd, i, j, xi, xj)
for k, v in sorted(best.items()):
    print(f"  {k}: {v[0]:.9f} (slits {v[1]},{v[2]}; x={v[3]:+.5f},{v[4]:+.5f})")
# start-point separation on the circle
al = np.sort(np.array([t[1] % (2 * np.pi) for t in S]))
gaps = np.diff(np.append(al, al[0] + 2 * np.pi))
print("\nstart angles on circle sorted (deg):", np.round(np.degrees(al), 4))
print("min gap between start points: %.5f deg = chord %.6f" % (np.degrees(gaps.min()), 2 * R * np.sin(gaps.min() / 2)))
print("min distance from 0 to slits:", min(np.abs(p).min() for p in pts))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(9, 9))
th = np.linspace(0, 2 * np.pi, 2000)
ax.plot(R * np.cos(th), R * np.sin(th), "k-", lw=1)
col = dict(A="C0", B="C1", C="C2", D="C3", E="C4")
for t, p in zip(S, pts):
    ax.plot(p.real, p.imag, color=col[t[0]], lw=1.5)
import json, os
if os.path.exists("maxima.json"):
    for m in json.load(open("maxima.json")):
        ax.add_patch(plt.Circle((m["x"], m["y"]), m["r"], fill=False, color="gray", lw=0.5))
ax.set_aspect("equal")
fig.savefig("domain.png", dpi=120)
