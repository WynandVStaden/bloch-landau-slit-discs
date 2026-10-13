"""Step 1: geometry sanity."""
import numpy as np
from scipy.spatial import cKDTree
from dom import R_CLAIM as R, claimed_slits, zslit, dzslit

S = claimed_slits()
print("idx cls alpha_deg  tip_mod  tip_arg_deg  a=log(R/tip)  arclength  max|turn|deg")
for i, (c, al, tip, b) in enumerate(S):
    a = np.log(R / tip)
    zt = zslit(-1.0, R, al, tip, b)
    x = np.linspace(-1, 0, 200001)
    d = dzslit(x, R, al, tip, b)
    L = np.trapezoid(np.abs(d), x)
    # angle between tangent and radial direction: atan(2 bend u)
    ang = np.degrees(np.arctan(2 * b * a))
    print(f"{i:2d} {c} {np.degrees(al):10.5f} {abs(zt):.13f} {np.degrees(np.angle(zt)) % 360:11.6f} "
          f"{a:.6f} {L:.6f} {ang:8.3f}")

# modulus is strictly monotone along each slit (|z| = R exp(a x)) => no self intersection,
# inside disc for x<0, touches circle only at x=0.
# Two slits can only meet at a common modulus rho; compare arguments there.
print("\nPairwise check: min over common moduli of angular gap (mod 2pi), and min Euclidean distance")
N = 20001
pts = [zslit(np.linspace(-1, 0, N), R, al, tip, b) for (c, al, tip, b) in S]
trees = [cKDTree(np.c_[p.real, p.imag]) for p in pts]
gmin = (9, None)
dmin = (9, None)
for i in range(len(S)):
    for j in range(i + 1, len(S)):
        ci, ali, ti, bi = S[i]
        cj, alj, tj, bj = S[j]
        rho = np.linspace(max(ti, tj), R, 200001)
        u = np.log(rho / R)
        gap = (ali + bi * u * u) - (alj + bj * u * u)
        gap = np.abs((gap + np.pi) % (2 * np.pi) - np.pi)
        k = gap.argmin()
        d, _ = trees[j].query(np.c_[pts[i].real, pts[i].imag])
        dd = d.min()
        if gap[k] < gmin[0]:
            gmin = (gap[k], (i, j, rho[k]))
        if dd < dmin[0]:
            dmin = (dd, (i, j))
        if dd < 1.0:
            print(f"  slits {i:2d}{ci}-{j:2d}{cj}: min ang gap {np.degrees(gap[k]):8.4f} deg at rho={rho[k]:.4f}; "
                  f"min dist {dd:.6f}")
print("global min angular gap (deg):", np.degrees(gmin[0]), gmin[1])
print("global min distance between distinct slits:", dmin)
print("min distance from 0 to slits:", min(np.abs(p).min() for p in pts))

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 8))
    th = np.linspace(0, 2 * np.pi, 2000)
    ax.plot(R * np.cos(th), R * np.sin(th), "k-", lw=1)
    col = dict(A="C0", B="C1", C="C2", D="C3")
    for (c, al, tip, b), p in zip(S, pts):
        ax.plot(p.real, p.imag, color=col[c], lw=1.5)
    import json, os
    if os.path.exists("inradius_maxima.json"):
        for m in json.load(open("inradius_maxima.json")):
            ax.add_patch(plt.Circle((m["x"], m["y"]), m["r"], fill=False, color="gray", lw=0.5))
    ax.set_aspect("equal")
    fig.savefig("domain.png", dpi=130)
    print("wrote domain.png")
except Exception as e:
    print("plot failed", e)
