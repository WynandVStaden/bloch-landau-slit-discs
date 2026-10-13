"""Referee checks for the paper (read-only w.r.t. the rest of the repository)."""
import json, sys, math
from fractions import Fraction
import numpy as np
sys.path.insert(0, r"C:\Users\wvs00\claude_chat\bloch_record\certify")
sys.path.insert(0, r"C:\Users\wvs00\claude_chat\bloch_record")
from flint import arb, ctx
ctx.prec = 200
import inradius_certificate as ic
CERT = r"C:\Users\wvs00\claude_chat\bloch_record\certify"
# 1 arithmetic
lam = 1 + Fraction(1, 10**9)
lo = Fraction("1.523726208329")
print("theorem ratio", float(lam / lo), lam / lo < Fraction("0.65628589"))
print("is double(1.523726208329) < decimal?", Fraction(1.523726208329) < lo, float(lo - Fraction(1.523726208329)))
print("B ratio", float(lam / Fraction("1.523673641792")))
print("float ratio 1/1.523726248805 =", 1 / 1.523726248805, " proved - float =", float(lam / lo) - 1 / 1.523726248805,
      " stated(0.65628589) - float =", 0.65628589 - 1 / 1.523726248805)
# 2 density file
d = json.load(open(CERT + r"\density_G_n64.json"))
for c, v in d["coefficients"].items():
    nz = [i for i, x in enumerate(v) if x != 0]
    print("class", c, "len", len(v), "max even coef", max(abs(x) for x in v[0::2]), "max index", nz[-1], "sum|c|", sum(abs(x) for x in v), "c1", v[1])
# 3 exact c and exp(c - delta)
big_r, classes = ic.DESIGNS["G"]
table = ic.slit_table(classes)
print("slits:", len(table), "distinct angles:", len({(round(float(ic.float_angle(r, s, a)) % (2 * math.pi), 9), float(t)) for r, s, a, t, _, _ in table}))
tips = [Fraction(c[1]) for c in classes]
R = ic.arb_frac(Fraction(big_r))
c1 = arb(0)
for r, sgn, al, tip, b1, b2 in table:
    c1 += (R / ic.arb_frac(tip)).log() * arb(d["coefficients"][str(tips.index(tip))][1])
c = R.log() - arb.pi() / 4 * c1
print("c =", c.str(20), " exp(c) =", c.exp().str(20))
cert = json.load(open(CERT + r"\conformal_certificate_G.json"))
dl = cert["delta_upper"]
print("delta_upper", dl, " exp(c-delta) =", (c - arb(dl)).exp().str(20))
# 4 geometry: min distance between slits, class 0-2 and 1-3 distance in w-plane
def wpts(row, n=4001, ext=False):
    r, sgn, al, tip, b1, b2 = row
    a = math.log(float(Fraction(big_r)) / float(tip))
    s = a * np.linspace(-1, 1 if ext else 0, n)
    return math.log(float(Fraction(big_r))) + s + 1j * (ic.float_angle(r, sgn, al) + float(b1) * s**2 + float(b2) * s**4)
W = [wpts(t) for t in table]
Z = [np.exp(w) for w in W]
best = (9, None)
for i in range(24):
    for j in range(i + 1, 24):
        dd = np.abs(Z[i][:, None] - Z[j][None, ::4]).min()
        if dd < best[0]: best = (dd, (i, j))
print("min z-distance between slits", best)
def wdist(i, j):
    dd = W[i][:, None] - W[j][None, ::4]
    im = (dd.imag + np.pi) % (2 * np.pi) - np.pi
    return np.hypot(dd.real, im).min()
print("w-dist class0-class2", wdist(0, 6), " class1-class3", wdist(3, 12), "class 4 nearest", min(wdist(18, j) for j in range(24) if j != 18))
# 5 solver: condition number and closed-form accuracy
import bloch_curved2 as bc
from bloch_residual import solve
from bloch_solver import phi
fl = [bc.Slit2(ic.float_angle(r, s, a), math.log(float(Fraction(big_r)) / float(t)), float(b1), float(b2)) for r, s, a, t, b1, b2 in table]
for n in (64, 96):
    # rebuild matrix as in solve()
    big_l = math.log(float(Fraction(big_r)))
    k = np.arange(1, n + 1); x = np.cos((2 * k - 1) * np.pi / (2 * n))
    psi, rows, lr = solve(float(Fraction(big_r)), fl, n)
    A = np.vstack([rows(i, x) for i in range(24)])
    print("n", n, "size", A.shape, "cond", np.linalg.cond(A), "rho", math.exp(lr))
