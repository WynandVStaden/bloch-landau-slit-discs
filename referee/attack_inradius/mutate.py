# Mutation tests of the authors' stage-2 verifier (run on a COPY) and of the referee verifier.
import sys, numpy as np
from fractions import Fraction
sys.path.insert(0, "run_copy")
import inradius_certificate_copy as P
import verify_leaves as V
z = np.load("stored_leaves_G.npz"); CELLS, PARS = z["cells"], z["pars"]
BIG_R = Fraction(P.DESIGNS["G"][0]); SLITS = P.slit_table(P.DESIGNS["G"][1])
BAD = "   ###### UNEXPECTED"
def authors(cells, pars, eps="1e-9", slits=SLITS, big_r=BIG_R):
    try:
        P.no_overlap(cells); P.verify(big_r, slits, Fraction(eps), cells, pars); return "ACCEPT"
    except BaseException as e:
        return "reject(%s: %s)" % (type(e).__name__, str(e)[:50])
def referee(cells, pars, eps="1e-9", T=None):
    np.savez("tmp_mut.npz", cells=cells, pars=pars, eps=1e-9)
    try:
        V.main("tmp_mut.npz", eps, T=T, quiet=True); return "ACCEPT"
    except BaseException as e:
        return "reject(%s: %s)" % (type(e).__name__, str(e)[:50])
def show(name, expect, a, r=None):
    ok = (a.startswith("ACCEPT") == (expect == "ACCEPT")) and (r is None or (r.startswith("ACCEPT") == (expect == "ACCEPT")))
    print("%-60s expect %-6s authors: %-66s referee: %s%s" % (name, expect, a, r, "" if ok else BAD), flush=True)
show("baseline stored leaves, eps=1e-9", "ACCEPT", authors(CELLS, PARS), referee(CELLS, PARS))
# ---- EPS mutations on the stored leaves
for e, exp in [("0", "reject"), ("-1e-9", "reject"), ("1e-10", "reject"), ("5e-10", "?"), ("9e-10", "?"), ("-1e-30", "reject"), ("-1", "reject"), ("-2", "reject"), ("-2.0000001", "reject"), ("-3", "reject"), ("-100", "reject")]:
    a = authors(CELLS, PARS, e)
    r = referee(CELLS, PARS, e)
    show("stored leaves, eps=%s" % e, exp if exp != "?" else ("ACCEPT" if a == "ACCEPT" else "reject"), a, r)
# ---- domain mutations checked against the stored leaves
def design(mod):
    R, cls = P.DESIGNS["G"]; cls = [list(c) for c in cls]; R, cls = mod(R, cls); return Fraction(R), P.slit_table([tuple(c) for c in cls]), cls
muts = {
 "class0 tip 1 -> 1.000001": lambda R, c: (R, [c[0][:1] + ["1.000001"] + c[0][2:]] + c[1:]),
 "class0 tip 1 -> 1.00000001": lambda R, c: (R, [c[0][:1] + ["1.00000001"] + c[0][2:]] + c[1:]),
 "class1 tip -> 2.00000001": lambda R, c: (R, c[:1] + [c[1][:1] + ["2.00000001"] + c[1][2:]] + c[2:]),
 "class2 tip -> 2.78934242 (+1e-8)": lambda R, c: (R, c[:2] + [c[2][:1] + ["2.78934242"] + c[2][2:]] + c[3:]),
 "class3 tip -> 3.47338126 (+1.2e-8)": lambda R, c: (R, c[:3] + [c[3][:1] + ["3.47338126"] + c[3][2:]] + c[4:]),
 "class4 tip -> 4.3428742 (+1.5e-8)": lambda R, c: (R, c[:4] + [c[4][:1] + ["4.3428742"] + c[4][2:]]),
 "class4 tip -> 4.9": lambda R, c: (R, c[:4] + [c[4][:1] + ["4.9"] + c[4][2:]]),
 "class4 removed": lambda R, c: (R, c[:4]),
 "class2 removed": lambda R, c: (R, c[:2] + c[3:]),
 "R -> 5.0012": lambda R, c: ("5.0012", c),
 "class2 b1 -> 0.9027": lambda R, c: (R, c[:2] + [c[2][:2] + ["0.9027"] + c[2][3:]] + c[3:]),
 "class3 theta -> 1.0158": lambda R, c: (R, c[:3] + [["1.0158"] + c[3][1:]] + c[4:]),
}
for name, mod in muts.items():
    R, sl, cls = design(mod)
    show("stored leaves vs mutated domain: " + name, "reject", authors(CELLS, PARS, "1e-9", sl, R), None)
# ---- tampered leaves
rng = np.random.default_rng(1)
i2 = np.nonzero(CELLS[:, 3] == 2)[0]; i0 = np.nonzero(CELLS[:, 3] == 0)[0]; i1 = np.nonzero(CELLS[:, 3] == 1)[0]
def T(name, f, n=1, expect="reject"):
    acc_a = acc_r = 0; last = ("", "")
    for _ in range(n):
        c, p = CELLS.copy(), PARS.copy(); out = f(c, p)
        if out is not None: c, p = out
        a = authors(c, p); r = referee(c, p); acc_a += a == "ACCEPT"; acc_r += r == "ACCEPT"
        if a != "ACCEPT": last = (a, r)
    flag = "" if (acc_a == 0 and acc_r == 0) or expect != "reject" else BAD
    print("%-60s trials %3d  authors accepted %3d  referee accepted %3d   last reject: %s | %s%s" % (name, n, acc_a, acc_r, last[0][:60], last[1][:60], flag), flush=True)
T("delete one leaf", lambda c, p: (np.delete(c, rng.integers(len(c)), 0), np.delete(p, 0)), 5)
def f(c, p): c[rng.choice(i2), 3] = 0
T("witness leaf relabelled outside-disc", f, 150)
def f(c, p): c[rng.choice(i2), 3] = 1
T("witness leaf relabelled outside-sector", f, 150)
def f(c, p): c[rng.choice(i0), 3] = 1
T("outside-disc leaf relabelled outside-sector (may be true)", f, 150, "some")
def f(c, p): c[rng.choice(i1), 3] = 0
T("outside-sector leaf relabelled outside-disc (may be true)", f, 150, "some")
def f(c, p):
    k = rng.choice(i2[CELLS[i2, 4] == 0]); p[k] = abs(p[k]) + 1e-3
T("slit witness parameter made positive", f, 20)
def f(c, p):
    k = rng.choice(i2[CELLS[i2, 4] == 0]); p[k] = -2.0
T("slit witness parameter below -a (beyond the tip)", f, 20)
def f(c, p):
    k = rng.choice(i2[CELLS[i2, 4] == 0]); p[k] = np.nan
T("slit witness parameter NaN", f, 5)
def f(c, p):
    k = rng.choice(i2[CELLS[i2, 4] == 2]); p[k] = np.inf
T("circle witness parameter inf", f, 3)
def f(c, p):
    k = rng.choice(i2); c[k, 5] = (c[k, 5] + rng.integers(1, 24)) % 24
T("witness slit index changed (may be true for circle kind)", f, 150, "some")
def f(c, p):
    k = rng.choice(i2[CELLS[i2, 4] != 2]); c[k, 5] = 24
T("witness slit index out of range (24)", f, 3)
def f(c, p):
    k = rng.choice(i2[CELLS[i2, 4] != 2]); c[k, 5] = -1 - rng.integers(0, 24)
T("witness slit index negative (may hit a valid slit)", f, 100, "some")
def f(c, p):
    k = rng.choice(i2); c[k, 4] = (c[k, 4] + rng.integers(1, 3)) % 3
T("witness kind changed slit/tip/circle (may be true)", f, 200, "some")
def f(c, p):
    k = rng.choice(i2); c[k, 4] = rng.choice([3, -1, 7])
T("witness kind invalid value (authors treat as slit point)", f, 100, "some")
def f(c, p):
    k = rng.choice(i2); c[k, 3] = rng.choice([3, -1, 5])
T("code invalid value", f, 50)
def f(c, p):
    k = rng.choice(len(c)); c[k, 1] += rng.choice([-1, 1])
T("leaf moved sideways by one cell (duplicate + gap)", f, 100)
def f(c, p):
    for _ in range(10000):
        k = rng.choice(i2); lv, ix, iy = c[k, :3]
        sib = [np.nonzero((c[:, 0] == lv) & (c[:, 1] == (ix | 1) - dx) & (c[:, 2] == (iy | 1) - dy))[0] for dx in (0, 1) for dy in (0, 1)]
        if all(len(s) == 1 for s in sib) and lv > 6: break
    sib = [int(s[0]) for s in sib]
    row = c[k].copy(); row[0] -= 1; row[1] >>= 1; row[2] >>= 1
    keep = np.ones(len(c), bool); keep[sib] = False
    return np.vstack([c[keep], row[None]]), np.concatenate([p[keep], p[k:k+1]])
T("4 siblings merged into parent, one child's witness (may be true)", f, 200, "some")
def f(c, p):
    k = rng.choice(np.nonzero(c[:, 0] > 6)[0]); big = rng.choice(np.nonzero(c[:, 0] < c[k, 0])[0])
    d = int(c[k, 0] - c[big, 0]); c[k, 1] = (int(c[big, 1]) << d) + int(rng.integers(0, 1 << min(d, 20))); c[k, 2] = (int(c[big, 2]) << d)
T("leaf moved inside a larger leaf (area preserved)", f, 100)
def f(c, p):
    k = rng.choice(i2[CELLS[i2, 0] < 20])
    lv, ix, iy, code, kind, j = c[k]
    kids = np.array([[lv + 1, 2*ix + dx, 2*iy + dy, 2, 2, 0] for dx in (0, 1) for dy in (0, 1)])
    return np.vstack([np.delete(c, k, 0), kids]), np.concatenate([np.delete(p, k), np.full(4, 3.0)])
T("leaf split, children given circle witness at angle 3.0", f, 30)
def f(c, p): c[rng.choice(len(c)), 0] = rng.choice([5, 0, -1, 90, 200])
T("level set to illegal value", f, 30)
def f(c, p):
    k = rng.choice(i2[CELLS[i2, 4] == 0]); p[k] = p[k]*(1 + rng.choice([-1, 1])*10.0**rng.integers(-9, -2))
T("slit witness parameter perturbed rel 1e-9..1e-3 (may be true)", f, 200, "some")
