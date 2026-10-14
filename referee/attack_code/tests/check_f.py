"""Independent mpmath evaluation of f_0(x) = U(w_0(x)) for design G, straight from definition (eq:U),
compared with the program's ball enclosures at the same nodes. Run from tests/ ; imports the COPY in ../certify."""
import sys, os, json, random, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "certify"))
import conformal_certificate as cc
from flint import arb, ctx
from mpmath import mp, mpf, mpc, pi, log, sinh, cos, acos, quad, fabs
mp.dps = 22
REP = int(sys.argv[1]) if len(sys.argv) > 1 else 0
NPTS = int(sys.argv[2]) if len(sys.argv) > 2 else 20
# ---- geometry written out independently from the paper's Table 1
R = mpf("5.001117628756668")
TAB = [("0", "1", "0", "0"), ("pi/3", "1.9999999974092333", "0", "0"),
       ("0.029930188727091613", "2.789342409403121", "0.9026921410641155", "-0.3301459494660587"),
       ("1.0157094972537686", "3.473381248030611", "-1.2656124700253997", "0.9794269477754911"),
       ("0.5757776442359243", "4.342874184977691", "0.03225521024585182", "2.9360081049550737")]
slits = []  # (class, theta, a, b1, b2)
for c, (th, t, b1, b2) in enumerate(TAB):
    th = pi / 3 if th == "pi/3" else mpf(th)
    signs = (1,) if c < 2 else (1, -1)
    for r in range(3):
        for e in signs:
            slits.append((c, 2 * pi * r / 3 + e * th, log(R / mpf(t)), e * mpf(b1), e * mpf(b2)))
assert len(slits) == 24
L = log(R)
dens = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "certify", "density_G_n64.json")))["coefficients"]
dens = {int(k): [mpf(v) for v in vals] for k, vals in dens.items()}
def wj(j, x):
    c, th, a, b1, b2 = slits[j]
    s = a * x
    return mpc(L + s, th + b1 * s**2 + b2 * s**4)
def psi_theta(c, th):
    return sum(v * cos(n * th) for n, v in enumerate(dens[c]) if v != 0)
def Phi(w):
    return -log(abs(2 * sinh(w / 2)))
def U_at(i, x):
    w = wj(i, x)
    tot = L - w.real
    th0 = acos(x)
    for j in range(24):
        c = slits[j][0]
        g = lambda th: Phi(w - wj(j, cos(th))) * psi_theta(c, th)
        if j == i:
            pts = sorted(set([mpf(0), th0 / 2, th0 - mpf("0.05"), th0, th0 + mpf("0.05"), (th0 + pi) / 2, pi]))
            pts = [p for p in pts if 0 <= p <= pi]
        else:
            pts = [pi * k / 12 for k in range(13)]
        v, err = quad(g, pts, error=True, maxdegree=9)
        tot += v
        U_at.err = max(U_at.err, err)
    return tot
U_at.err = 0
# ---- the program's enclosures
ctx.prec = cc.PREC
geo = cc.Geometry("G")
coef = {int(c): v for c, v in json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "certify", "density_G_n64.json")))["coefficients"].items()}
# geometry cross-check: program's table vs mine
i_rep = geo.reps[REP]
for j in range(24):
    assert geo.cls[j] == slits[j][0]
    for xx in (-1.0, -0.37, 0.61):
        d = complex(geo.w(j, arb(xx))) - complex(wj(j, mpf(xx)))
        assert abs(d) < 1e-14, (j, d)
print("geometry of the program agrees with the table of the paper (24 arcs, 3 sample points each)")
t0 = time.time()
plan_log = []
pl = cc.plan(geo, coef, plan_log.append)
print("plan:", plan_log[REP].strip(), f"[{time.time()-t0:.0f}s]")
M = pl[i_rep]["M"]
random.seed(12345 + REP)
ks = sorted(random.sample(range(M // 2, M), NPTS - 2)) + [M // 2, M - 1]
cc._init("G", coef)
_, res = cc._eval_nodes((i_rep, M, ks, pl[i_rep]["quad"]))
worst = 0
print(" k      x_k           program mid      program rad   mpmath U         diff       inside?  (quad err est)")
for k, mid, rad in res:
    x = cos(pi * (2 * k + 1) / (2 * M))
    U_at.err = 0
    t0 = time.time()
    u = U_at(i_rep, x)
    diff = float(u - mpf(mid))
    ok = abs(diff) <= rad + float(U_at.err)
    worst = max(worst, abs(float(u)))
    print(f"{k:5d} {float(x):+.9f}  {mid:+.6e}  {rad:.2e}  {float(u):+.6e}  {diff:+.2e}  {ok}  ({float(U_at.err):.1e}) [{time.time()-t0:.0f}s]", flush=True)
print("max |U| at sampled points (mpmath):", worst)
