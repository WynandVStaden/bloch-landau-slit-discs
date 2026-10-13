"""Independent evaluation of U for design G (own geometry code; mpmath for the singular self integral, straight from Phi)."""
import json, sys, math, random
import numpy as np, mpmath as mp
sys.path.insert(0, "../../certify")
from inradius_certificate import DESIGNS
name = "G"
Rs, classes = DESIGNS[name]
R = float(Rs); L = math.log(R)
cert = json.load(open(f"../../certify/conformal_certificate_{name}.json"))
coefs = {int(k): np.array(v) for k, v in cert["coefficients"].items()}
slits = []  # (theta, a, b1, b2, class)
for c, (alpha, tip, b1, b2) in enumerate(classes):
    al = math.pi / 3 if alpha == "pi/3" else float(alpha)
    axial = alpha in ("0", "pi/3")
    for r in range(3):
        for sgn in ((1,) if axial else (1, -1)):
            slits.append((2 * math.pi * r / 3 + sgn * al, math.log(R / float(tip)), sgn * float(b1), sgn * float(b2), c))
print(len(slits), "slits; classes", [s[4] for s in slits])
def wj(j, x):
    th, a, b1, b2, _ = slits[j]; s = a * x
    return L + s + 1j * (th + b1 * s**2 + b2 * s**4)
def psi(c, y):  # numpy
    t = np.arccos(y); return sum(v * np.cos(n * t) for n, v in enumerate(coefs[c]) if v != 0)
def far(w, j, N):
    y = np.cos((2 * np.arange(1, N + 1) - 1) * np.pi / (2 * N))
    return np.sum(-np.log(np.abs(2 * np.sinh((w - wj(j, y)) / 2))) * psi(slits[j][4], y)) * np.pi / N
mp.mp.dps = 25
def self_mp(i, x):
    th, a, b1, b2, c = slits[i]
    cf = [mp.mpf(float(v)) for v in coefs[c]]
    def wm(xx):
        s = mp.mpf(a) * xx; return mp.mpf(L) + s + 1j * (mp.mpf(th) + mp.mpf(b1) * s**2 + mp.mpf(b2) * s**4)
    wx = wm(mp.mpf(x)); t0 = mp.acos(mp.mpf(x))
    def g(t):
        return -mp.log(abs(2 * mp.sinh((wx - wm(mp.cos(t))) / 2))) * sum(v * mp.cos(n * t) for n, v in enumerate(cf) if v != 0)
    return float(mp.quad(g, mp.linspace(0, t0, 5)[:-1] + mp.linspace(t0, mp.pi, 5)))
def U_on_slit(i, x, N=6000):
    w = wj(i, x)
    return (L - w.real) + self_mp(i, x) + sum(far(w, j, N) for j in range(len(slits)) if j != i)
def U_off(w, N=6000):
    return (L - w.real) + sum(far(w, j, N) for j in range(len(slits)))
# constant c
c_formula = L - math.pi / 4 * sum(s[1] * coefs[s[4]][1] for s in slits)
w = -30 + 0.7j
print("c: formula", repr(c_formula), " U + Re w at Re w=-30:", repr(U_off(w) + w.real), " exp(c) =", repr(math.exp(c_formula)))
random.seed(7)
print("U on circle (should be 0):", max(abs(U_off(L + 1j * random.uniform(0, 2 * math.pi), 20000)) for _ in range(8)))
print("U inside disc at a few points (must be > 0 if U ~ Green):", [round(U_off(math.log(r) + 1j * t), 6) for r, t in [(0.5, 0.3), (2.0, 0.6), (4.5, 0.3)]])
bound = {0: 2.656e-8, 1: 3.548e-9, 2: 1.758e-8, 3: 3.758e-9, 4: 7.478e-10}
worst = {c: 0.0 for c in bound}
for i in range(len(slits)):
    c = slits[i][4]
    xs = [-0.999999, -0.9999, -0.99, -0.5, -1e-3] + [-random.random() for _ in range(3)]
    vals = [U_on_slit(i, x) for x in xs]
    chk = U_on_slit(i, xs[5], 12000) - vals[5]
    worst[c] = max(worst[c], max(abs(v) for v in vals))
    print(f"slit {i:2d} class {c}: max|U| {max(abs(v) for v in vals):.3e} at x={xs[int(np.argmax(np.abs(vals)))]:.6f}  (certified sup {bound[c]:.3e}; N-doubling change {chk:.1e})", flush=True)
print("worst per class", worst, "all below certified:", all(worst[c] <= bound[c] for c in bound))
