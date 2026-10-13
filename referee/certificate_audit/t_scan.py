"""Dense float scan of f_i on the representative slits of G (own code), compared with the certified sup."""
import json, sys, math
import numpy as np
sys.path.insert(0, "../../certify")
from inradius_certificate import DESIGNS
Rs, classes = DESIGNS["G"]; R = float(Rs); L = math.log(R)
coefs = {int(k): np.array(v) for k, v in json.load(open("../../certify/conformal_certificate_G.json"))["coefficients"].items()}
slits = []
for c, (alpha, tip, b1, b2) in enumerate(classes):
    al = math.pi / 3 if alpha == "pi/3" else float(alpha)
    for r in range(3):
        for sgn in ((1,) if alpha in ("0", "pi/3") else (1, -1)):
            slits.append((2 * math.pi * r / 3 + sgn * al, math.log(R / float(tip)), sgn * float(b1), sgn * float(b2), c))
def wj(j, x):
    th, a, b1, b2, _ = slits[j]; s = a * x
    return L + s + 1j * (th + b1 * s**2 + b2 * s**4)
N = 5000
tq = (2 * np.arange(1, N + 1) - 1) * np.pi / (2 * N); yq = np.cos(tq)
psi_q = {c: sum(v * np.cos(n * tq) for n, v in enumerate(cf) if v != 0) for c, cf in coefs.items()}
bound = {0: 2.656e-8, 1: 3.548e-9, 2: 1.758e-8, 3: 3.758e-9, 4: 7.478e-10}
for i in [0, 3, 6, 12, 18] + [9, 15]:
    th, a, b1, b2, c = slits[i]
    tx = np.concatenate([np.linspace(np.pi / 2, np.pi, 2500), np.arccos(-np.random.default_rng(0).random(500))]); x = np.cos(tx)
    f = -a * x + sum(v * np.pi / n * np.cos(n * tx) for n, v in enumerate(coefs[c]) if v != 0)
    wx = wj(i, x)
    for j in range(len(slits)):
        for s0 in range(0, len(x), 500):
            xs = x[s0:s0 + 500, None]; ws = wx[s0:s0 + 500, None]
            if j == i:
                dd = a + 1j * (b1 * a**2 * (xs + yq) + b2 * a**4 * (xs**3 + xs**2 * yq + xs * yq**2 + yq**3))
                u = (xs - yq) * dd / 2
                k = dd * np.where(np.abs(u) < 1e-4, 1 + u * u / 6, np.sinh(u) / np.where(u == 0, 1, u))
            else:
                k = 2 * np.sinh((ws - wj(j, yq)) / 2)
            f[s0:s0 + 500] -= (np.log(np.abs(k)) * psi_q[slits[j][4]]).sum(axis=1) * np.pi / N
    m = int(np.argmax(np.abs(f)))
    print(f"slit {i:2d} class {c}: sup|f| over {len(x)} points = {abs(f[m]):.4e} at x = {x[m]:.5f}; max f = {f.max():.3e}, min f = {f.min():.3e}; certified {bound[c]:.3e} -> {'ok' if abs(f[m]) <= bound[c] else '*** EXCEEDS ***'}", flush=True)
