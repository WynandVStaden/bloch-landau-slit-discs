"""Extra independent checks of r(D) for the 24-slit domain."""
import numpy as np
from bie import conformal_radius, G0
from dom import R_G as R, slits_G, zslit
CL = 1.523726248805
S = slits_G()
sl = [t[1:] for t in S]
# (1) cube-map quotient: zeta = z^3 maps D (3-fold symmetric) onto disc radius R^3 minus 8 slits;
#     r(D)^3 = r(D^3).  In log coords s'=3s: alpha'=3alpha, b1'=b1/3, b2'=b2/27, tip'=tip^3.
q = {}
for c, al, tip, b1, b2 in S:
    key = (round(np.cos(3 * al), 9), round(np.sin(3 * al), 9), c)
    q[key] = (3 * al, tip ** 3, b1 / 3, b2 / 27)
qs = list(q.values())
print("quotient domain slits:", len(qs), " R^3 =", R ** 3)
for n, p in [(8, 16), (12, 16), (16, 16), (12, 20)]:
    r3 = conformal_radius(R ** 3, qs, n=n, p=p, m=2)
    print(f"  cube quotient n={n} p={p}: r^(1/3)={r3 ** (1/3):.14f}  minus claim {r3 ** (1/3) - CL:+.3e}")
# (2) residual of the full solve: g at interior points, boundary values between nodes, r from g near 0
r, info = conformal_radius(R, sl, n=12, p=16, m=2, return_all=True)
Z, sig = info["Z"], info["sig"]; W = np.tile(info["w"], len(sl))
print("full solve r = %.14f ; total charge (harmonic measure of the slits from 0) = %.12f" % (r, np.sum(sig * W)))
def g(z):
    return -np.log(abs(z)) + np.log(R) - np.sum(G0(z, Z, R) * sig * W)
for z in (1e-3, 1e-3j, -1e-4 + 2e-4j, 0.3 * np.exp(1j)):
    print("  h(z)=g+log|z| at z=%s: %.13f  (log r = %.13f; harmonic, equal only as z->0)" % (z, g(z) + np.log(abs(z)), np.log(r)))
rng = np.random.default_rng(7)
mn = 1e9; cnt = 0
for _ in range(40000):
    z = R * np.sqrt(rng.random()) * np.exp(2j * np.pi * rng.random())
    if np.abs(Z - z).min() < 0.04 or abs(z) < 1e-3: continue
    v = g(z); mn = min(mn, v); cnt += 1
print("  min of g over %d random interior points >=0.04 from slits: %.3e (must be > 0)" % (cnt, mn))
# off-node boundary residual: evaluate g on slit points 0.5 from other slits' influence -> use points on a
# parallel offset curve at tiny normal distance is near-singular; instead check mean value property:
# log r = mean of h on small circle
th = np.linspace(0, 2 * np.pi, 64, endpoint=False)
for rad in (0.2, 0.5, 0.9):
    hm = np.mean([g(rad * np.exp(1j * t)) + np.log(rad) for t in th])
    print("  mean of g+log|z| on |z|=%.1f: exp = %.13f" % (rad, np.exp(hm)))
# narrow channel between A (0 deg) and C: g must be positive & tiny there, and positive between B and D
for (i, j) in [(0, 6), (3, 12)]:
    for x in (-0.3, -0.1, -0.03, -0.01):
        zc = 0.5 * (zslit(x, R, *sl[i]) + zslit(x * np.log(R / sl[i][1]) / np.log(R / sl[j][1]), R, *sl[j]))
        print(f"  channel {S[i][0]}-{S[j][0]} midpoint |z|={abs(zc):.4f}: gap={abs(zslit(x, R, *sl[i]) - zc)*2:.4f} g={g(zc):.3e}")
# (3) sensitivity: how big a perturbation of r would be needed
print("needed r >= %.9f to beat 0.656319277272; margin = %.3e" % (1 / 0.656319277272, r - 1 / 0.656319277272))
print("1/r = %.12f" % (1 / r))
