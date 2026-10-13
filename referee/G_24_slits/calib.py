import sys, time
import numpy as np
from bie import conformal_radius
from dom import R_18, REF_18, slits_18
d2r = np.pi / 180
print("== cal1: unit disc minus [s,1]")
for s in (0.5, 0.1, 0.9, 0.01):
    ex = 4 * s / (1 + s) ** 2
    for n in (4, 8):
        r = conformal_radius(1.0, [(0.0, s, 0.0, 0.0)], n=n)
        print(f"  s={s} n={n}: r={r:.15f} exact={ex:.15f} err={r-ex:+.2e}")
print("== cal2: radius 3, slits [1,3] x3, [2,3] x3; exact 1.5207430363125")
sl = [(k * 120 * d2r, 1.0, 0, 0) for k in range(3)] + [((60 + k * 120) * d2r, 2.0, 0, 0) for k in range(3)]
for n, p, m in [(4, 16, 1), (8, 16, 1), (8, 16, 2), (12, 16, 1)]:
    r = conformal_radius(3.0, sl, n=n, p=p, m=m)
    print(f"  n={n} p={p} m={m}: r={r:.14f} err={r-1.5207430363125:+.2e}")
print("== cal3: curved 18-slit domain, ref", REF_18)
sl = [t[1:] for t in slits_18()]
for n, p, m, pw in [(6, 16, 1, 2), (10, 16, 1, 2), (10, 16, 2, 2), (10, 12, 2, 2), (14, 16, 1, 2), (10, 16, 2, 3)]:
    t0 = time.time()
    r = conformal_radius(R_18, sl, n=n, p=p, m=m, power=pw)
    print(f"  n={n} p={p} m={m} power={pw}: r={r:.14f} err={r-REF_18:+.2e} ({time.time()-t0:.0f}s)", flush=True)
