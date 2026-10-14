"""Validation of Method B (and of the explicit two-slit map) on exact cases."""
import numpy as np
from methodB import solve, Arc, rot, two_slit_map, arcs_B0, log_rule, gl01, bary_matrix
from geom import R0, CLASSES0

# 0. unit test of the product log rule
for x in (0.3, 0.0001, 0.9999, -0.01, -1e-4, 1.2, 1.0005):
    u, om = log_rule(x, 0.0, 1.0)
    num = om @ np.cos(u)
    import mpmath as mp
    pts = [0, x, 1] if 0 < x < 1 else [0, 1]
    ex = mp.quad(lambda y: mp.log(abs(x - y)) * mp.cos(y), np.linspace(0, 1, 2).tolist() if not (0 < x < 1) else pts, maxdegree=10) \
        if not (0 < x < 1 or abs(x) < 0.02 or abs(x - 1) < 0.02) else None
    if ex is None:
        mp.mp.dps = 30
        ex = mp.quad(lambda y: mp.log(abs(x - y)) * mp.cos(y), [0, x, 1] if 0 < x < 1 else [0, 1e-3, 1e-2, 0.1, 0.9, 0.99, 0.999, 1])
    print(f'log rule x={x:8.4f}: err={float(num - ex):+.2e}')


def run(name, arcs, back, exact, cfgs=((2, 8), (4, 8), (4, 12), (8, 16))):
    for P, q in cfgs:
        rho_eta, hm, info = solve(arcs, P, q)
        print(f'{name:38s} P={P} q={q:2d} rho={back(rho_eta):.15f} exact={exact:.15f} err={back(rho_eta) - exact:+.2e}')
    return hm


for s in (0.5, 0.1, 0.01):
    hm = run(f'unit disc minus [{s},1)', [Arc(np.log(1 / s), 0.0, 0.0, 0.0, [lambda z: z])], lambda r: r, 4 * s / (1 + s)**2)
    k = 4 * s / (1 + s)**2
    print(f'    harmonic measure of slit {hm[0]:.14f} exact {np.arccos(2 * k - 1) / np.pi:.14f}')
R = R0
for t in (1.0, 3.0):
    s = (t / R)**3
    arcs, back = arcs_B0(R, [(0.0, t, 0.0, 0.0)])
    run(f'R-disc minus 3 slits [{t},R)', arcs, back, R * (4 * s / (1 + s)**2)**(1 / 3))
t1, t2 = CLASSES0[0][1], CLASSES0[1][1]
p, q = (t1 / R)**3, (t2 / R)**3
Aa, Bb = (p + 1 / p) / 2, (q + 1 / q) / 2
arcs, back = arcs_B0(R, CLASSES0[:2])
run('R-disc minus 6 radial slits (cl.1+2)', arcs, back, R * (2 / (Aa + Bb))**(1 / 3))

# explicit map checks: boundary of the two-slit disc goes to the unit circle, F(z)/z -> (A+B)/2
F, F0 = two_slit_map(p, q)
th = np.linspace(0.01, np.pi - 0.01, 7)
print('|F| on circle      :', np.abs(F(np.exp(1j * th))[0]))
xs = np.array([p + 0.3 * (1 - p), p * 1.0001, 0.999]) + 1e-15j
print('|F| on slit [p,1)  :', np.abs(F(xs)[0]))
xs = -np.array([q + 0.3 * (1 - q), q * 1.0001, 0.999]) + 1e-15j
print('|F| on slit (-1,-q]:', np.abs(F(xs)[0]))
z = 1e-7 * np.exp(0.3j)
print('F(z)/z at 0, (A+B)/2:', F(np.array([z]))[0][0] / z, F0)
z = np.array([0.3 + 0.2j]); hh = 1e-6
print('F\' check:', (F(z + hh)[0] - F(z - hh)[0]) / (2 * hh), F(z)[1])
