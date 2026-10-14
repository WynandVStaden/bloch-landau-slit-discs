"""Geometry sanity, Green-function residuals (max principle a-posteriori bound), harmonic measures."""
import numpy as np
from geom import full_problem, slit_z, apply_image, all_slits
from methodA import solve, postprocess, eval_G, residual_on_slits

R, classes = full_problem()

# ---------------- geometry: pairwise distances between the 24 physical slits (fine polylines)
M = 4001
sl = all_slits(R, classes)
pts = []
for ic, img in sl:
    c = classes[ic]
    a = np.log(R / c['t'])
    pts.append(apply_image(slit_z(R, c, -a * np.linspace(0, 1, M)), img))
print('number of slits:', len(pts))
base_ang = sorted((np.angle(p[0]) % (2 * np.pi), sl[k][0] + 1) for k, p in enumerate(pts))
print('base angles (rad) and class:', [(round(a, 5), c) for a, c in base_ang])
gaps = np.diff([a for a, c in base_ang] + [base_ang[0][0] + 2 * np.pi])
print('min/max angular gap between consecutive bases:', gaps.min(), gaps.max())
dmin = np.inf
rows = []
for i in range(len(pts)):
    for j in range(i + 1, len(pts)):
        D = np.abs(pts[i][::8, None] - pts[j][None, ::8])
        k = np.unravel_index(np.argmin(D), D.shape)
        rows.append((D[k], sl[i][0] + 1, sl[j][0] + 1, abs(pts[i][::8][k[0]]), abs(pts[j][::8][k[1]])))
rows.sort()
print('closest slit pairs (distance, classes, |z| at closest points):')
for r in rows[:8]:
    print('   d=%.5f  classes %d-%d  |z|=%.4f, %.4f' % r)
# closest approach away from the outer circle (|z|<0.98R)
dd = []
for i in range(len(pts)):
    for j in range(i + 1, len(pts)):
        pi_, pj_ = pts[i][::8], pts[j][::8]
        pi_, pj_ = pi_[np.abs(pi_) < 0.98 * R], pj_[np.abs(pj_) < 0.98 * R]
        dd.append((np.abs(pi_[:, None] - pj_[None, :]).min(), sl[i][0] + 1, sl[j][0] + 1))
dd.sort()
print('closest pairs restricted to |z|<0.98R:', [(round(d, 4), a, b) for d, a, b in dd[:4]])
# self-intersection / monotonic radius is automatic (|z| = R e^s strictly monotone); slit tips:
for ic, c in enumerate(classes):
    a = np.log(R / c['t'])
    zt = slit_z(R, c, -a)
    print(f'class {ic + 1}: a={a:.6f} tip |z|={abs(zt):.12f} arg={np.angle(zt):+.6f}  base arg={c["theta"]:+.6f} '
          f'max|dv/ds|={np.max(np.abs(2 * c["b1"] * np.linspace(-a, 0, 999) + 4 * c["b2"] * np.linspace(-a, 0, 999)**3)):.4f}')

# ---------------- Green function residuals
for N in (48, 96, 192):
    f, C, A, rhs = solve(R, classes, N)
    o = postprocess(R, classes, N, f, C)
    res = residual_on_slits(R, classes, N, f, C)
    # circle |z|=R : G must vanish (holds through antisymmetry of f)
    th = np.linspace(0, 2 * np.pi / 3, 2000, endpoint=False) + 1e-3
    Gc = eval_G(R, classes, N, f, C, R * np.exp(1j * th))
    # pole normalisation: G(z)+log|z| -> log rho
    z0 = 1e-3 * np.exp(1j * np.array([0.3, 1.7, 4.0]))
    G0 = eval_G(R, classes, N, f, C, z0) + np.log(np.abs(z0))
    print(f'N={N:3d} rho={o["rho"]:.15f} max|G| on slits (off-node, both halves)={np.abs(res).max():.2e} per class='
          f'{np.abs(res).max(axis=1)}  max|G| on circle={np.abs(Gc).max():.2e}  G+log|z| near 0 - log rho={G0 - np.log(o["rho"])}')
    print(f'       => max-principle bound on |log rho_N - log rho| ~ {max(np.abs(res).max(), np.abs(Gc).max()):.2e},'
          f' i.e. |d rho| <~ {o["rho"] * max(np.abs(res).max(), np.abs(Gc).max()):.2e}')
print('positivity of density on physical half (f>0 for q>pi/2):', all((f[i][N // 2:] > 0).all() for i in range(5)))
print('harmonic measure from 0, per class (all slits, both sides):', o['hm_class'])
print('per single slit:', o['hm_per_slit'])
print('outer circle arcs:', 1 - o['hm_class'].sum())
# positivity of G inside D at random points (G>0 in D)
rng = np.random.default_rng(1)
zr = R * np.sqrt(rng.random(4000)) * np.exp(2j * np.pi * rng.random(4000))
Gr = eval_G(R, classes, N, f, C, zr)
print('min G at 4000 random interior points:', Gr.min(), ' (should be >= 0 up to near-slit quadrature error)')
