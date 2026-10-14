"""Method C (sanity only, ~3 digits): walk-on-spheres.  log rho = E_0[ log|B_tau| ]  (B_tau = exit point of
Brownian motion from 0), since G(z,0) = -log|z| + h with h = log|z| on the boundary.  Also exit statistics
(harmonic measure) per slit class.  Uses only the problem statement (polyline slits), no integral equation."""
import sys, time
import numpy as np
from geom import full_problem, slit_z, apply_image, all_slits

R, classes = full_problem()
NSEG = 120


def build(slit_list):
    A, B, lab = [], [], []
    for ic, img in slit_list:
        c = classes[ic]
        a = np.log(R / c['t'])
        n = 1 if (c['b1'] == 0 and c['b2'] == 0) else NSEG
        p = apply_image(slit_z(R, c, -a * np.linspace(0, 1, n + 1)), img)
        A += list(p[:-1]); B += list(p[1:]); lab += [ic] * n
    return np.array(A), np.array(B), np.array(lab)


def dist(P, A, B):
    AB = B - A
    t = np.clip(((P[:, None] - A[None, :]) * np.conj(AB)[None, :]).real / (np.abs(AB)**2)[None, :], 0, 1)
    Q = A[None, :] + t * AB[None, :]
    D = np.abs(P[:, None] - Q)
    k = np.argmin(D, axis=1)
    idx = np.arange(len(P))
    return D[idx, k], Q[idx, k], k


def fold(z):
    """map to fundamental sector 0<=arg<=pi/3 by the D3 symmetry"""
    th = np.mod(np.angle(z), 2 * np.pi / 3)
    th = np.where(th > np.pi / 3, 2 * np.pi / 3 - th, th)
    return np.abs(z) * np.exp(1j * th)


full = build(all_slits(R, classes))
# slits relevant for the fundamental sector: cl1 (0), cl3, conj cl3, cl5, cl4, cl2 (pi/3), mirror of cl4
local = build([(0, (0.0, False)), (2, (0.0, False)), (2, (0.0, True)), (4, (0.0, False)), (3, (0.0, False)),
               (1, (0.0, False)), (3, (2 * np.pi / 3, True))])

if __name__ == '__main__':
    ntot = int(float(sys.argv[1])) if len(sys.argv) > 1 else 200000
    eps = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-6
    rng = np.random.default_rng(int(sys.argv[3]) if len(sys.argv) > 3 else 12345)
    # test of the folding shortcut against the full 24-slit distance
    zt = R * np.sqrt(rng.random(5000)) * np.exp(2j * np.pi * rng.random(5000))
    d1 = dist(zt, *full[:2])[0]
    d2 = dist(fold(zt), *local[:2])[0]
    print('fold-distance test, max |d_full - d_fold| =', np.abs(d1 - d2).max(), flush=True)
    batch = 20000
    vals, labs = [], []
    t0 = time.time()
    done = 0
    while done < ntot:
        z = np.zeros(batch, dtype=complex)
        out = np.full(batch, np.nan)
        lab = np.full(batch, -1)
        act = np.arange(batch)
        steps = 0
        while len(act):
            zf = fold(z[act])
            ds, Q, k = dist(zf, local[0], local[1])
            dc = R - np.abs(zf)
            d = np.minimum(ds, dc)
            hit = d < eps
            if hit.any():
                hs = hit & (ds <= dc)
                hc = hit & (ds > dc)
                out[act[hs]] = np.log(np.abs(Q[hs]))
                lab[act[hs]] = local[2][k[hs]]
                out[act[hc]] = np.log(R)
                lab[act[hc]] = 5
            keep = ~hit
            act = act[keep]
            z[act] = zf[keep] + d[keep] * np.exp(2j * np.pi * rng.random(len(act)))
            steps += 1
        vals.append(out); labs.append(lab)
        done += batch
        v = np.concatenate(vals); l = np.concatenate(labs)
        m, se = v.mean(), v.std() / np.sqrt(len(v))
        print(f'n={done:7d} steps={steps:5d} log rho={m:.5f}+-{se:.5f}  rho={np.exp(m):.5f}+-{np.exp(m) * se:.5f}  '
              f'hm classes 1..5, circle: {np.round([np.mean(l == i) for i in range(6)], 4)}  t={time.time() - t0:.0f}s', flush=True)
