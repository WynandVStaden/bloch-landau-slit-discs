"""Method A in extended precision (mpmath, 34 digits): matrix assembled in mp arithmetic, linear system solved
by float64-LU-preconditioned iterative refinement with mp residuals.  Independent re-coding of methodA.py."""
import sys, time
import numpy as np
import mpmath as mp
from scipy.linalg import lu_factor, lu_solve
from geom import CLASSES0, R0

mp.mp.dps = 34
R = mp.mpf(R0)
cls = [(mp.mpf(th), mp.mpf(t), mp.mpf(b1), mp.mpf(b2)) for th, t, b1, b2 in CLASSES0]
cls[1] = (mp.pi / 3, cls[1][1], mp.mpf(0), mp.mpf(0))
nimg = [3, 3, 6, 6, 6]
logabs = np.frompyfunc(lambda d: mp.log(abs(d)), 1, 1)
tofloat = np.frompyfunc(float, 1, 1)


def run(N):
    t0 = time.time()
    q = [(mp.mpf(k) + mp.mpf(1) / 2) * mp.pi / N for k in range(N)]
    cq = np.array([mp.cos(x) for x in q], dtype=object)
    Cm = np.array([[mp.cos(m * x) for m in range(1, N)] for x in q], dtype=object)     # (N, N-1)
    Cs = Cm * np.array([mp.mpf(2) / m for m in range(1, N)], dtype=object)[None, :]
    h = mp.pi / N
    W = (Cs @ Cm.T + mp.log(2)) * h
    zs, dzs, aa = [], [], []
    for th, t, b1, b2 in cls:
        a = mp.log(R / t)
        s = a * cq
        z = np.array([R * mp.exp(si + 1j * (th + b1 * si**2 + b2 * si**4)) for si in s], dtype=object)
        dz = np.array([zi * (1 + 1j * (2 * b1 * si + 4 * b2 * si**3)) * a for zi, si in zip(z, s)], dtype=object)
        zs.append(z); dzs.append(dz); aa.append(a)
    rots = [mp.expjpi(mp.mpf(2 * r) / 3) for r in range(3)]
    conj = np.frompyfunc(mp.conj, 1, 1)
    n = 5 * N
    A = np.empty((n + 1, n + 1), dtype=object)
    A[:] = mp.mpf(0)
    rhs = np.empty(n + 1, dtype=object)
    rhs[:] = mp.mpf(0)
    eye = np.eye(N, dtype=bool)
    for i in range(5):
        zi = zs[i]
        for j in range(5):
            blk = np.empty((N, N), dtype=object)
            blk[:] = mp.mpf(0)
            imgs = [(r, c) for r in range(3) for c in ((False,) if j < 2 else (False, True))]
            for (r, c) in imgs:
                zj = (conj(zs[j]) if c else zs[j]) * rots[r]
                d = zi[:, None] - zj[None, :]
                if i == j and r == 0 and not c:
                    dc = cq[:, None] - cq[None, :]
                    dc[eye] = mp.mpf(1)
                    d[eye] = dzs[i]
                    blk = blk + W - h * logabs(d / dc)
                else:
                    blk = blk - h * logabs(d)
            A[i * N:(i + 1) * N, j * N:(j + 1) * N] = blk
        rhs[i * N:(i + 1) * N] = -aa[i] * cq
        A[i * N:(i + 1) * N, n] = mp.mpf(-1)
        A[n, i * N:(i + 1) * N] = h * nimg[i]
    t1 = time.time()
    Af = tofloat(A).astype(float)
    lu = lu_factor(Af)
    x = np.array([mp.mpf(v) for v in lu_solve(lu, tofloat(rhs).astype(float))], dtype=object)
    for it in range(6):
        r = rhs - A @ x
        rn = max(abs(v) for v in r)
        print(f'   N={N} refinement {it}: max residual {mp.nstr(rn, 3)}  rho = {mp.nstr(R * mp.exp(-x[-1]), 30)}', flush=True)
        if rn < mp.mpf(10)**(-30):
            break
        dx = lu_solve(lu, tofloat(r).astype(float))
        x = x + np.array([mp.mpf(v) for v in dx], dtype=object)
    rho = R * mp.exp(-x[-1])
    # identity check and tails
    f = x[:-1].reshape(5, N)
    C2 = sum(nimg[i] * h * sum(-aa[i] * cq * f[i]) for i in range(5)) / 2
    tails = []
    for i in range(5):
        cm = (Cm[:, -6:].T @ f[i]) * 2 / N
        tails.append(float(max(abs(v) for v in cm) / max(abs(v) for v in f[i])))
    print(f'N={N:4d} rho={mp.nstr(rho, 30)}  C-C_id={mp.nstr(x[-1] - C2, 3)} tails={["%.1e" % t for t in tails]} '
          f'(assembly {t1 - t0:.0f}s, total {time.time() - t0:.0f}s)', flush=True)
    return rho


if __name__ == '__main__':
    prev = None
    for N in [int(a) for a in sys.argv[1:]]:
        rho = run(N)
        if prev is not None:
            print(f'      diff to previous N: {mp.nstr(rho - prev, 5)}', flush=True)
        prev = rho
