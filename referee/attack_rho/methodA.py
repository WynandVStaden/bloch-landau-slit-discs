"""Method A: global spectral (Chebyshev / cosine) single-layer method on the DOUBLED arcs.

Each slit z(s)=R exp(s+i v(s)), s in [-a,0], v even in s, continues analytically through |z|=R to its
mirror image in the circle (s in [0,a]).  The Green function G of D with pole 0, continued oddly by reflection,
is   G(z) = log(R/|z|) - (U(z) - C),   U(z) = sum_j int_0^pi -log|z - z_j(a_j cos p)| f_j(p) dp
with zero total charge.  Unknowns: f_j at midpoint nodes, and the constant C.   log rho = log R - C.
The free-space log kernel is split as -log|cos p - cos q| (integrated exactly against cosines) + analytic.
"""
import sys
import numpy as np
from geom import full_problem, slit_z, slit_dz, apply_image


def nodes(N):
    return (np.arange(N) + 0.5) * np.pi / N


def logcos_weights(p_t, N):
    """W[i,k]: int_0^pi -log|cos p_t[i] - cos q| f(q) dq ~ sum_k W[i,k] f(q_k), exact for cos(mq), m<N."""
    q = nodes(N)
    m = np.arange(1, N)
    Ct = np.cos(np.outer(p_t, m)) * (2.0 / m)
    Cq = np.cos(np.outer(q, m))
    return (np.pi / N) * (np.log(2.0) + Ct @ Cq.T)


def assemble(R, classes, N):
    nc = len(classes)
    q = nodes(N)
    cq = np.cos(q)
    n = nc * N
    A = np.zeros((n + 1, n + 1))
    rhs = np.zeros(n + 1)
    W = logcos_weights(q, N)
    h = np.pi / N
    zs = [slit_z(R, c, np.log(R / c['t']) * cq) for c in classes]
    for i, ci in enumerate(classes):
        ai = np.log(R / ci['t'])
        zi = zs[i]
        for j, cj in enumerate(classes):
            blk = np.zeros((N, N))
            for img in cj['images']:
                zj = apply_image(zs[j], img)
                if i == j and img == cj['images'][0]:
                    d = zi[:, None] - zj[None, :]
                    dc = cq[:, None] - cq[None, :]
                    np.fill_diagonal(dc, 1.0)
                    M = -np.log(np.abs(d / dc), where=~np.eye(N, dtype=bool), out=np.zeros((N, N)))
                    dz = slit_dz(R, ci, ai * cq) * ai
                    np.fill_diagonal(M, -np.log(np.abs(dz)))
                    blk += W + h * M
                else:
                    blk += -h * np.log(np.abs(zi[:, None] - zj[None, :]))
            A[i * N:(i + 1) * N, j * N:(j + 1) * N] = blk
        rhs[i * N:(i + 1) * N] = -ai * cq          # log(R/|z|) = -s
        A[i * N:(i + 1) * N, n] = -1.0
        A[n, i * N:(i + 1) * N] = h * len(ci['images'])
    return A, rhs


def solve(R, classes, N):
    A, rhs = assemble(R, classes, N)
    x = np.linalg.solve(A, rhs)
    nc = len(classes)
    f = x[:-1].reshape(nc, N)
    C = x[-1]
    return f, C, A, rhs


def postprocess(R, classes, N, f, C):
    q = nodes(N)
    h = np.pi / N
    out = {}
    out['rho'] = R * np.exp(-C)
    # C as half of sum_j int (-a cos q) f dq  (independent identity)
    C2 = 0.0
    hm = []
    m = np.arange(1, N)
    for i, c in enumerate(classes):
        a = np.log(R / c['t'])
        C2 += 0.5 * len(c['images']) * h * np.sum(-a * np.cos(q) * f[i])
        # harmonic measure of ONE slit (both sides) = int_{pi/2}^{pi} f
        cm = (2.0 / N) * (np.cos(np.outer(m, q)) @ f[i])
        c0 = np.sum(f[i]) / N
        hm.append(c0 * np.pi / 2 + np.sum(cm * (-np.sin(m * np.pi / 2) / m)))
    out['C'] = C
    out['C_id'] = C2
    out['hm_per_slit'] = np.array(hm)
    out['hm_class'] = np.array([hm[i] * len(c['images']) for i, c in enumerate(classes)])
    out['antisym'] = max(np.max(np.abs(f[i] + f[i][::-1])) / np.max(np.abs(f[i])) for i in range(len(classes)))
    out['tail'] = [float(np.abs((2.0 / N) * (np.cos(np.outer(np.arange(N - 4, N), q)) @ f[i])).max()
                         / np.abs(f[i]).max()) for i in range(len(classes))]
    return out


def eval_G(R, classes, N, f, C, zt):
    """Green function approximation at arbitrary points zt NOT on slits (plain quadrature)."""
    q = nodes(N)
    h = np.pi / N
    U = np.zeros(zt.shape)
    for i, c in enumerate(classes):
        a = np.log(R / c['t'])
        zb = slit_z(R, c, a * np.cos(q))
        for img in c['images']:
            zj = apply_image(zb, img)
            U += -h * (np.log(np.abs(zt[:, None] - zj[None, :])) @ f[i])
    return np.log(R / np.abs(zt)) - (U - C)


def residual_on_slits(R, classes, N, f, C, M=None, K=None):
    """Residual of the boundary condition G=0 on the doubled arcs at off-node parameters, using the
    trigonometric interpolant of f (degree < N) and a K-point oversampled rule for the analytic kernel parts.
    By the maximum principle max|residual| bounds |log rho_N - log rho| (up to the oversampling error)."""
    if K is None:
        K = 4 * N
    if M is None:
        M = 3 * N + 7
    qN = nodes(N)
    qK = nodes(K)
    # interpolate f to K nodes by cosine series
    m = np.arange(N)
    coef = (2.0 / N) * (np.cos(np.outer(m, qN)) @ f.T)        # (N, nc)
    coef[0] *= 0.5
    fK = (np.cos(np.outer(qK, m)) @ coef).T                    # (nc, K)
    pt = (np.arange(M) + 0.37) * np.pi / M
    hK = np.pi / K
    W = logcos_weights(pt, K)
    res = []
    cK = np.cos(qK)
    ct = np.cos(pt)
    for i, ci in enumerate(classes):
        ai = np.log(R / ci['t'])
        zi = slit_z(R, ci, ai * ct)
        U = np.zeros(M)
        for j, cj in enumerate(classes):
            aj = np.log(R / cj['t'])
            zb = slit_z(R, cj, aj * cK)
            for img in cj['images']:
                zj = apply_image(zb, img)
                if i == j and img == cj['images'][0]:
                    Mm = -np.log(np.abs((zi[:, None] - zj[None, :]) / (ct[:, None] - cK[None, :])))
                    U += (W + hK * Mm) @ fK[j]
                else:
                    U += -hK * (np.log(np.abs(zi[:, None] - zj[None, :])) @ fK[j])
        res.append(U - C + ai * ct)
    return np.array(res)


if __name__ == '__main__':
    R, classes = full_problem()
    Ns = [int(a) for a in sys.argv[1:]] or [16, 24, 32, 48, 64, 96, 128, 192, 256, 384]
    prev = None
    for N in Ns:
        f, C, A, rhs = solve(R, classes, N)
        o = postprocess(R, classes, N, f, C)
        d = '' if prev is None else f'{o["rho"] - prev:+.3e}'
        print(f'N={N:4d} rho={o["rho"]:.15f} diff={d:>11s} C-C_id={o["C"] - o["C_id"]:+.2e} antisym={o["antisym"]:.1e}'
              f' cond={np.linalg.cond(A):.2e}', flush=True)
        prev = o['rho']
    print('harmonic measure per class (all slits of class, both sides):', o['hm_class'], 'circle:', 1 - o['hm_class'].sum())
    print('cos-coefficient tails:', o['tail'])
