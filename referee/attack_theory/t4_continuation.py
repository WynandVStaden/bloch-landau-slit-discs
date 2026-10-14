"""Step 1 of Section 4.4: is h_i really analytic on E_r with Re h_i = f_i on [-1,1], and |Re h_i| <= B_i on E_r?"""
import numpy as np, sys, time
from scipy.fft import dct
from ubuild import *
TAB = {0:(1.00934,9.4,3196), 3:(1.01733,6.6,1676), 6:(1.00805,8.7,3716), 12:(1.01548,5.7,1872), 18:(1.34403,4.1,88)}
def h_on_path(i, xpath, N=6000, chunk=256):
    """complex h_i along a path of complex x (continuous log, cumulative phase)."""
    y, wt = nodes(N, 0.5)
    y = y[:N//2]; wt = 2*wt   # Gauss-Chebyshev nodes (each cos value appears twice)
    order = np.argsort(y); y = y[order]
    c = COEF[SL[i][4]]; n = np.arange(len(c)); cc = np.zeros(len(c)); cc[1:] = c[1:]*np.pi/n[1:]
    out = -SL[i][1]*xpath + Ch.chebval(xpath, cc)
    kmin, kmax = np.inf, 0
    for j in range(NS):
        ker = (lambda x: Kself(i, x[:,None], y[None,:])) if j == i else (lambda x, j=j: 2*np.sinh(0.5*(wj(i,x)[:,None]-wj(j,y)[None,:])))
        k0 = ker(xpath[:1])[0]
        lg0 = np.log(np.abs(k0)) + 1j*np.unwrap(np.angle(k0))      # continuous in y at first x
        prev, lg = k0, lg0
        res = np.empty(len(xpath), complex)
        for s in range(0, len(xpath), chunk):
            k = ker(xpath[s:s+chunk])
            kmin = min(kmin, np.abs(k).min()); kmax = max(kmax, np.abs(k).max())
            ext = np.vstack([prev[None,:], k])
            inc = np.log(ext[1:]/ext[:-1])                         # principal log of small steps
            assert np.abs(inc.imag).max() < 1.0, "x-steps too coarse"
            lgs = lg[None,:] + np.cumsum(inc, axis=0)
            res[s:s+chunk] = wt*(lgs*psi(j,y)[None,:]).sum(axis=1)
            prev, lg = k[-1], lgs[-1]
        out = out - res
    return out, kmin, kmax
for i in [int(a) for a in sys.argv[1:]] or [18, 0, 6]:
    r, B, M = TAB[i]; t0 = time.time()
    P = 4096
    th = 2*np.pi*np.arange(P)/P
    z = r*np.exp(1j*th); xe = (z+1/z)/2
    H, kmin, kmax = h_on_path(i, np.concatenate([xe, xe[:1]]))
    print("slit %d: closure of h around the ellipse (single-valuedness): |h(end)-h(start)| = %.2e" % (i, abs(H[-1]-H[0])))
    H = H[:-1]
    print("   min|kernel| %.3e max|kernel| %.3e on dE_r x [-1,1];  max|Re h| on dE_r = %.4e  (paper B_i = %.1f);  max|h| = %.3e" % (kmin, kmax, np.abs(H.real).max(), B, np.abs(H).max()))
    F = np.fft.fft(H)/P                  # F[k] = beta_k r^k (k>=0), F[-k] = beta_{-k} r^{-k}
    # real-axis Chebyshev coefficients of f_i
    n = 8192; tt = (np.arange(n)+0.5)*np.pi/n; xr = np.cos(tt)
    f = np.concatenate([f_slit(i, xr[s:s+256], 6000) for s in range(0, n, 256)])
    phi = dct(f, type=2)/n; phi[0] /= 2
    # also h on the real axis from the same routine starting branch: Re must equal f
    ks = np.arange(1, 400)
    beta_out = F[ks]/r**ks; beta_in = F[-ks]*r**ks
    print("   Laurent symmetry beta_k = beta_{-k}: max|diff| = %.2e (max|beta| %.2e)" % (np.abs(beta_out-beta_in)[r**ks < 1e3].max(), np.abs(beta_out).max()))
    d = np.abs(2*beta_out.real - phi[ks])
    print("   max_k<400 |2 Re beta_k (from ellipse) - phi_k (from [-1,1])| = %.2e ; max|phi_k| = %.2e ; sum|phi_k| = %.4e" % (d.max(), np.abs(phi).max(), np.abs(phi).sum()))
    kk = np.arange(1, len(phi))
    with np.errstate(over='ignore'):
        ratio = np.abs(phi[1:])*(r**kk + r**-kk)/(4*B)
    sig = np.abs(phi[1:]) > 3e-16
    if not sig.any(): sig[0] = True
    print("   max |phi_k|(r^k+r^-k)/(4B) over coefficients above noise: %.3e at k=%d ; even-k max |phi_k| = %.1e" % (ratio[sig].max(), kk[sig][np.argmax(ratio[sig])], np.abs(phi[0::2]).max()))
    print("   numerical tail sum_{k>=M}|phi_k| (noise level) = %.2e ; sup|f| on 8192 nodes = %.4e   [%.0f s]" % (np.abs(phi[M:]).sum() if M < len(phi) else 0, np.abs(f).max(), time.time()-t0), flush=True)
