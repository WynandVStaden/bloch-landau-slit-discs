import numpy as np, time
from ubuild import *
rng = np.random.default_rng(7)
# c
c = L - np.pi/4*sum(SL[j][1]*COEF[SL[j][4]][1] for j in range(NS))
print("c =", repr(c), "exp(c) =", repr(np.exp(c)))
# even coefficients zero?
print("even coefs all zero:", all(np.all(COEF[k][0::2]==0) for k in COEF), "len", [len(COEF[k]) for k in COEF])
# (ii) asymptotics
for re in (-1.0, -3.0, -6.0, -12.0, -30.0):
    w = re + 1j*rng.uniform(0, 2*np.pi, 6)
    v = U(w, 4000) + w.real - c
    print("Re w=%6.1f  U+Re w-c: min %.3e max %.3e" % (re, v.min(), v.max()))
# (i) U = 0 on Re w = L (unsymmetric rule shift .3), away from arcs by >=0.004 and anywhere
im = rng.uniform(0, 2*np.pi, 400)
ths = np.array([s[0] for s in SL]) % (2*np.pi)
d = np.min(np.abs(((im[:,None]-ths[None,:]+np.pi) % (2*np.pi))-np.pi), axis=1)
for N in (4000, 16000, 64000):
    v = U(L+1j*im, N, 0.3)
    for lo,hi in ((0.02,9),(0.004,0.02),(0.0,0.004)):
        m = (d>=lo)&(d<hi)
        print("N=%6d dist-to-arc in [%.3f,%.3f): n=%3d max|U(L+it)| = %.3e" % (N, lo, hi, m.sum(), np.abs(v[m]).max() if m.any() else 0))
# oddness under sigma at random points, unsymmetric rule
w = rng.uniform(L-3, L, 200) + 1j*rng.uniform(0, 2*np.pi, 200)
allarc = np.concatenate([wj(j, np.linspace(-1,1,4001)) for j in range(NS)])
def dcyl(w):
    dd = w[:,None]-allarc[None,:]
    imd = (dd.imag+np.pi)%(2*np.pi)-np.pi
    return np.hypot(dd.real, imd).min(axis=1)
dw = dcyl(w); m = dw > 0.02
v1 = U(w[m], 16000, 0.3); v2 = U(2*L-np.conj(w[m]), 16000, 0.3)
print("oddness: max|U(w)+U(sigma w)| =", np.abs(v1+v2).max(), " (n=%d, max|U|=%.3f)" % (m.sum(), np.abs(v1).max()))
# symmetry (d)
v3 = U(w[m]+2j*np.pi/3, 16000, 0.3); v4 = U(np.conj(w[m]), 16000, 0.3); v5 = U(w[m]+2j*np.pi,16000,0.3)
print("rot sym:", np.abs(v1-v3).max(), "conj sym:", np.abs(v1-v4).max(), "2pi i:", np.abs(v1-v5).max())
# positivity inside D (U should be ~ Green function >= -delta)
print("min U at random interior points:", v1.min())
# (iii) mean value test in the w-plane (harmonic there too)
bad = 0; worst = 0
t = 2*np.pi*(np.arange(64)+0.5)/64
for k in np.where(m)[0][:60]:
    rad = 0.5*dw[k]
    ring = w[k] + rad*np.exp(1j*t)
    vals = U(np.concatenate([[w[k]], ring]), 24000, 0.3)
    worst = max(worst, abs(vals[0]-vals[1:].mean()))
print("mean-value defect (60 discs avoiding arcs): max", worst)
# power check: disc crossing a slit
k = 6; cen = wj(k, np.array([-0.5]))[0]
ring = cen + 0.05*np.exp(1j*t)
vals = U(np.concatenate([[cen+0.0], ring]), 48000, 0.3)
print("mean-value defect for disc centred ON slit 6 (should be nonzero):", vals[0]-vals[1:].mean(), " U(centre)=", vals[0])
