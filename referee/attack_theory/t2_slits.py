import numpy as np, time, mpmath as mp
from ubuild import *
t0=time.time()
# A. independent check of formula (7) incl. Lemma 4.5: raw potential by mpmath tanh-sinh (self term), cross terms by GL in theta
mp.mp.dps = 30
def raw_self(i, x):
    th,a,b1,b2,cl = SL[i]; c = [mp.mpf(float(v)) for v in COEF[cl]]
    a=mp.mpf(a); b1=mp.mpf(b1); b2=mp.mpf(b2)
    def wv(t):
        s=a*t; return s + 1j*(b1*s**2+b2*s**4)
    wx = wv(mp.mpf(x)); ph = mp.acos(mp.mpf(x))
    def g(t):
        y = mp.cos(t); p = mp.fsum(c[n]*mp.cos(n*t) for n in range(1,64,2))
        return -mp.log(abs(2*mp.sinh((wx-wv(y))/2)))*p
    return mp.quad(g, [0, ph, mp.pi]) if 0 < ph < mp.pi else mp.quad(g, [0, mp.pi])
def split_self(i, x, N=6000):
    y, wt = nodes(N, 0.5); c = COEF[SL[i][4]]
    n = np.arange(len(c)); cc = np.zeros(len(c)); cc[1:] = c[1:]*np.pi/n[1:]
    return Ch.chebval(x, cc) - wt*(np.log(np.abs(Kself(i, x, y)))*psi(i,y)).sum()
print("A. self term: mpmath raw integral vs split formula (Lemma 4.5 + K_i)")
for i in (0, 6, 7, 12, 18):
    for x in (-1.0, -0.999999, -0.73, -0.2, -1e-7, 0.0, 0.4, 1.0):
        r = raw_self(i, x); s = split_self(i, x)
        print("  slit %2d x=%10.7f raw=%+.15f split=%+.15f diff=%.2e" % (i, x, float(r), s, float(r)-s))
print("time", time.time()-t0)
# B. dense scan of f_i on class representatives
def xs_dense(n):
    t = (np.arange(n)+0.5)/n*np.pi/2
    x = -np.sin(t)**2    # clusters at 0 and -1
    ends = np.concatenate([-np.logspace(-15,-3,40), -1+np.logspace(-15,-3,40), [0.0,-1.0]])
    return np.sort(np.concatenate([x, ends]))
x = xs_dense(1500)
print("B. dense scan, N=6000 (and N=9000 convergence check)")
for i in (0,3,6,12,18):
    f = np.concatenate([f_slit(i, x[k:k+200], 6000) for k in range(0,len(x),200)])
    f2 = f_slit(i, x[::37], 9000)
    k = np.argmax(np.abs(f))
    print("  slit %2d class %d: sup|f| = %.4e at x=%.6f ; f(-1)=%.3e f(0)=%.3e ; near tip max %.3e, near circle max %.3e; N-conv %.1e" % (
        i, SL[i][4], np.abs(f).max(), x[k], f[0], f[-1], np.abs(f[x<-0.999]).max(), np.abs(f[x>-0.001]).max(), np.abs(f[::37]-f2).max()), flush=True)
print("time", time.time()-t0)
# C. all 24 slits, coarser
xc = xs_dense(60)
print("C. all 24 slits (symmetry claim): sup|f_j| and difference from class representative")
rep = {}
for j in range(NS):
    f = f_slit(j, xc, 6000)
    cl = SL[j][4]
    if cl not in rep: rep[cl] = f
    print("  slit %2d class %d sup|f|=%.4e  max|f_j - f_rep|=%.2e" % (j, cl, np.abs(f).max(), np.abs(f-rep[cl]).max()), flush=True)
# D. oddness of f_i
xo = np.linspace(-1,1,41)
for i in (0,6,12,18):
    f = f_slit(i, xo, 6000); print("  oddness slit", i, np.abs(f+f[::-1]).max())
print("time", time.time()-t0)
