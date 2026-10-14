import numpy as np, mpmath as mp
from numpy.polynomial import chebyshev as Ch
rng = np.random.default_rng(3)
def chebcoef(g, n=4096):
    th = (np.arange(n)+0.5)*np.pi/n; v = g(np.cos(th))
    from scipy.fft import dct
    a = dct(v, type=2)/n; a[0]/=2; return a[:n//2]
def test46(h, r, name, n=32768):
    th = np.linspace(0, 2*np.pi, 20001); z = r*np.exp(1j*th); x = (z+1/z)/2
    A = np.abs(np.real(h(x))).max()
    al = chebcoef(lambda t: np.real(h(t+0j)), n)
    k = np.arange(1, len(al))
    with np.errstate(over='ignore'):
        ratio = np.abs(al[1:])*(r**k + r**-k)/(4*A)
    ok = np.abs(al[1:]) > 1e-13*max(1,np.abs(al).max())   # ignore rounding noise
    if not ok.any():
        print("  4.6 %-38s r=%.3f: Re h = 0 on [-1,1], trivial" % (name, r)); return 0.0
    print("  4.6 %-38s r=%.3f A=%.3e max ratio |alpha_k|(r^k+r^-k)/(4A) = %.4f at k=%d" % (name, r, A, ratio[ok].max(), k[ok][np.argmax(ratio[ok])]))
    return ratio[ok].max()
print("Lemma 4.6")
worst = 0
for r in (1.01, 1.1, 1.5, 3.0):
    ax, bx = (r+1/r)/2, (r-1/r)/2
    for eps in (1e-3, 1e-1):
        for nm, x0 in (("pole real axis", ax*(1+eps)), ("pole imag axis", 1j*bx*(1+eps)), ("pole generic", (1+eps)*(ax*np.cos(1)+1j*bx*np.sin(1)))):
            worst = max(worst, test46(lambda x: 1/(x-x0), r, "1/(x-x0) %s eps=%g" % (nm, eps)))
            worst = max(worst, test46(lambda x: 1j/(x-x0), r, "i/(x-x0) %s eps=%g" % (nm, eps)))
            worst = max(worst, test46(lambda x: (0.3-0.8j)/(x-x0)**2, r, "c/(x-x0)^2 %s eps=%g" % (nm, eps)))
            worst = max(worst, test46(lambda x: np.log(1 - x/x0), r, "log(1-x/x0) %s eps=%g" % (nm, eps)))
    worst = max(worst, test46(lambda x: np.exp(3j*x), r, "exp(3ix)"))
    for kk in (1,2,7):
        c = np.zeros(kk+1); c[kk]=1
        worst = max(worst, test46(lambda x: Ch.chebval(x,c), r, "T_%d" % kk))
        worst = max(worst, test46(lambda x: (1+2j)*Ch.chebval(x,c), r, "(1+2i)T_%d" % kk))
        worst = max(worst, test46(lambda x: 1j*Ch.chebval(x,c), r, "i T_%d (real part 0 on [-1,1])" % kk)) if False else worst
# adversarial random search: complex polynomials in zeta-symmetric form; maximise ratio for a given k
print("  random search over complex polynomials (deg<=12), 20000 trials per r")
for r in (1.01, 1.3, 2.5):
    th = np.linspace(0, 2*np.pi, 2001); z = r*np.exp(1j*th); x = (z+1/z)/2
    best = 0
    V = np.array([Ch.chebval(x, np.eye(13)[q]) for q in range(13)])   # T_q on ellipse
    for trial in range(20000):
        cf = (rng.standard_normal(13)+1j*rng.standard_normal(13))*rng.random(13)**3
        A = np.abs(np.real(cf@V)).max()
        k = np.arange(1,13)
        best = max(best, (np.abs(cf.real[1:])*(r**k+r**-k)/(4*A)).max())
    # near-extremal: Re H(r e^{it}) ~ square wave sign(cos kt)
    print("    r=%.2f best ratio found %.4f" % (r, best)); worst = max(worst, best)
print("  WORST ratio overall (must be <= 1):", worst)

print("Lemma 4.5")
mp.mp.dps = 25
for n in (1,2,3,7,20):
    for x in (-1, -0.3, 0, 0.77, 1):
        ph = mp.acos(x)
        g = lambda t: -mp.log(abs(mp.mpf(x)-mp.cos(t)))*mp.cos(n*t)
        val = mp.quad(g, [0, ph, mp.pi]) if 0<ph<mp.pi else mp.quad(g, [0, mp.pi])
        print("  n=%2d x=%5.2f  integral=%+.14f  (pi/n)T_n=%+.14f" % (n, x, val, mp.pi/n*mp.cos(n*ph)))
g0 = lambda t: -mp.log(abs(mp.mpf('0.3')-mp.cos(t)))
print("  n=0: integral =", mp.quad(g0,[0,mp.acos(0.3),mp.pi]), " pi log 2 =", mp.pi*mp.log(2), "(n=0 would need this term; not needed since c_{j,0}=0)")

print("Lemma 4.7 (quadrature)")
worst = 0
for trial in range(300):
    r = 1 + 10**rng.uniform(-2, 0.3)
    ax = (r+1/r)/2
    x0 = ax*(1+10**rng.uniform(-3,-0.5))*np.exp(1j*rng.uniform(0,0.2))  # pole just outside
    F = lambda y: np.real(1/(y-x0)) + np.real(np.exp(2j*y))
    al = chebcoef(F, 32768); k = np.arange(1,len(al))
    with np.errstate(over='ignore', invalid='ignore'):
        C = np.nanmax(np.where(np.abs(al[1:])>1e-13, np.abs(al[1:])*r**k, 0))
    d = int(rng.integers(0, 30)); pc = rng.standard_normal(d+1); 
    N = int(rng.integers((d+1)//2+1, 80))
    kap = 2*N-1-d
    if kap < 0: continue
    psin = np.abs(Ch.chebval(np.linspace(-1,1,20001), pc)).max()
    def rule(n):
        y = np.cos((2*np.arange(1,n+1)-1)*np.pi/(2*n)); return np.pi/n*(F(y)*Ch.chebval(y,pc)).sum()
    err = abs(rule(20000)-rule(N)); bd = 2*np.pi*C*psin/(r-1)*r**(-kap)
    if bd > 1e-11: worst = max(worst, err/bd)
print("  max err/bound over random trials (must be <= 1):", worst)

print("Lemma 4.8 (interpolation)")
# aliasing identity
bad = 0
for M in (2, 5, 8, 31):
    xl = np.cos((2*np.arange(M)+1)*np.pi/(2*M))
    for k in range(0, 6*M+3):
        m = (k+M)//(2*M); kp = abs(k-2*M*m)   # k = 2Mm +- kp, 0<=kp<=M
        tk = np.cos(k*np.arccos(xl)); ref = (-1)**m*np.cos(kp*np.arccos(xl))
        if kp == M: ref = 0*ref
        bad = max(bad, np.abs(tk-ref).max())
print("  aliasing identity max defect:", bad)
worst = 0; worst2 = 0
for trial in range(300):
    r = 1 + 10**rng.uniform(-2, 0.3); ax = (r+1/r)/2
    x0 = ax*(1+10**rng.uniform(-3,-0.5))*np.exp(1j*rng.uniform(0,3.1))
    f = lambda y: np.real((0.4+1j)/(y-x0))
    al = chebcoef(f, 32768); k = np.arange(1,len(al))
    with np.errstate(over='ignore', invalid='ignore'):
        C = np.nanmax(np.where(np.abs(al[1:])>1e-13, np.abs(al[1:])*r**k, 0))
    M = int(rng.integers(1, 120))
    xl = np.cos((2*np.arange(M)+1)*np.pi/(2*M))
    q = np.arange(M); gam = (2/M)*np.cos(np.outer(q, np.arccos(xl)))@f(xl); gam[0]/=2
    xx = np.cos(np.linspace(0,np.pi,40001))
    e = np.abs(f(xx)-Ch.chebval(xx,gam)).max(); bd = 2*C*r**(1-M)/(r-1)
    tail = 2*np.abs(al[M:]).sum()
    if bd > 1e-10: worst = max(worst, e/bd)
    if tail > 1e-10: worst2 = max(worst2, e/tail)
    assert np.abs(f(xx)).max() <= np.abs(gam).sum() + bd + 1e-12
print("  max ||f-p|| / (2C r^(1-M)/(r-1)) =", worst, " ; max ||f-p|| / (2 sum_{k>=M}|phi_k|) =", worst2, "(both must be <= 1)")
# sharpness: f = T_M -> interpolant 0, error 1 ; f = T_{M+1}: interpolant -T_{M-1}, error up to 2
M=10; xl=np.cos((2*np.arange(M)+1)*np.pi/(2*M)); f=lambda y: np.cos((M+1)*np.arccos(y))
gam=(2/M)*np.cos(np.outer(np.arange(M),np.arccos(xl)))@f(xl); gam[0]/=2
print("  T_{M+1}: interpolant coefficients nonzero:", {int(i): round(float(g),6) for i,g in enumerate(gam) if abs(g)>1e-9}, " ||f-p|| =", np.abs(f(xx)-Ch.chebval(xx,gam)).max())
