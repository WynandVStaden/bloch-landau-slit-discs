import numpy as np
from ubuild import *
from flint import acb, arb
print("flint sinc(1) =", acb(1).sinc(), " sin(1) =", np.sin(1.0))
rng = np.random.default_rng(1)
# divided-difference identity for complex x, y
worst = 0
for i in (0, 6, 13, 18, 23):
    x = rng.uniform(-1.2,1.2,200)+1j*rng.uniform(-.4,.4,200); y = rng.uniform(-1.2,1.2,200)+1j*rng.uniform(-.4,.4,200)
    lhs = Kself(i,x,y)*(x-y); rhs = 2*np.sinh(0.5*(wj(i,x)-wj(i,y)))
    worst = max(worst, np.abs(lhs-rhs).max())
print("max |K_i (x-y) - 2 sinh((w_i(x)-w_i(y))/2)| for complex x,y:", worst)
# sigma identity
x = rng.uniform(-1,1,50)
print("max |sigma(w_j(x)) - w_j(-x)|:", max(np.abs(2*L-np.conj(wj(j,x))-wj(j,-x)).max() for j in range(NS)))
# two-sided approach to a slit and approach to a tip: U -> f ~ 0 continuously
def Uhi(w, N):
    return U(w, N, 0.5)
for i, x0 in ((0, -0.5), (6, -0.5), (6, -0.97), (12, -0.05)):
    th,a,b1,b2,_ = SL[i]
    w0 = wj(i, np.array([x0]))[0]
    s = a*x0; tang = a*(1 + 1j*(2*b1*s + 4*b2*s**3)); nrm = 1j*tang/abs(tang)
    for eps in (1e-2, 1e-3, 1e-4):
        N = int(min(2_000_000, 40/eps*4))
        vals = Uhi(np.array([w0+eps*nrm, w0-eps*nrm]), N)
        print("slit %2d x=%.2f eps=%.0e (N=%d): U(+)=%.4e U(-)=%.4e ; U/eps = %.4f, %.4f" % (i, x0, eps, N, vals[0], vals[1], vals[0]/eps, vals[1]/eps))
for i in (0, 6, 18):
    tip = wj(i, np.array([-1.0]))[0]
    for eps in (1e-2, 1e-3, 1e-4):
        N = int(min(2_000_000, 40/eps*4))
        ang = np.array([np.pi, np.pi*0.6, -np.pi*0.6, 0.3, -0.3])  # pi = straight ahead of the tip (roughly), 0.3 = alongside the slit
        vals = Uhi(tip + eps*np.exp(1j*ang), N)
        print("tip of slit %2d, eps=%.0e: U/sqrt(eps) at 5 directions:" % (i, eps), np.round(vals/np.sqrt(eps), 4))
