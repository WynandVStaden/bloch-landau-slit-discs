"""Independent construction of U (referee, attack_theory). numpy only; own code."""
import json, numpy as np
from numpy.polynomial import chebyshev as Ch
R = 5.001117628756668
L = np.log(R)
CLASSES = [(0.0, 1.0, 0.0, 0.0), (np.pi/3, 1.9999999974092333, 0.0, 0.0),
 (0.029930188727091613, 2.789342409403121, 0.9026921410641155, -0.3301459494660587),
 (1.0157094972537686, 3.473381248030611, -1.2656124700253997, 0.9794269477754911),
 (0.5757776442359243, 4.342874184977691, 0.03225521024585182, 2.9360081049550737)]
SL = []  # (theta, a, b1, b2, cls)
for c,(th,t,b1,b2) in enumerate(CLASSES):
    for r in range(3):
        for e in ((1,) if c < 2 else (1,-1)):
            SL.append((2*np.pi*r/3+e*th, np.log(R/t), e*b1, e*b2, c))
NS = len(SL)
with open(r"C:\Users\wvs00\claude_chat\bloch_record\certify\density_G_n64.json") as fh:
    COEF = {int(k): np.array(v) for k,v in json.load(fh)["coefficients"].items()}

def wj(j, x):
    th,a,b1,b2,_ = SL[j]; s = a*x
    return L + s + 1j*(th + b1*s**2 + b2*s**4)
def psi(j, y):
    return Ch.chebval(y, COEF[SL[j][4]])
def Phi(w):
    w = np.where(w.real < 0, -w, w)
    return -0.5*w.real - np.log(np.abs(1 - np.exp(-w)))
def nodes(N, shift=0.5):
    """periodic trapezoid in theta over [0,2pi): int g dmu ~ (pi/N) sum g(cos th_l). shift=.5 with N even = Gauss-Chebyshev doubled"""
    th = 2*np.pi*(np.arange(N)+shift)/N
    return np.cos(th), np.pi/N
def U(w, N=8000, shift=0.3):
    w = np.atleast_1d(np.asarray(w, complex))
    y, wt = nodes(N, shift)
    out = L - w.real
    for j in range(NS):
        out = out + wt*(Phi(w[:,None]-wj(j,y)[None,:])*psi(j,y)[None,:]).sum(axis=1)
    return out
def shc(t):
    t = np.asarray(t, complex); small = np.abs(t) < 1e-3
    ts = np.where(small, 1.0, t)
    return np.where(small, 1 + t*t/6 + t**4/120, np.sinh(ts)/ts)
def Kself(i, x, y):
    th,a,b1,b2,_ = SL[i]
    dd = a + 1j*(b1*a**2*(x+y) + b2*a**4*(x**3+x*x*y+x*y*y+y**3))
    return dd*shc(0.5*(x-y)*dd)
def f_slit(i, x, N=6000, parts=False):
    """formula (7) of the paper, x real array (or complex: returns h-like real part using log|.|)"""
    x = np.atleast_1d(x)
    y, wt = nodes(N, 0.5)
    c = COEF[SL[i][4]]
    n = np.arange(len(c)); cc = np.zeros(len(c)); cc[1:] = c[1:]*np.pi/n[1:]
    out = -SL[i][1]*x + Ch.chebval(x, cc)
    out = np.real(out) if np.iscomplexobj(out) else out
    for j in range(NS):
        if j == i:
            k = Kself(i, x[:,None], y[None,:])
        else:
            k = 2*np.sinh(0.5*(wj(i,x)[:,None]-wj(j,y)[None,:]))
        out = out - wt*(np.log(np.abs(k))*psi(j,y)[None,:]).sum(axis=1)
    return out
