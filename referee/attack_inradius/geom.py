# Independent geometry of the slit disc (referee code; written before reading the repository code)
import numpy as np, math
from mpmath import mp, mpf, mpc
mp.dps = 50
R_S = "5.001117628756668"
CLS_S = [("0", "1", "0", "0"),
         ("pi/3", "1.9999999974092333", "0", "0"),
         ("0.029930188727091613", "2.789342409403121", "0.9026921410641155", "-0.3301459494660587"),
         ("1.0157094972537686", "3.473381248030611", "-1.2656124700253997", "0.9794269477754911"),
         ("0.5757776442359243", "4.342874184977691", "0.03225521024585182", "2.9360081049550737")]
import os, json
if os.environ.get('GEOM_CLS'): CLS_S = [tuple(c) for c in json.loads(os.environ['GEOM_CLS'])]  # referee mutation experiments only
if os.environ.get('GEOM_R'): R_S = os.environ['GEOM_R']
def slits_mp():
    R = mpf(R_S); out = []
    for ci, (th, t, b1, b2) in enumerate(CLS_S):
        th = mp.pi/3 if th == "pi/3" else mpf(th)
        t, b1, b2 = mpf(t), mpf(b1), mpf(b2)
        eps = [1] if ci < 2 else [1, -1]
        for r in range(3):
            for e in eps:
                out.append(dict(cls=ci, r=r, eps=e, th=2*mp.pi*r/3 + e*th, t=t, b1=e*b1, b2=e*b2, a=mp.log(R/t)))
    return R, out
R_MP, SL_MP = slits_mp()
R_F = float(R_MP)
SL_F = [dict(cls=s['cls'], r=s['r'], eps=s['eps'], th=float(s['th']), t=float(s['t']), b1=float(s['b1']), b2=float(s['b2']), a=float(s['a'])) for s in SL_MP]

def gam_f(sl, s):
    return R_F*np.exp(s + 1j*(sl['th'] + sl['b1']*s**2 + sl['b2']*s**4))
def gam(sl, s, d=0):
    g = R_MP*mp.exp(mpc(s, sl['th'] + sl['b1']*s**2 + sl['b2']*s**4))
    if d == 0: return g
    p1 = mpc(1, 2*sl['b1']*s + 4*sl['b2']*s**3)
    if d == 1: return g, g*p1
    return g, g*p1, g*(p1*p1 + mpc(0, 2*sl['b1'] + 12*sl['b2']*s**2))

NS = 3001
SAMP = [gam_f(sl, np.linspace(-sl['a'], 0, NS)) for sl in SL_F]
SGRID = [np.linspace(-sl['a'], 0, NS) for sl in SL_F]

def near_f(k, p):
    """double: nearest point on slit k to p: returns (dist, s)"""
    sl = SL_F[k]; d = np.abs(SAMP[k] - p); j = int(np.argmin(d)); s = SGRID[k][j]
    lo, hi = -sl['a'], 0.0
    # all local minima among samples? take global sample min then Newton
    for _ in range(40):
        g = gam_f(sl, s); p1 = 1 + 1j*(2*sl['b1']*s + 4*sl['b2']*s**3)
        g1 = g*p1; g2 = g*(p1*p1 + 1j*(2*sl['b1'] + 12*sl['b2']*s*s))
        f = ((g - p)*np.conj(g1)).real; fp = abs(g1)**2 + ((g - p)*np.conj(g2)).real
        if fp <= 0: break
        sn = min(hi, max(lo, s - f/fp))
        if abs(sn - s) < 1e-16: s = sn; break
        s = sn
    dd = abs(gam_f(sl, s) - p)
    if dd > d[j]: dd, s = d[j], SGRID[k][j]
    return dd, s
def dist_f(p):
    best = (R_F - abs(p), -1, 0.0)
    for k in range(len(SL_F)):
        d, s = near_f(k, p)
        if d < best[0]: best = (d, k, s)
    return best
def alld_f(p):
    L = [(R_F - abs(p), -1, 0.0)]
    for k in range(len(SL_F)):
        d, s = near_f(k, p); L.append((d, k, s))
    return sorted(L)

def near_mp(k, c, s0=None):
    """mp: nearest point on slit k to c (global via double sampling start), returns (dist, s, kind)"""
    sl = SL_MP[k]
    cands = []
    d = np.abs(SAMP[k] - complex(c))
    # all sample local minima
    idx = [j for j in range(NS) if (j == 0 or d[j] <= d[j-1]) and (j == NS-1 or d[j] <= d[j+1])]
    idx = sorted(idx, key=lambda j: d[j])[:3]
    for j in idx:
        s = mpf(float(SGRID[k][j])); lo, hi = -sl['a'], mpf(0)
        for _ in range(100):
            g, g1, g2 = gam(sl, s, 2)
            f = mp.re((g - c)*mp.conj(g1)); fp = abs(g1)**2 + mp.re((g - c)*mp.conj(g2))
            if fp <= 0: break
            sn = min(hi, max(lo, s - f/fp))
            if abs(sn - s) < mpf(10)**(-mp.dps + 5): s = sn; break
            s = sn
        kind = 'tip' if s == lo else ('root' if s == hi else 'int')
        cands.append((abs(gam(sl, s) - c), s, kind))
    return min(cands, key=lambda t: t[0])
def dist_mp(c):
    L = [(R_MP - abs(c), -1, mpf(0), 'circle')]
    for k in range(len(SL_MP)):
        d, s, kind = near_mp(k, c); L.append((d, k, s, kind))
    return sorted(L, key=lambda t: t[0])
