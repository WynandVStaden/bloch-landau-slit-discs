# Independent verifier of a stored inradius certificate (referee code).
# Exact rationals (fractions / Python ints) for cells; mpmath interval arithmetic (iv, 200 bits) for witnesses.
# Does NOT import the authors' program.  Usage: python verify_leaves.py leaves.npz [eps as decimal text]
import sys, numpy as np
from fractions import Fraction as Fr
from mpmath import iv
iv.prec = 200
R_S = "5.001117628756668"
CLS = [("0", "1", "0", "0"), ("pi/3", "1.9999999974092333", "0", "0"),
       ("0.029930188727091613", "2.789342409403121", "0.9026921410641155", "-0.3301459494660587"),
       ("1.0157094972537686", "3.473381248030611", "-1.2656124700253997", "0.9794269477754911"),
       ("0.5757776442359243", "4.342874184977691", "0.03225521024585182", "2.9360081049550737")]
def I(q):
    q = Fr(q); return iv.mpf(q.numerator)/iv.mpf(q.denominator)
def table(cls=CLS, drop=None, shorten=None):
    T = []
    for ci, (al, t, b1, b2) in enumerate(cls):
        for r in range(3):
            for e in ((1,) if ci < 2 else (1, -1)):
                T.append((ci, r, e, al, Fr(t), e*Fr(b1), e*Fr(b2)))
    return T
def main(fn, eps_text="1e-9", T=None, quiet=False):
    T = T or table()
    R = Fr(R_S); lam = 1 + Fr(eps_text); assert lam > 0
    lam2 = I(lam)*I(lam); Ri = I(R)
    z = np.load(fn, allow_pickle=False)
    cells, pars = z["cells"], z["pars"]
    assert cells.ndim == 2 and cells.shape[1] == 6 and len(pars) == len(cells) and cells.dtype.kind == 'i'
    cl = [tuple(int(v) for v in row) for row in cells.tolist()]; pl = [float(p) for p in pars.tolist()]
    leaf = {}
    for n, (lv, ix, iy, code, kind, j) in enumerate(cl):
        assert lv >= 6 and 0 <= ix < 2**lv and 0 <= iy < 2**lv, "bad cell"
        assert (lv, ix, iy) not in leaf, "duplicate cell"
        leaf[(lv, ix, iy)] = n
    maxlv = max(c[0] for c in cl)
    # ---- tiling of [0,8]^2 by explicit recursive descent: each leaf must be reached exactly once
    visited = 0
    stack = [(6, i, j) for i in range(64) for j in range(64)]
    while stack:
        lv, i, j = stack.pop()
        if (lv, i, j) in leaf: visited += 1; continue
        assert lv < maxlv, "uncovered region: cell (%d,%d,%d) is neither a leaf nor subdivided" % (lv, i, j)
        stack += [(lv + 1, 2*i + di, 2*j + dj) for di in (0, 1) for dj in (0, 1)]
    assert visited == len(cl), "leaves overlap (some leaf lies under another leaf)"
    # ---- each leaf
    cache = {}; cnt = {'a': 0, 'b': 0, 'c': 0}; worst = None
    def wit(kind, j, par):
        if kind == 2:
            th = I(Fr(par)); return Ri*iv.cos(th), Ri*iv.sin(th)
        ci, r, e, al, t, b1, b2 = T[j] if 0 <= j < len(T) else (_ for _ in ()).throw(AssertionError("bad slit index"))
        a = iv.log(Ri/I(t))
        if kind == 1: s = -a
        else:
            assert kind == 0, "unknown witness kind"
            s = I(Fr(par)); assert par <= 0 and ((s + a).a > 0) is True, "parameter not in (-a,0]"
        al_i = iv.pi/3 if al == "pi/3" else I(al)
        ph = 2*iv.pi*r/3 + e*al_i + I(b1)*s**2 + I(b2)*s**4
        rad = Ri*iv.exp(s)
        return rad*iv.cos(ph), rad*iv.sin(ph)
    for (lv, ix, iy, code, kind, j), par in zip(cl, pl):
        h = Fr(8, 2**lv); x0, y0 = ix*h, iy*h; x1, y1 = x0 + h, y0 + h
        if code == 0:
            assert x0*x0 + y0*y0 >= R*R, "(a) fails"; cnt['a'] += 1
        elif code == 1:
            assert y0*y0 > 3*x1*x1, "(b) fails"; cnt['b'] += 1
        else:
            assert code == 2, "unknown code"
            key = (kind, j, par)
            if key not in cache: cache[key] = wit(kind, j, par)
            qx, qy = cache[key]
            for cx in (x0, x1):
                for cy in (y0, y1):
                    d2 = (I(cx) - qx)**2 + (I(cy) - qy)**2
                    ma = (lam2 - d2).a  # lower end of the interval lam^2 - d^2 (degenerate interval)
                    assert (ma >= 0) is True, "(c) fails: corner too far from witness, cell %s" % ((lv, ix, iy),)
                    m = float(ma)
                    if worst is None or m < worst[0]: worst = (m, lv, ix, iy)
            cnt['c'] += 1
    if not quiet:
        print("VERIFIED %s: %d leaves tile [0,8]^2; (a) %d, (b) %d, (c) %d; %d distinct witnesses; max level %d; smallest margin lam^2-d^2 = %.3g at %s; claim: inradius <= 1 + %s"
              % (fn, len(cl), cnt['a'], cnt['b'], cnt['c'], len(cache), maxlv, worst[0], worst[1:], eps_text))
    return cnt
if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "1e-9")
