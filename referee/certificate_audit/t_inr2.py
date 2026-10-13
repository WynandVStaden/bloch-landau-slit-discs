import sys, io, contextlib, math
from fractions import Fraction
import numpy as np
sys.path.insert(0, "../../certify")
import inradius_certificate as ic
G = ic.DESIGNS["G"]
def tips(des):
    R = float(des[0]); out = []
    for r, sgn, alpha, tip, b1, b2 in ic.slit_table(des[1]):
        a = math.log(R / float(tip)); th = ic.float_angle(r, sgn, alpha) + float(b1) * a**2 + float(b2) * a**4
        out.append(float(tip) * np.exp(1j * th))
    return out
def circumradius(p, q, r):
    a, b, c = abs(p - q), abs(q - r), abs(r - p); s = (a + b + c) / 2
    return a * b * c / (4 * math.sqrt(s * (s - a) * (s - b) * (s - c)))
for db in ("0.9026921410641155", "0.9026922410641155", "0.9026920410641155"):
    des = (G[0], list(G[1][:2]) + [(G[1][2][0], G[1][2][1], db, G[1][2][3])] + list(G[1][3:]))
    t = tips(des); cr = circumradius(t[6], t[12], t[18])
    sl = ic.slit_table(des[1]); br = Fraction(des[0])
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            c2, p2 = ic.search(float(br), sl, 1e-9, max_level=36)
            ic.no_overlap(c2); res = ic.verify(br, sl, Fraction(1e-9), c2, p2)
        v = "CERTIFIED <= 1+1e-9"
    except (SystemExit, AssertionError) as ex:
        v = "not certified"
    print(f"class-2 b1 = {db}: circumradius of tips 6C,12D,18E = 1 + {cr-1:.3e}  -> {v}", "" if (cr <= 1 + 1e-9) or v == "not certified" else "*** FALSE CERTIFICATE ***")
