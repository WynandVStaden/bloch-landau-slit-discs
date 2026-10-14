import sys, os, random, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "certify"))
import conformal_certificate as cc
from flint import arb, acb, ctx
from fractions import Fraction
ctx.prec = 100
box = [(-1.0, 1.0), (-0.5, 0.5), (-1.0, 1.0), (0.0, 0.0)]
sc = [1, 1, 1, 0]
for z0 in (0.3 + 0.2j, 1.0 + 0.5j, 1.0 + 0.5j + 1e-9, 1.0 + 0.5000001j, 1.001 + 0.5j, -1 - 0.5j, 1.2 + 0.7j):
    r = cc.bound_abs(lambda x, y, z0=z0: x - acb(z0.real, z0.imag) + 0 * y, box, sc)
    true_min = min(abs(complex(a, b) - z0) for a in (-1, 1) for b in (-.5, .5)) if not (abs(z0.real) <= 1 and abs(z0.imag) <= .5) else 0
    if abs(z0.real) > 1 and abs(z0.imag) <= .5: true_min = abs(z0.real) - 1
    print(z0, "->", None if r is None else (float(r[0]), float(r[1])), "true min dist", true_min)
# zero on y-dimension boundary & function depending on both
r = cc.bound_abs(lambda x, y: x * x + y * y + 1.0000001, [(-0.1, 0.1), (-1.0, 1.0), (-0.01, 0.01), (0.0, 0.0)], sc)
print("x^2+y^2+1.0000001 (zero just outside at x=+-i(1+5e-8)):", None if r is None else (float(r[0]), float(r[1])))
r = cc.bound_abs(lambda x, y: x * x + y * y + 1, [(-0.1, 0.1), (-1.0, 1.0), (-0.01, 0.01), (0.0, 0.0)], sc)
print("x^2+y^2+1 (zero at corner):", r)
# exactness / containment of float round trips used between processes
random.seed(1)
bad = 0
for _ in range(20000):
    f = (arb(random.uniform(-1, 1)) * arb.pi() * arb(10) ** random.randint(-30, 3)) + arb(0, random.uniform(0, 1) * 10.0 ** random.randint(-40, -5))
    mid = float(f.mid()); rad = cc.up(f.rad()) + abs(mid) * 2.0**-50
    g = arb(mid, rad)
    if not g.contains(f): bad += 1
print("arb -> (mid, rad) -> arb containment failures:", bad, "of 20000")
x = arb(1) / 3
print("float(arb) rounding:", Fraction(float(x)) - Fraction(1, 3), "; up(x) >= x:", all(arb(cc.up(arb.pi() * k)) > arb.pi() * k for k in range(1, 2000)))
print("up(0) =", cc.up(arb(0)), " up(tiny) =", cc.up(arb(1e-320)))
# self kernel divided difference vs direct formula, and sinc semantics
geo = cc.Geometry("G")
mx = 0
for _ in range(2000):
    i = random.randrange(24)
    x = acb(random.uniform(-1.3, 1.3), random.uniform(-0.3, 0.3)); y = acb(random.uniform(-1.3, 1.3), random.uniform(-0.3, 0.3))
    d = geo.self_kernel(i, x, y) - 2 * ((geo.w(i, x) - geo.w(i, y)) / 2).sinh() / (x - y)
    assert d.contains(0)
    mx = max(mx, float(abs(d).upper()))
print("self_kernel == 2 sinh((w(x)-w(y))/2)/(x-y) on 2000 random complex pairs, max enclosure width", mx)
k = geo.self_kernel(6, acb(arb(-0.5, 0.01)), acb(arb(-0.5, 0.01)))
kd = geo.self_kernel(6, acb(-0.5), acb(-0.5)); wp = geo.a[6] * (1 + acb(0, 1) * (2 * geo.b1[6] * (geo.a[6] * arb(-0.5)) + 4 * geo.b2[6] * (geo.a[6] * arb(-0.5)) ** 3))
print("diagonal value K(x,x) = w'(x):", kd, wp, k.contains(kd))
# class structure of G
print("G classes", geo.cls, "reps", geo.reps, "tips distinct:", len({t[3] for t in geo.table}) == 5)
# symmetry: every slit is the image of its class representative under rotation / conjugation
import cmath
for j in range(24):
    rep = geo.reps[geo.cls[j]]
    r, sgn = geo.table[j][0], geo.table[j][1]
    for xv in (-1.0, -0.4, 0.7):
        wr = complex(geo.w(rep, arb(xv))); wj = complex(geo.w(j, arb(xv)))
        img = (wr if sgn == geo.table[rep][1] else wr.conjugate())
        d = (wj - img).imag / (2 * math.pi / 3)
        assert abs(wj.real - img.real) < 1e-13 and abs(d - round(d)) < 1e-12, (j, d)
print("each of the 24 arcs = rotation (multiple of 2pi/3) of its representative or of its conjugate: OK")
