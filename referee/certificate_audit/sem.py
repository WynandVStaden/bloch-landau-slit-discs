from flint import arb, acb, ctx, fmpq
from fractions import Fraction
import flint
print(flint.__version__)
ctx.prec=100
z=acb(0.3,0.4)
print("sinc", z.sinc(), z.sin()/z)
b=arb(1.0,1e-20); print("rad", b.rad(), Fraction(float(b.rad()))>=Fraction(1e-20))
b=arb(0,3e-13); print(b.rad(), Fraction(float(b.rad()))>=Fraction(3e-13))
x=arb(1,0.5)
print("cmp", x>0.4, x>0.6, x<2, x<1.2, (x>0.6) is True)
v=acb(arb(1,0.5),arb(0,0.5)); print(v.abs_lower(), v.abs_upper(), type(v.abs_lower()))
v=acb(arb(0,0.5),arb(0,0.5)); print(v.abs_lower(), v.abs_lower()>0)
print(arb.cos_pi_fmpq(fmpq(1,3)), arb.cos_pi_fmpq(fmpq(7,3)))
u=arb(1)/3; print(u.upper(), type(u.upper()), u.lower(), float(u.upper()))
print(float(Fraction(1e-9)-Fraction(1,10**9)), 2.0**-69, 2**0.5*2.0**-70)
print(u.rad()<arb(2)**(-72))
w=arb(1,0.25); print(w.mid(), w.rad(), abs(arb(-1,0.5)).lower(), arb(3).max(arb(4,1)))
import mpmath as mp
mp.mp.dps=30
print((1+mp.mpf('1e-9'))/mp.mpf('1.523726208330'), (1+mp.mpf(1e-9))/mp.mpf(1.5237262083295537))
