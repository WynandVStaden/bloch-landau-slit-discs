import json, mpmath as mp, numpy as np
from fractions import Fraction
mp.mp.dps = 40
cert = json.load(open(r"C:\Users\wvs00\claude_chat\bloch_record\certify\conformal_certificate_G.json"))
dens = json.load(open(r"C:\Users\wvs00\claude_chat\bloch_record\certify\density_G_n64.json"))
print("coefficients in certificate == density file:", cert["coefficients"] == dens["coefficients"])
print("delta_upper =", repr(cert["delta_upper"]), " rho_lower =", repr(cert["conformal_radius_lower"]))
R = mp.mpf(Fraction("5.001117628756668").numerator)/Fraction("5.001117628756668").denominator
tips = ["1","1.9999999974092333","2.789342409403121","3.473381248030611","4.342874184977691"]; mult = [3,3,6,6,6]
c = mp.log(R)
for k,(t,m) in enumerate(zip(tips,mult)):
    a = mp.log(R/(mp.mpf(Fraction(t).numerator)/Fraction(t).denominator))
    c1 = mp.mpf(dens["coefficients"][str(k)][1])
    c -= mp.pi/4*m*a*c1
    print("  class", k, "a =", mp.nstr(a, 8), "c_1 =", mp.nstr(c1, 17))
print("c =", mp.nstr(c, 20), " exp(c) =", mp.nstr(mp.exp(c), 20))
d = mp.mpf(cert["delta_upper"])
print("exp(c - delta_upper) =", mp.nstr(mp.exp(c-d), 20), " >= 1.523726208329 ?", mp.exp(c-d) >= mp.mpf("1.523726208329"))
print("exp(c - 2.656e-8)    =", mp.nstr(mp.exp(c-mp.mpf("2.656e-8")), 20), "; exp(c - 2.66e-8) =", mp.nstr(mp.exp(c-mp.mpf("2.66e-8")), 20))
q = (1+mp.mpf(10)**-9)/mp.mpf("1.523726208329")
print("(1+1e-9)/1.523726208329 =", mp.nstr(q, 20), " < 0.65628589 ?", q < mp.mpf("0.65628589"))
print("float ratio 1/exp(c) =", mp.nstr(1/mp.exp(c), 15), "; proved - float =", mp.nstr(q-1/mp.exp(c), 5))
# six-slit closed form
s = mp.mpf(1)/27; k = lambda z: z/(1+z)**2
target = k(mp.mpf(-8)/27)/(4*s/(1+s)**2)
beta = -mp.findroot(lambda w: k(w)-target, -0.5)
print("six-slit: beta =", beta, " rho =", mp.nstr(3*(4*s/(1+s)**2)**(mp.mpf(1)/3)*(4*beta/(1+beta)**2)**(mp.mpf(1)/3), 14), " (paper 1.520743036313); 8/9 =", mp.nstr(mp.mpf(8)/9, 13))
# pairwise slit distances (z-plane) and min modulus
import ubuild as ub
xs = np.linspace(-1, 0, 4001)
Z = [np.exp(ub.wj(j, xs)) for j in range(ub.NS)]
dm = min(np.abs(Z[i][:,None]-Z[j][None,::4]).min() for i in range(ub.NS) for j in range(i))
print("min pairwise slit distance (z-plane):", dm, "; min |z| on K:", min(np.abs(z).min() for z in Z))
# w-plane: arcs of different slits on the cylinder
W = [ub.wj(j, np.linspace(-1,1,4001)) for j in range(ub.NS)]
def dc(a,b):
    d = a[:,None]-b[None,::4]; return np.hypot(d.real, (d.imag+np.pi)%(2*np.pi)-np.pi).min()
print("min pairwise arc distance on the cylinder (w-plane):", min(dc(W[i],W[j]) for i in range(ub.NS) for j in range(i)))
# self-distance modulo 2 pi i (zeros of K_i need w_i(x)-w_i(y) in 2 pi i Z\{0}) : Im range of each arc
print("max Im-extent of an arc:", max(np.ptp(w.imag) for w in W), " (< 2pi so no self-overlap mod 2 pi i on the real parameter)")
