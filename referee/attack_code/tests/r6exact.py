from mpmath import mp, mpf, sqrt, cbrt
mp.dps = 40
R, t = mpf(3), mpf(2)
k = lambda z: z / (1 - z) ** 2
s = R ** -3
# conformal map g of D \ [s,1) onto D, g(0)=0: k(g(z)) = c * zeta(z), derived independently:
# k maps D onto C \ (-inf,-1/4]; D\[s,1) -> C \ ((-inf,-1/4] U [k(s),inf)); Moebius zeta = w/(1 - w/k(s)) sends k(s)->inf,
# inf -> -k(s), so the image is C \ (-inf... ] one slit from -k(s) through inf to -1/4/(1+1/(4k(s))) -- single slit (-inf, -ks/(4ks+1)]
ks = k(s)
x = (t / R) ** 3
w = k(-x)
zeta = w / (1 - w / ks)
tip = -ks / (4 * ks + 1)          # image slit is (-inf, tip]; rescale so tip -> -1/4
kk = abs(zeta) * (mpf(1) / 4) / (-tip)

beta = ((1 - 2 * kk) - sqrt(1 - 4 * kk)) / (2 * kk)
rho1 = 4 * s / (1 + s) ** 2      # g'(0)^-1 for first map: check: k(g) ~ g, zeta ~ z, scale (1/4)/(-tip) = (4ks+1)/(4ks)
print("check first factor", rho1, 1 / ((4 * ks + 1) / (4 * ks)))
rho = R * cbrt(rho1) * cbrt(4 * beta / (1 + beta) ** 2)
print("R6 exact", rho)
print("ONE exact", mpf(8) / 9)
