"""Validation of Method A on exact cases."""
import numpy as np
from methodA import solve, postprocess
from geom import R0, CLASSES0


def run(name, R, classes, exact, Ns=(8, 16, 32, 64)):
    for N in Ns:
        f, C, A, rhs = solve(R, classes, N)
        o = postprocess(R, classes, N, f, C)
        print(f'{name:38s} N={N:3d} rho={o["rho"]:.15f} exact={exact:.15f} err={o["rho"] - exact:+.2e}')
    return o


rad = lambda th, t, imgs: dict(theta=th, t=t, b1=0.0, b2=0.0, images=imgs)
# 1. unit disc minus [s,1)
for s in (0.5, 0.1, 0.01):
    o = run(f'unit disc minus [{s},1)', 1.0, [rad(0.0, s, [(0.0, False)])], 4 * s / (1 + s)**2)
    # exact harmonic measure of the slit from 0: Koebe-type map sends slit to arc of angle 2*acos((1-6s+s^2)/-(1+s)^2)...
    # computed below independently via the explicit map
    k = 4 * s / (1 + s)**2          # k z/(1-z)^2 composition: slit endpoints map to e^{+-i alpha}, cos(alpha/2)...
    # w/(1+w)^2 = k z/(1+z)^2 ; z=1 -> w/(1+w)^2 = k/4 -> w = e^{i al}: 1/(2+2cos al) = k/4
    al = np.arccos(2 / k - 1) if abs(2 / k - 1) <= 1 else None
    # slit = arc between -al and al THROUGH w=1?  w=1 gives 1/4 >= k/4, image of z=1 is w=e^{+-i al}; slit tip z=s -> w=1
    hm_exact = (2 * np.arccos(2 / k - 1)) / (2 * np.pi) if al is not None else float('nan')
    print(f'    harmonic measure of slit: {o["hm_per_slit"][0]:.14f}  exact {hm_exact:.14f}')

# 2. radius R minus 3 symmetric radial slits [t,R)
R = R0
rot3 = [(2 * np.pi * r / 3, False) for r in range(3)]
for t in (1.0, 3.0):
    s = (t / R)**3
    run(f'R-disc minus 3 slits [{t},R)', R, [rad(0.0, t, rot3)], R * (4 * s / (1 + s)**2)**(1 / 3))

# 3. radius R minus 6 radial slits: [t1,R) at angles 2pi k/3 and [t2,R) at pi/3+2pi k/3  (z^3 + two-slit map)
t1, t2 = CLASSES0[0][1], CLASSES0[1][1]
p, q = (t1 / R)**3, (t2 / R)**3
Aa, Bb = (p + 1 / p) / 2, (q + 1 / q) / 2
exact6 = R * (2 / (Aa + Bb))**(1 / 3)
run('R-disc minus 6 radial slits (cl.1+2)', R, [rad(0.0, t1, rot3), rad(np.pi / 3, t2, rot3)], exact6)
# 4. same, unsymmetric tips, no symmetry reduction used (each slit its own class)
cl = [rad(2 * np.pi * r / 3, t1, [(0.0, False)]) for r in range(3)] + [rad(np.pi / 3 + 2 * np.pi * r / 3, t2, [(0.0, False)]) for r in range(3)]
run('same, 6 independent classes', R, cl, exact6, Ns=(16, 32))
