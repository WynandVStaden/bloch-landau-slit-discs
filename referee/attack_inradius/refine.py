# Part A stage 2: mp refinement (50 digits) of the candidate maxima
from geom import *
import pickle, sys
res = pickle.load(open("search_res.pkl", "rb"))
cands = []
for c in res:
    for f in c['found']:
        if not any(abs(complex(f[1], f[2]) - complex(g[1], g[2])) < 1e-6 for g in cands): cands.append(f)
cands.sort(reverse=True)
print("distinct NM optima:", len(cands))
def solve(c0, contacts):
    # contacts: list of (k, s, kind)
    nint = [i for i, ct in enumerate(contacts) if ct[2] == 'int']
    def F(*v):
        c = mpc(v[0], v[1]); r = v[2]; out = []; j = 3
        for (k, s, kind, _d) in contacts:
            if kind == 'circle': out.append(abs(c) - (R_MP - r)); continue
            sl = SL_MP[k]
            if kind == 'int':
                s = v[j]; j += 1
                g, g1 = gam(sl, s, 1)
                out.append(mp.re((g - c)*mp.conj(g1)))
            else:
                g = gam(sl, -sl['a'] if kind == 'tip' else mpf(0))
            out.append(abs(g - c)**2 - r*r)
        return out
    x0 = [mp.re(c0), mp.im(c0), contacts[0][3]] + [contacts[i][1] for i in nint]
    return mp.findroot(F, x0, tol=mpf(10)**-44, maxsteps=80, verify=False), nint
out = []
for (d0, x, y) in cands:
    c = mpc(x, y)
    L = dist_mp(c)
    near = [t for t in L if t[0] - L[0][0] < 1e-5]
    desc0 = [(t[1], t[3], float(t[0] - L[0][0])) for t in near]
    if len(near) < 3:
        print("TWO-CONTACT candidate d=%.12f c=(%.9f,%.9f)" % (d0, x, y), desc0); 
        out.append(dict(c=c, r=L[0][0], contacts=near, ok=False, note='2-contact')); continue
    best = None
    import itertools
    for trio in itertools.combinations(near, 3):
        contacts = [(t[1], t[2], t[3], t[0]) for t in trio]
        try:
            sol, nint = solve(c, contacts)
        except Exception as e:
            print('   solve exception', repr(e)[:200]); continue
        cc = mpc(sol[0], sol[1]); rr = sol[2]
        LL = dist_mp(cc)
        if abs(LL[0][0] - rr) < mpf(10)**-38:
            best = (cc, rr, LL); break
    if best is None:
        print("REFINE FAIL d=%.12f c=(%.9f,%.9f)" % (d0, x, y), desc0); out.append(dict(c=c, r=L[0][0], contacts=near, ok=False, note='fail')); continue
    cc, rr, LL = best
    con = [t for t in LL if t[0] - rr < mpf(10)**-30]
    # local max test: contact directions surround the centre
    angs = sorted(float(mp.arg((gam(SL_MP[t[1]], t[2]) if t[1] >= 0 else cc/abs(cc)*R_MP) - cc)) for t in con)
    gaps = [angs[(i + 1) % len(angs)] - angs[i] + (2*math.pi if i == len(angs) - 1 else 0) for i in range(len(angs))]
    out.append(dict(c=cc, r=rr, contacts=con, ok=True, maxgap=max(gaps), slack=float(LL[len(con)][0] - rr)))
out_ok = [o for o in out if o['ok']]
# dedupe
uniq = []
for o in sorted(out_ok, key=lambda o: -o['r']):
    if not any(abs(o['c'] - u['c']) < mpf(10)**-20 for u in uniq): uniq.append(o)
print("distinct refined critical points:", len(uniq))
for o in uniq:
    cc = o['c']; ismax = o['maxgap'] < math.pi
    print("r=%s  r-1=%s  c=(%s, %s) |c|=%s arg=%s  %s maxgap=%.4f next-slack=%.3g" % (mp.nstr(o['r'], 32), mp.nstr(o['r'] - 1, 8), mp.nstr(cc.real, 32), mp.nstr(cc.imag, 32), mp.nstr(abs(cc), 12), mp.nstr(mp.arg(cc), 10), "LOCALMAX" if ismax else "not-max", o['maxgap'], o['slack']))
    for t in o['contacts']:
        if t[1] < 0: print("      circle")
        else:
            sl = SL_MP[t[1]]; g = gam(sl, t[2])
            print("      slit#%d cls%d r%d eps%+d %s s/a=%s  pt=(%s,%s)" % (t[1], sl['cls'], sl['r'], sl['eps'], t[3], mp.nstr(t[2]/sl['a'], 15), mp.nstr(g.real, 25), mp.nstr(g.imag, 25)))
pickle.dump(uniq, open("maxima.pkl", "wb"))
print("count localmax r>0.999:", sum(1 for o in uniq if o['maxgap'] < math.pi and o['r'] > mpf('0.999')))
print("max r - 1 =", mp.nstr(max(o['r'] for o in uniq) - 1, 20))
