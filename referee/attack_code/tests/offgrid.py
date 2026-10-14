import sys, os, json, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "certify"))
import conformal_certificate as cc
from flint import ctx
ctx.prec = cc.PREC
geo = cc.Geometry("G")
coef = {int(c): v for c, v in json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "certify", "density_G_n64.json")))["coefficients"].items()}
pl = cc.plan(geo, coef, lambda s: None)
cc._init("G", coef)
claimed = {0: 2.656e-8, 3: 3.548e-9, 6: 1.758e-8, 12: 3.758e-9, 18: 7.478e-10}
random.seed(3)
for i in geo.reps:
    Mp = 20011  # prime: nodes of T_Mp are not nodes of T_M
    ks = random.sample(range(Mp // 2 + 1, Mp), 150) + list(range(Mp // 2 + 1, Mp // 2 + 40)) + list(range(Mp - 40, Mp))
    _, res = cc._eval_nodes((i, Mp, ks, pl[i]["quad"]))
    mx = max(abs(m) + r for _, m, r in res)
    print(f"rep {i}: max |f| over {len(ks)} off-grid points <= {mx:.3e}; certified sup {claimed[i]:.3e}; ok {mx <= claimed[i]}", flush=True)
