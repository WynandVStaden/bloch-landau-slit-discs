"""Latent-bug demo: classes with EQUAL tips are merged by `tips.index(tip)` (line 90), so only one of them is checked."""
import os, sys, json, shutil, subprocess, math
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = r"C:\Users\wvs00\claude_chat\constants\.venv\Scripts\python.exe"
SRC = os.path.join(ROOT, "certify")
ORDER = sys.argv[1]  # "AB" or "BA"
A, B = '("0", "0.5", "0", "0")', '("0.3", "0.5", "0", "0")'
classes = f"[{A}, {B}]" if ORDER == "AB" else f"[{B}, {A}]"
ADD = ('DESIGNS["ONE"] = ("1", [("0", "0.5", "0", "0")])', 'DESIGNS["ONE"] = ("1", [("0", "0.5", "0", "0")])\nDESIGNS["BAD"] = ("1", ' + classes + ')')
PR = ("    assert arb(lower) < rho_lower\n", "    assert arb(lower) < rho_lower\n    log('RHO_LOWER_BALL ' + rho_lower.lower().str(25) + ' DELTA ' + delta.upper().str(12) + ' classes ' + str(geo.cls) + ' reps ' + str(geo.reps))\n")
FIX = ("            self.cls.append(tips.index(tip))\n", "            self.cls.append([k for k, c in enumerate(classes) if c[0] == alpha and Fraction(c[1]) == tip][0])\n")
def build(tag, subs):
    d = os.path.join(ROOT, "v_" + tag); os.makedirs(d, exist_ok=True)
    src = open(os.path.join(SRC, "conformal_certificate.py")).read()
    for old, new in subs:
        assert src.count(old) == 1, old
        src = src.replace(old, new)
    open(os.path.join(d, "conformal_certificate.py"), "w").write(src)
    shutil.copy(os.path.join(SRC, "inradius_certificate.py"), d)
    return d
def go(d):
    p = subprocess.run([PY, "-u", "conformal_certificate.py", "BAD", "48", "2"], cwd=d, capture_output=True, text=True, timeout=1200)
    out = (p.stdout + p.stderr).splitlines()
    for l in out:
        if "RHO_LOWER" in l or "CERTIFIED" in l or "Error" in l or "sup |U|" in l or "density" in l: print("   ", l[:250])
# 1. honest run with the class assignment fixed: solver density per true class, gives float rho and a certified bound
dfix = build("merged_fix_" + ORDER, [ADD, PR, FIX])
for f in os.listdir(dfix):
    if f.startswith("density_BAD"): os.remove(os.path.join(dfix, f))
print("--- fixed class assignment (truth reference)"); go(dfix)
# 2. adversarial density for the ORIGINAL program: one psi for all 9 slits, tuned so that U = 0 on slit 0 only
dbad = build("merged_bad_" + ORDER, [ADD, PR])
sys.path.insert(0, dbad); os.chdir(dbad)
import conformal_certificate as cc
geo = cc.Geometry("BAD")
print("original program: classes", geo.cls, "reps", geo.reps)
m, N = 24, 600
ks = np.arange(m, 2 * m); xs = np.cos((2 * ks + 1) * np.pi / (4 * m))
ys = np.cos((2 * np.arange(1, N + 1) - 1) * np.pi / (2 * N)); wt = np.pi / N
ns = np.arange(1, 2 * m, 2)
Ty = np.cos(np.outer(ns, np.arccos(ys)))  # T_n(y_l)
Tx = np.cos(np.outer(ns, np.arccos(xs)))
i = 0
wx = np.array([geo.fw(i, x) for x in xs])
Mat = (np.pi / ns)[None, :] * Tx.T
for j in range(geo.n_slits):
    wy = np.array([geo.fw(j, y) for y in ys])
    if j == i:
        t = geo.fa[i] * (xs[:, None] - ys[None, :]) / 2
        lk = np.log(np.abs(geo.fa[i] * np.sinc(1j * t / np.pi)))
    else:
        lk = np.log(np.abs(2 * np.sinh((wx[:, None] - wy[None, :]) / 2)))
    Mat -= wt * lk @ Ty.T
c = np.linalg.solve(Mat, geo.fa[i] * xs)
coef = [0.0] * (2 * m)
for n, v in zip(ns, c): coef[n] = float(v)
json.dump({"design": "BAD", "n": 48, "coefficients": {"0": coef}, "floating_point_conformal_radius": 0.0}, open("density_BAD_n48.json", "w"), indent=1)
print("--- original program with the adversarial (odd, single-class) density"); go(dbad)
