"""Mutation harness: builds attack_code/v_<tag>/ with a patched copy and runs it."""
import os, shutil, subprocess, sys, json
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = r"C:\Users\wvs00\claude_chat\constants\.venv\Scripts\python.exe"
SRC = os.path.join(ROOT, "certify")
ALWAYS = [("    assert arb(lower) < rho_lower\n",
           "    assert arb(lower) < rho_lower\n    log('RHO_LOWER_BALL ' + rho_lower.lower().str(25) + ' DELTA ' + delta.upper().str(12))\n")]
def run(tag, subs, design, n, procs=2, density_edit=None, timeout=1500, keep=5):
    d = os.path.join(ROOT, "v_" + tag)
    os.makedirs(d, exist_ok=True)
    src = open(os.path.join(SRC, "conformal_certificate.py")).read()
    for old, new in ALWAYS + list(subs):
        assert src.count(old) == 1, (tag, old, src.count(old))
        src = src.replace(old, new)
    open(os.path.join(d, "conformal_certificate.py"), "w").write(src)
    shutil.copy(os.path.join(SRC, "inradius_certificate.py"), d)
    f = f"density_{design}_n{n}.json"
    dens = json.load(open(os.path.join(SRC, f)))
    if density_edit:
        density_edit(dens)
    json.dump(dens, open(os.path.join(d, f), "w"), indent=1)
    try:
        p = subprocess.run([PY, "-u", "conformal_certificate.py", design, str(n), str(procs)], cwd=d, capture_output=True, text=True, timeout=timeout)
        out = (p.stdout + p.stderr).strip().splitlines()
    except subprocess.TimeoutExpired as e:
        out = ["TIMEOUT"] + ((e.stdout or b"").decode(errors="replace").splitlines() if isinstance(e.stdout, bytes) else (e.stdout or "").splitlines())
    open(os.path.join(d, f"out_{design}.log"), "w").write("\n".join(out))
    res = [l for l in out if l.startswith("RHO_LOWER") or "CERTIFIED" in l or "Error" in l or "exact value" in l or l.startswith("TIMEOUT")]
    print(f"=== {tag} [{design}]")
    for l in (res if res else out[-keep:]):
        print("   ", l[:230])
    sys.stdout.flush()
    return out
