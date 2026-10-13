"""Audit helper: imports the certificate module without touching certify/, adds test designs (also in worker processes)."""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "certify"))
import conformal_certificate as cc
EXTRA = json.loads(os.environ.get("AUD_DESIGNS", "{}"))
for k, v in EXTRA.items():
    cc.DESIGNS[k] = (v[0], [tuple(c) for c in v[1]])
_orig_geo_init = cc.Geometry.__init__
def _geo_init(self, name):
    # treat any design whose name starts with ONE as a single-slit design, like the built-in ONE
    if name.startswith("ONE") and name != "ONE":
        cc.DESIGNS["ONE"], saved = cc.DESIGNS[name], cc.DESIGNS["ONE"]
        try: _orig_geo_init(self, "ONE")
        finally: cc.DESIGNS["ONE"] = saved
    else:
        _orig_geo_init(self, name)
cc.Geometry.__init__ = _geo_init
def init(name, coef):
    cc._init(name, coef)
