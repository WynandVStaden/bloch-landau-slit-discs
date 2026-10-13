# Univalent Bloch–Landau constant: slit-disc upper bounds

**Private working repository. Nothing here has been submitted or published.**

Numerical evidence that the univalent Bloch–Landau constant satisfies U ≤ 0.656285866, against the best proved bound 0.6563937 (Carroll and Ortega-Cerdà, 2009) and the best published numerical value 0.656319277 (Ortega-Cerdà and Pridhnani, 2010). A computer-assisted proof (ball arithmetic) of U ≤ 0.65628589 is in `certify/`; it runs in under a minute and is **still being audited**.

| File | Purpose |
|---|---|
| `NOTE.md` | draft write-up |
| `notes/WALKTHROUGH.md` | the method step by step, and how to run every check |
| `notes/LESSONS.md` | the mathematics from the basics |
| `notes/LITERATURE.md` | earlier bounds, searches done, possible venues |
| `designs.py` | the domains A, B and G; recomputes their ratios and draws `domains.png` |
| `bloch_solver.py` | conformal radius of a slit disc (spectral integral-equation method); self-test against exact formulas |
| `bloch_explicit.py` | exact conformal radius for six slits on the symmetry axes |
| `bloch_design.py`, `bloch_curved.py`, `bloch_curved2.py` | inradius for straight slits, one-bend slits and two-parameter curved slits |
| `bloch_verify.py`, `bloch_verify_curved.py`, `bloch_verify_curved2.py` | checks of a domain (brute force, convergence, rotation) |
| `bloch_residual.py` | a posteriori error bound for the conformal radius |
| `bloch_evolve.py`, `bloch_evolve2.py` | the differential-evolution searches |
| `certify/inradius_certificate.py` | proof that the inradius is at most 1 + 1e-9 |
| `certify/conformal_certificate.py` | proof of a lower bound for the conformal radius (try `ONE` and `R6` first: exactly solvable tests) |
| `notes/CERTIFICATE.md` | the mathematics of the certificate |
| `referee/` | scripts written by independent reviewers |
| `logs/` | search and verification logs |

Requirements: Python 3 with numpy, scipy, matplotlib (`certify/` also needs python-flint).

```
python designs.py
cd certify
python inradius_certificate.py G 1e-9
python conformal_certificate.py G 64
```

Prepared with AI assistance (Claude); see the disclosure in `NOTE.md`.
