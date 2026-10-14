# Reproducing the results

This guide lists every computation behind the result, in order of importance, with the command to run and the output to expect. Run the commands from the repository root unless stated otherwise.

## 0. Set-up

Python 3.10 or later.

```
python -m venv .venv
.venv\Scripts\activate          (Windows)      source .venv/bin/activate   (Linux, macOS)
pip install -r requirements.txt
```

`requirements.txt` installs numpy, scipy, matplotlib and python-flint (the interface to the Arb ball-arithmetic library).

## 1. The proof (about one minute)

These two programs are the whole proof. Nothing else in the repository is needed for the theorem.

### 1.1 Inradius

```
cd certify
python inradius_certificate.py G 1e-9
```

Expected:

```
stage 2: exact verification of every leaf
  outside disc 2799, outside sector 464, covered by a witness 1730 (610 distinct witnesses)
  smallest cell side 2^-29
CERTIFIED: inradius of design G <= 1 + 1e-9
```

Stage 1 is a floating-point search for a covering by cells and is not trusted. Stage 2 re-checks every cell in exact arithmetic; only stage 2 matters. The argument is Section 3 of the paper.

### 1.2 Conformal radius

```
python conformal_certificate.py G 64
```

The second argument is the number of coefficients per slit; a third optional argument sets the number of processes (default 8). Expected last lines:

```
delta <= 2.6563e-08   (sup of U over the slits)
c = L - (pi/4) sum a_j c_j1 >= 0.4211588143613
CERTIFIED: conformal radius of design G >= 1.523726208329
```

The program reads the density coefficients from `density_G_n64.json`. If that file is deleted it recomputes them with the floating-point solver; the proof is valid for any coefficients, which only influence how small `delta` comes out. The argument is Section 4 of the paper.

### 1.3 Conclusion

(1 + 10⁻⁹) / 1.523726208329 = 0.6562858836… < 0.65628589.

### 1.4 Testing the certificate itself

The conformal-radius program can be run on two domains whose conformal radius is known exactly. A correct program must return a lower bound just below the exact value.

```
python conformal_certificate.py ONE 32
python conformal_certificate.py R6 48
```

| Design | Domain | Exact value | Expected certified bound |
|---|---|---|---|
| `ONE` | unit disc minus the segment from 1/2 to 1 | 8/9 = 0.888888888889 | 0.888888888327 |
| `R6` | disc of radius 3 with six radial slits | 1.520743036313 | 1.520743034966 |

Things worth trying, to convince yourself that the program cannot be fooled:

- run design G with fewer coefficients (`python conformal_certificate.py G 32`): `delta` grows and the bound drops, but stays below the true value;
- edit a coefficient in `density_G_n64.json`: `delta` grows and the bound stays below the true value;
- ask the inradius program for an impossible bound (`python inradius_certificate.py G -1e-6`): after a minute or two it stops with `NOT CERTIFIED`.

## 2. Floating-point values of the three domains (a few seconds)

```
python designs.py
```

Expected:

```
A: 18 straight slits: ... ratio 0.656364185 ...
B: 18 slits, 12 curved: ... ratio 0.656308522 ...
G: 24 slits, 18 curved: R = 5.001117628756668, inradius 1.0000000000, conformal radius 1.523726248805, ratio 0.656285866 ...
figure written to domains.png
```

Domain G is the one in the theorem. A and B are simpler domains found earlier; B also beats the earlier values (it certifies as U < 0.65630853 with `python conformal_certificate.py B 64`).

## 3. The solver against exact formulas (a few seconds)

```
python bloch_solver.py
```

Prints the computed conformal radius next to the exact one for a disc with one slit and with six slits, at several resolutions. The errors should fall to about 1e-16.

## 4. Extended floating-point checks of domain G (tens of minutes)

```
python bloch_verify_curved2.py 5.001117628756668 1.9999999974092333 0.029930188727091613 2.789342409403121 0.9026921410641155 -0.3301459494660587 1.0157094972537686 3.473381248030611 -1.2656124700253997 0.9794269477754911 0.5757776442359243 4.342874184977691 0.03225521024585182 2.9360081049550737
```

The arguments are the outer radius, the tip radius of the slits at 60°, then (angle, tip radius, first bend, second bend) for the three curved classes. The script reports:

- that no two slits cross, and the smallest distance between slits (0.14842);
- the inradius from two grids and from 30 million random points;
- the conformal radius at four resolutions and for a rotated copy of the domain (all 1.523726248805);
- the largest value of the computed Green function on the slits.

A saved run is in `logs/verify_G_final.log`.

## 5. Independent recomputations

The folder `referee/` contains programs written independently of the main code, from a description of the claim only, using different methods.

| Folder | What it checks | Method |
|---|---|---|
| `referee/G_24_slits/` | domain G: geometry, inradius, conformal radius | branch and bound over the whole disc; boundary integral equation in the z-plane |
| `referee/B_18_slits/` | domain B, in the same way | as above |
| `referee/straight_18_slits/` | an earlier straight-slit domain | Schwarz–Christoffel parameter problem |
| `referee/certificate_audit/` | the proof and the two certificate programs | re-derivation, perturbation attacks, independent high-precision evaluation |
| `referee/attack_rho/` | the conformal radius, from the description of the domain only | three boundary-integral solvers, one in 34-digit arithmetic (value 1.523726248805208911047994), and a random-walk estimate |
| `referee/attack_inradius/` | the inradius and its certificate | 50-digit search for all locally largest discs; an independent exact verifier of the stored covering; mutation tests |
| `referee/attack_code/` | soundness of the conformal-radius program | line-by-line audit; corrupted inputs and injected bugs |
| `referee/attack_theory/` | the lemmas of the paper | re-derivation and numerical tests of each lemma |

These scripts were written for one-off use and expect to be run from their own folder; some refer to file locations on the machine where they were written.

## 6. Re-running the search (hours)

The search that found G:

```
set DE_WORKERS=4            (Windows; use  export DE_WORKERS=4  elsewhere)
python -u bloch_curved2.py 3 0 3000
```

Arguments: number of free slit classes, random seed, maximum number of generations. Set the environment variables `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS` and `MKL_NUM_THREADS` to 1 when using several workers. The search is stochastic and need not return the same domain; any domain it returns can be checked with steps 1 and 4 after adding its parameters to the `DESIGNS` table in `certify/inradius_certificate.py`.

`bloch_evolve.py` and `bloch_evolve2.py` are the earlier searches over straight slits and slits with one bend.

## 7. Rebuilding the paper

```
cd paper
pdflatex bloch_slit_discs.tex
pdflatex bloch_slit_discs.tex
```
