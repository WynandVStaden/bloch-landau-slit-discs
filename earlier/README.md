# Earlier models

Superseded versions of the slit model and the search, kept so that the earlier results can be reproduced. None of them is needed for the proof or for domain G.

| File | Contents |
|---|---|
| `bloch_curved.py` | slits with one bend parameter: inradius routine (domain B was found with this model) |
| `bloch_evolve.py` | differential-evolution search over straight radial slits (domain A; the six- and twelve-slit checks) |
| `bloch_evolve2.py` | the same search over slits with one bend |
| `bloch_verify.py` | floating-point checks of a straight-slit domain |
| `bloch_verify_curved.py` | floating-point checks of a domain with one-bend slits |

Run them from inside this folder; they use the shared modules (`bloch_solver.py`, `bloch_design.py`, `bloch_residual.py`) in the parent folder.
