# Literature notes

Compiled on 13 October 2026. Entries marked "read" were checked against the source text; the others come from abstracts, reviews or secondary accounts. **Repeat the search immediately before any submission.**

## Upper bounds for the univalent Bloch–Landau constant

| Bound | Source | Nature | Checked |
|---|---|---|---|
| π/4 ≈ 0.7854 | infinite strip (attributed to Szegő) | closed form | secondary |
| 0.658 | Robinson, Bull. Amer. Math. Soc. 41 (1935) | closed form | read |
| 0.65647 | Goodman, Bull. Amer. Math. Soc. 51 (1945) | closed form | read |
| 0.6564155 | Beller and Hummel, Complex Variables 4 (1985) | numerical | **not read** (not freely available) |
| 0.6563937 | Carroll and Ortega-Cerdà, J. Math. Pures Appl. 92 (2009), arXiv:0806.2282 | closed form, stated correct to 10 decimal places | read |
| 0.656319277272 | Ortega-Cerdà and Pridhnani, Contemp. Math. 505 (2010), arXiv:0809.2483 | numerical; the authors state that no rigorous error analysis was done | read (arXiv version 1) |

Best lower bound found: 0.5708858 (Skinner, 2009).

The 2010 value is rarely quoted. Later surveys and preprints, including a preprint of March 2026 (arXiv:2603.05047), give 0.6563937 as the best upper bound. The catalogue of optimisation constants maintained at github.com/teorth/optimizationproblems lists only the trivial bound 1 for this constant (entry 57c) as of the date above.

## The domains behind the earlier bounds

- **Robinson (1935).** Disc of radius about 2.9562 with radial slits from radius 1 (at 0° and ±120°) and from radius 2 (at ±60° and 180°) out to the circle.
- **Goodman (1945).** Radial slits added in generations, each generation bisecting the sectors of the previous one, with tips at radius 1, 2, 2.9557 and then 1 + 1/sin(π/(3·2^(n−1))).
- **Beller and Hummel (1985).** As Goodman for the first two generations; the third-generation slits are moved from argument π/6 to about 0.3931 (according to later accounts by Carroll). The tip radius and other details were not established because the paper was not available.
- **Carroll and Ortega-Cerdà (2009).** Disc of radius 4.0546358 with radial slits from radius 1 and 2 on the six axes, plus six non-radial arcs satisfying a "harmonic symmetry" condition that is necessary for an extremal domain. The conformal radius comes from an explicit formula.
- **Ortega-Cerdà and Pridhnani (2010).** The same axes in a disc of radius 5.1195152501, with two harmonically symmetric arcs in each half-sector; the conformal radius is computed numerically.

## Searches that found nothing newer

arXiv full text ("univalent Bloch", "Bloch–Landau"); zbMATH; the citation graphs of the papers above (Semantic Scholar, OpenAlex); GitHub repositories, code and issues; Zenodo; OpenReview; MathOverflow and Mathematics Stack Exchange; the pull requests of the optimisation-constants catalogue; the Hayman–Lingham problem list (problem 5.10).

Not consulted: the full text of Beller and Hummel; MathSciNet; Google Scholar.

## Possible venues and their AI policies (to be re-checked)

- **arXiv, category math.CV.** AI tools may not be listed as authors, and significant use must be reported.
- **Computational Methods and Function Theory** (Springer). Covers complex analysis and computational methods; AI use must be declared.
- **Experimental Mathematics**, **Complex Variables and Elliptic Equations** (Taylor & Francis). Use of generative AI must be acknowledged, naming the tool and its purpose; code should be made available.
- **The optimisation-constants catalogue.** Accepts pull requests in Markdown; AI use must be disclosed and the contributor must have verified the content. An entry should also record the 2009 and 2010 values, which the table currently lacks.
