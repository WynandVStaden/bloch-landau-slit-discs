# Literature notes

Compiled on 13 October 2026 and fact-checked on 14 October 2026. Entries marked "read" were checked against the source text; the others come from abstracts, reviews or secondary accounts. All bibliographic details (authors, titles, volumes, pages) were confirmed against publisher or Crossref records. **Repeat the search immediately before any submission.**

## Upper bounds for the univalent Bloch–Landau constant

| Bound | Source | Nature | Checked |
|---|---|---|---|
| π/4 ≈ 0.7854 | infinite strip (attributed to Szegő) | closed form | secondary |
| 0.658 | Robinson, Bull. Amer. Math. Soc. 41 (1935) | not checked | value confirmed from later papers; original **not read** |
| 0.65647 | Goodman, Bull. Amer. Math. Soc. 51 (1945) | not checked | value confirmed from later papers; original **not read** |
| 0.6564155 | Beller and Hummel, Complex Variables 4 (1985) | numerical | **not read** (not freely available) |
| 0.6563937 | Carroll and Ortega-Cerdà, J. Math. Pures Appl. 92 (2009), arXiv:0806.2282 | explicit formula evaluated numerically, stated correct to 10 decimal places | read |
| 0.656319277272 | Ortega-Cerdà and Pridhnani, Contemp. Math. 505 (2010), arXiv:0809.2483 | numerical; the authors state that no rigorous error analysis was done | read (arXiv version 1) |

Best lower bound found: 0.5708858 (Skinner, 2009).

The 2010 value is rarely quoted. Later surveys and preprints, including a preprint of March 2026 (arXiv:2603.05047), give 0.6563937 as the best upper bound. The catalogue of optimisation constants maintained at github.com/teorth/optimizationproblems lists only the trivial bound 1 for this constant (entry 57c) as of the date above.

## The domains behind the earlier bounds

- **Robinson (1935).** Disc of radius about 2.9562 with radial slits from radius 1 (at 0° and ±120°) and from radius 2 (at ±60° and 180°) out to the circle.
- **Goodman (1945).** Radial slits added in generations, each generation bisecting the sectors of the previous one, with the first tips at radius 1, 2 and about 2.9557 (as described by Carroll and Ortega-Cerdà).
- **Beller and Hummel (1985).** As Goodman for the first two generations; the third-generation slits are no longer bisectors (Carroll and Ortega-Cerdà describe them as "far from being bisectors"). An earlier search reported the argument as about 0.3931, which matches the 0.393 found by our own twelve-slit search (`logs/evolve_3_twelve_straight_slits.log`), but this figure was not re-verified. The paper itself was not available.
- **Carroll and Ortega-Cerdà (2009).** Disc of radius 4.0546358 with radial slits from radius 1 and 2 on the six axes, plus six separate curved arcs running from an interior point to the outer circle. The arcs are "harmonically symmetric", a condition shown by Carroll (Comput. Methods Funct. Theory 8 (2008), extending Jenkins) to be necessary for boundary arcs of an extremal domain that do not touch an extremal disc. The conformal radius comes from an explicit theta-function formula, evaluated numerically ("correct to 10 decimal places").
- **Ortega-Cerdà and Pridhnani (2010).** The same axes in a disc of radius 5.1195152501, with two harmonically symmetric arcs in each 60° sector, one attached to the radial slit from radius 1 and one ending on the outer circle; the conformal radius is computed numerically. The authors write that they are confident of the first 10 digits but "have not done a rigorous error analysis".

## Searches that found nothing newer

arXiv full text ("univalent Bloch", "Bloch–Landau"); zbMATH; the citation graphs of the papers above (Semantic Scholar, OpenAlex); GitHub repositories, code and issues; Zenodo; OpenReview; MathOverflow and Mathematics Stack Exchange; the pull requests of the optimisation-constants catalogue; the Hayman–Lingham problem list (problem 5.10).

Not consulted: the full texts of Robinson, Goodman, and Beller and Hummel; MathSciNet; Google Scholar.

## Related recent activity

- F. Wikström posted two computer-assisted papers in August and September 2026 on neighbouring constants (arXiv:2608.17660 on the classical Bloch constant, arXiv:2609.23236 on the Landau and locally univalent constants). They do not concern the univalent upper bound.
- Ali and Azim (arXiv:2603.05047, March 2026) quote 0.6563937 and Skinner's 0.5708858 as the best available bounds.
- The names "univalent Bloch constant" and "schlicht Bloch constant" are used in the literature as often as "univalent Bloch–Landau constant"; search for all three.

## Possible venues and their AI policies (to be re-checked)

- **arXiv, category math.CV.** AI tools may not be listed as authors, and significant use must be reported.
- **Computational Methods and Function Theory** (Springer). Covers complex analysis and computational methods. Its AI policy was not verified.
- **Experimental Mathematics**, **Complex Variables and Elliptic Equations** (Taylor & Francis). The publisher's policy requires naming the AI tool with its version and explaining how and why it was used, and it restricts the use of generative AI for producing text, code and research outputs without rigorous revision. Given how this work was produced, read the current policy carefully before choosing these journals.
- **The optimisation-constants catalogue.** Accepts pull requests in Markdown; AI use must be disclosed and the contributor must have verified the content. An entry should also record the 2009 and 2010 values, which the table currently lacks.
