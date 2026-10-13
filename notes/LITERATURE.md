# Literature notes (compiled 2026-10-13)

Compiled from an AI-run literature review; items marked "read" were checked against the source text, the rest come from abstracts, reviews or secondary accounts. **Re-run the search immediately before any submission.**

## Upper bounds for the univalent Bloch–Landau constant

| Bound | Source | Nature |
|---|---|---|
| $\pi/4\approx0.7854$ | strip (Szegő, cited by Landau and Robinson) | closed form |
| $0.658$ | Robinson, Bull. AMS 41 (1935) — read | closed form |
| $0.65647$ | Goodman, Bull. AMS 51 (1945) — read | Koebe compositions |
| $0.6564155$ | Beller & Hummel, Complex Variables 4 (1985) — **not read** (paywalled) | numerical |
| $0.6563937$ | Carroll & Ortega-Cerdà, JMPA 92 (2009), arXiv:0806.2282 — read | closed form, "correct to 10 decimal places" |
| $0.656319277272$ | Ortega-Cerdà & Pridhnani, Contemp. Math. 505 (2010), arXiv:0809.2483 — read (arXiv v1) | numerical; "we have not done a rigorous error analysis" |

Lower bound: Skinner 2009, $\mathcal U>0.5708858$.

The 2010 value is almost never quoted. Finch's *Mathematical Constants* addenda, arXiv:2603.05047 (March 2026) and the optimization-constants catalogue all give $0.6563937$ (the catalogue entry 57c lists only the trivial bound 1).

## The classical domains

* **Robinson.** Disc of radius $a\approx2.95617$ with radial slits from radius 1 (angles $0,\pm120^\circ$) and from radius 2 (angles $\pm60^\circ,180^\circ$) to the circle.
* **Goodman.** Radial slits along the $3\cdot2^{n-1}$-th roots of unity, generation by generation, each bisecting the previous sectors; tip radii $1$, $2$, $2.9557$, then $1+\csc(\pi/(3\cdot2^{n-1}))$ ($\approx4.864,\ 8.661,\dots$). Proves $0.65646<B_\infty<0.65647$ for the limiting domain.
* **Beller–Hummel.** Same first two generations; third-generation radial slits at argument about $0.3931$ instead of $\pi/6$ (per Carroll, Irish Math. Soc. Bull. 47 (2001)). zbMATH review: domain "bounded by two circular arcs about 0 and four radial segments", nonlinear parameter problem solved numerically (Schwarz–Christoffel). Tip radius and truncation not established.
* **Carroll–Ortega-Cerdà.** Disc $|z|<R$, $R=4.0546358$, slits $[1,R]$ and $[2,R]$ on the six axes plus six harmonically symmetric non-radial arcs; conformal radius from Fedorov's four-point minimal-capacity continuum.
* **Ortega-Cerdà–Pridhnani.** Same axes in a disc of radius $5.1195152501$, two harmonically symmetric arcs per half-sector, capacity computed numerically.

## Searches that found nothing newer

arXiv full text ("univalent Bloch", "Bloch-Landau"), zbMATH, citation graphs of the papers above (Semantic Scholar, OpenAlex), GitHub repositories/code/issues, Zenodo, OpenReview, MathOverflow, math.stackexchange, teorth/optimizationproblems (no pull request touches 57c as of 2026-10-13), Hayman–Lingham problem 5.10 ("No progress on this problem has been reported to us").

Not accessed: Beller–Hummel full text; Brown, *Constructing mappings onto radial slit domains*, Rocky Mountain J. Math. 37 (2007); MathSciNet; Google Scholar.

## Possible venues

* arXiv math.CV: AI tools cannot be authors; significant use must be reported. Recent related preprints carry a "tool disclosure" paragraph.
* *Computational Methods and Function Theory* (Springer): complex analysis and computational methods; AI use to be declared.
* *Experimental Mathematics* and *Complex Variables and Elliptic Equations* (Taylor & Francis): generative-AI use must be acknowledged, naming tool and purpose; code should be available.
* The optimization-constants catalogue: Markdown-only pull request, AI use disclosed, contributor verifies everything; the entry should also record the 2009 and 2010 values.
