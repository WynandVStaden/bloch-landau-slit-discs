# Primer: the problem and the method, from the beginning

This guide explains what the repository proves and how, assuming only school mathematics and some patience. The precise statements and proofs are in the [paper](../paper/bloch_slit_discs.pdf); this is the intuition behind them.

Formulas are written in plain text so that they display everywhere. `log` is the natural logarithm and `|z|` is the distance of the point z from the centre.

## 1. The game

Think of a game with a player and a referee.

1. **The player draws a shape**: any region of the plane without holes, with a marked centre point.
2. **The referee measures two numbers:**
   - **R**, the radius of the biggest coin (disc) that fits anywhere inside the shape;
   - **ρ** (rho), the size of the shape *as felt from the centre* (explained in section 2).
3. **The score is R / ρ.** Lower is better.

A plain disc of radius 3 has R = 3 and ρ = 3, so it scores 1. That is the worst possible score.

The **univalent Bloch–Landau constant U** is the lowest score any shape can achieve. It is known that U ≥ 0.5708. Nobody knows U exactly. Each explicit shape shows that U is at most that shape's score, and this repository contains a shape with score below 0.65628589, lower than any found before.

Two things follow from this set-up:

- Only the *shape* matters, not its size: enlarging a shape multiplies R and ρ by the same factor.
- To prove an upper bound for U you only need **one** shape and a rigorous measurement of its two numbers. How the shape was found is irrelevant.

## 2. What ρ measures

Release a particle at the centre and let it wander at random until it first touches the edge of the shape. Record how far from the centre it landed. Repeat many times.

**ρ is the geometric mean of those landing distances** (average the logarithms, then exponentiate).

- Plain disc of radius 2: every landing point is at distance 2, so ρ = 2.
- A shape whose edge is partly near and partly far: ρ lies in between, weighted by how often the particle lands on each part.

Two facts about wandering particles drive everything that follows:

- Edge that is **close to the centre** is hit often and pulls ρ down strongly.
- Edge that is **thin or hidden** (a narrow slit, the far end of a corridor) is hit rarely and hardly affects ρ.

ρ is called the *conformal radius*. Its textbook definition uses conformal maps; the wandering-particle description is equivalent.

## 3. Why cuts help

A **cut** (or slit) is a thin line sawn into the shape from the outer edge inward. It counts as edge.

- A cut **blocks coins completely**: no coin may lie across it. So cuts make R small.
- A cut **catches only some particles**, because it is thin. So cuts reduce ρ only a little, provided they stay away from the centre.

Worked example, a disc of radius 3:

| Shape | R | ρ | Score |
|---|---|---|---|
| no cut | 3 | 3 | 1.00 |
| one cut, tip at distance 1.5 from the centre | 2.25 | 2.67 | 0.84 |
| the same cut pushed in to distance 0.5 | 1.75 | 1.47 | 1.19 |
| two opposite cuts, tips at distance 1.73 | 2 | 2.60 | 0.77 |

One cut helps, but pushing it too close to the centre backfires: the particle hits the tip almost at once and ρ collapses. More cuts block the coin from more sides, as long as each keeps its distance. The whole problem is this trade-off, tuned very finely:

| Year | Arrangement | Score |
|---|---|---|
| 1935 | 6 straight cuts at two depths | 0.658 |
| 1945–1985 | more straight cuts | 0.65647, then 0.6564155 |
| 2009 | cuts that fork like a tree | 0.6563937 |
| 2010 | a larger tree (numerical value only) | 0.656319277 |
| here | 24 cuts, 18 of them curved, found by computer search | 0.65628589 |

## 4. Measuring ρ: the heat picture

Simulating wandering particles would be far too slow and inaccurate. There is an equivalent problem that can be solved to twelve digits in a fraction of a second.

Imagine the shape as a metal plate. Hold its whole edge, cuts included, at temperature 0 and place a small heater at the centre. Wait until the temperature settles.

**Settled temperature has one defining property:** the value at every point equals the average of the values on a small circle around it. (Heat flows in from the warmer side and out on the cooler side in equal amounts.) Functions with this property are called *harmonic*. An immediate consequence, the **maximum principle**: the hottest and coldest points are on the edge, never inside.

For a plain disc of radius r the settled temperature is known exactly:

```
temperature at z  =  log(r) − log|z|
```

It is 0 on the edge (where |z| = r) and rises without bound at the heater. It splits into a heater part `−log|z|`, the same for every shape, and a constant `log(r)` that depends on the plate.

For any shape, the temperature near the heater has the same form:

```
temperature at z  ≈  −log|z| + c
```

and the constant c is exactly log ρ. So

```
ρ = e^c
```

A roomy plate keeps more heat near the heater (large c, large ρ). Cuts held at 0 drain heat (smaller c, smaller ρ). **Finding ρ means finding the settled temperature and reading off one constant.**

## 5. The solver

How do we find the settled temperature of a disc with 24 curved cuts? By building it out of simple pieces.

**Idea 1: sources.** A heater at a point p warms a point z by `−log|z − p|`, which depends only on the distance between them. A *cooler* is the same with the opposite sign. Temperatures simply add. A cut held at 0 acts like a line of coolers, and we look for the temperature in the form

```
temperature = (heater at the centre) − (coolers along the cuts, with unknown strengths)
```

These coolers are a mathematical device, not physical objects: we need their strengths only because ρ can be computed from them.

**Idea 2: mirror images.** The outer circle must be at 0. For every source inside the disc, place an imaginary source of opposite sign at the mirror-image point outside the disc. Each pair cancels exactly on the circle, so the circle is at 0 automatically, whatever the strengths are. Call the combined effect of a unit source at p and its mirror image `G(z, p)`.

**Idea 3: the equations.** What remains is to make the temperature 0 on the cuts. Choose check points on the cuts, as many as there are unknown strengths. At each check point,

```
s1·G(check point, p1) + s2·G(check point, p2) + …  =  log(R) − log|check point|
```

where s1, s2, … are the unknown strengths and R is the outer radius. Every G is a known number computed from distances. These are simultaneous linear equations, solved directly by the computer.

**Idea 4: read off the answer.**

```
c = log(R) − s1·log(R/d1) − s2·log(R/d2) − …        ρ = e^c
```

where d1, d2, … are the distances of the coolers from the centre. A cooler near the outer circle costs almost nothing; one near the centre costs a lot. This is the formula behind the table in section 3.

**Refinements that give the accuracy.** The actual program (`bloch_solver.py`):

- works in the coordinate `w = log z`, in which the outer circle becomes a straight line and the mirror image is an ordinary reflection;
- describes the cooling strength along each cut by a smooth formula with 64 coefficients (Chebyshev polynomials with the correct behaviour at the tip of the cut built in), which is far more efficient than separate point coolers.

With 24 cuts and 64 coefficients each, this is a system of 1,536 linear equations. It reproduces the two cases with exact textbook answers to 16 digits.

## 6. Measuring R

R is the largest distance from any point of the shape to the nearest edge. The distances to the outer circle and to each cut are computed from their formulas, not from a picture.

- **During the search:** sample points on a grid, keep those that are further from the edge than their neighbours, then solve for the exact point that is equally far from its three nearest pieces of edge. (A coin that cannot grow touches the edge in three places.)
- **In the proof:** see section 8.

## 7. Finding the shape

A shape is a short list of numbers: the outer radius and, for each family of cuts, the angle at which it leaves the circle, the radius of its tip, and two numbers describing how it bends. The final shape has 14 such numbers.

The search uses **differential evolution**: keep a population of candidate lists; create new candidates by combining and perturbing existing ones; keep a candidate when its score is lower than the one it replaces. The referee for each candidate is the pair of routines above. Hundreds of thousands of shapes are scored.

The search does not need to be rigorous. Only the final shape is proved.

## 8. From computation to proof

The solver's temperature is exactly 0 at the check points, but only *nearly* 0 between them. The leftover on the cuts is at most

```
δ = 0.0000000266
```

**Why the leftover is harmless.** Lower the computed temperature everywhere by δ. It is now at or below 0 on the whole edge, while the true temperature is exactly 0 there. Both have the same heater. A plate whose edge is colder cannot be warmer inside (maximum principle), so the lowered temperature is nowhere above the true one. Comparing them next to the heater gives

```
true c  ≥  computed c − δ
```

The errors do not accumulate over the 24 cuts: the value at the centre is an average of edge values, and an average of numbers of size at most δ has size at most δ.

**What must be proved.** That δ really is a ceiling at every point of every cut, not only at the points sampled. The program `certify/conformal_certificate.py` does this by evaluating the temperature at several thousand points per cut with guaranteed error bars, and by using its smoothness to bound how far it can rise in between. All arithmetic is *ball arithmetic*: each number carries a guaranteed error radius.

**The coin.** The program `certify/inradius_certificate.py` divides the shape into about 5,000 squares. For each square it names one point of the edge and checks, in exact whole-number arithmetic, that all four corners of the square are within 1.000000001 of it. Then every point of the square is, so no coin centred there can be larger. Squares are small where a coin of radius almost 1 just fits, and large elsewhere.

**The result.**

```
R ≤ 1.000000001        ρ ≥ 1.523726208329        score ≤ 0.65628589
```

## 9. Check your understanding

1. A plain disc of radius 4: what are R, ρ and the score?
2. A disc of radius 3 has one cut. Why does the biggest coin move away from the cut, and what stops it growing?
3. Why does extending a cut toward the centre eventually make the score worse?
4. The solver forces the temperature to be 0 at 1,536 points. Is it then 0 everywhere on the cuts?
5. The search uses a fast, non-rigorous coin finder. Why does that not weaken the final theorem?

**Answers.**
1. R = 4, ρ = 4, score 1.
2. A coin may not lie across the cut, so it slides to the open side; it stops when it touches both the outer circle and the tip of the cut.
3. The tip comes close to the centre, wandering particles hit it almost immediately, and ρ falls faster than R.
4. No. Between those points it is off by up to δ, which is why the comparison argument of section 8 is needed.
5. The theorem is about one shape, and that shape is checked rigorously afterwards. The search only chooses which shape to check.
