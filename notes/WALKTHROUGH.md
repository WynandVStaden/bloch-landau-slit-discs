# Walkthrough: what was computed and how to check it

This is for the author, to be able to explain and re-run every step. It assumes a first course in complex analysis.

## 1. The statement

Let $f$ be univalent on the unit disc $\mathbb D$ with $f(0)=0$. Write $R(D)$ for the inradius of a domain $D$ (the supremum of radii of discs contained in $D$). The univalent Bloch–Landau constant is

$$\mathcal U=\inf_f \frac{R(f(\mathbb D))}{|f'(0)|}.$$

By the Riemann mapping theorem every simply connected proper subdomain $D\ni0$ of $\mathbb C$ is $f(\mathbb D)$ for such an $f$, unique up to rotation, and $|f'(0)|$ is the conformal radius $\rho(D)$ of $D$ at $0$. So

$$\mathcal U\le\frac{R(D)}{\rho(D)}\quad\text{for every simply connected } D\ni0 .$$

An upper bound for $\mathcal U$ is therefore one domain together with two numbers: an upper bound for its inradius and a lower bound for its conformal radius. All published upper bounds (Robinson 1935, Goodman 1945, Beller–Hummel 1985, Carroll–Ortega-Cerdà 2009) are of this form.

## 2. The domains

Take the disc $|z|<R$ and remove finitely many slits, each starting on the circle and ending at a free tip. The complement in the sphere is connected, so the domain is simply connected. The slits stop large discs from fitting, at little cost in conformal radius.

We keep the symmetry of the classical examples: invariance under $z\mapsto e^{2\pi i/3}z$ and $z\mapsto\bar z$. The two reported domains are in `designs.py`:

* **straight**: $R=5.23368$; radial slits with tips at radius $1$ (angles $0,120^\circ,240^\circ$), $2$ ($60^\circ,180^\circ,300^\circ$), $2.83043$ ($\pm18.903^\circ$ and rotations), $3.74887$ ($\pm45.979^\circ$ and rotations). 18 slits.
* **curved**: $R=4.7477\ldots$; the same six axial slits and twelve curved ones (parameters in `designs.py`).

A curved slit is described in the coordinate $w=\log z$ by

$$w(x)=\log R+a x+i\big(\alpha+\beta (a x)^2\big),\qquad -1\le x\le0,\quad a=\log(R/t),$$

where $\alpha$ is the angle at which it leaves the circle, $t$ the modulus of the tip and $\beta$ the bend. $\beta=0$ is a radial slit.

## 3. The conformal radius

Let $g$ be the Green function of $D$ with pole at $0$: harmonic in $D\setminus\{0\}$, zero on $\partial D$, and $g(z)=-\log|z|+\log\rho(D)+o(1)$ at $0$. So $\log\rho$ is the constant term of $g$ at its pole.

**Change of variable.** Put $u(w)=g(e^w)$. With $L=\log R$, $u$ is harmonic in the half-plane $\{\operatorname{Re}w<L\}$ minus the slits, $2\pi i$-periodic, zero on $\operatorname{Re}w=L$ and on the slits, and $u(w)=-\operatorname{Re}w+\log\rho+o(1)$ as $\operatorname{Re}w\to-\infty$.

**Reflection.** Since $u=0$ on the line $\operatorname{Re}w=L$, it extends to the whole plane as an odd function under $w\mapsto 2L-\bar w$. Each slit and its mirror image form a single arc $\Gamma_j$, parametrised by $x\in[-1,1]$ with the formula above. Because the slit meets the circle at a right angle and the imaginary part is even in $x$, $\Gamma_j$ is analytic: the junction with the circle has disappeared, and the only singular points are the two tips.

**Representation.** With $\Phi(w)=-\log|2\sinh(w/2)|$, which is $2\pi i$-periodic, harmonic off $2\pi i\mathbb Z$, has a logarithmic pole at $0$ and satisfies $\Phi(w)=-|\operatorname{Re}w|/2+o(1)$ as $|\operatorname{Re}w|\to\infty$, write

$$U(w)=(L-\operatorname{Re}w)+\sum_j\int_{-1}^{1}\Phi\big(w-w_j(y)\big)\,\frac{\psi_j(y)}{\sqrt{1-y^2}}\,dy .$$

The factor $1/\sqrt{1-y^2}$ carries the inverse-square-root singularity of the charge density at a slit tip, so $\psi_j$ is smooth. By symmetry $\psi_j$ is odd, so each arc carries total charge zero.

**The equation.** $U=0$ on $\Gamma_i$ reads

$$\sum_j\int_{-1}^{1}\Phi\big(w_i(x)-w_j(y)\big)\frac{\psi_j(y)}{\sqrt{1-y^2}}\,dy=a_i x,\qquad -1\le x\le1 .$$

**The answer.** As $\operatorname{Re}w\to-\infty$, $\Phi(w-\omega)=\tfrac12(\operatorname{Re}w-\operatorname{Re}\omega)+o(1)$; the terms in $\operatorname{Re}w$ cancel because the total charge is zero, so

$$\log\rho=L-\frac12\sum_j\int_{-1}^{1}\operatorname{Re}w_j(y)\frac{\psi_j(y)}{\sqrt{1-y^2}}dy=L-\frac{\pi}{4}\sum_j a_j\,c_{j,1},$$

where $\psi_j=\sum_n c_{j,n}T_n$ is the Chebyshev expansion (use $\operatorname{Re}w_j=L+a_jy$ and $\int_{-1}^1 yT_n(y)(1-y^2)^{-1/2}dy=\tfrac\pi2\delta_{n1}$).

**Discretisation.** On one arc, $\Phi(w_j(x)-w_j(y))=-\log|x-y|+B_j(x,y)$ with $B_j$ analytic, and

$$\int_{-1}^{1}-\log|x-y|\,\frac{T_n(y)}{\sqrt{1-y^2}}dy=\begin{cases}\pi\log2,&n=0\\ \dfrac{\pi}{n}T_n(x),&n\ge1.\end{cases}$$

So the singular part is diagonal in the Chebyshev basis and everything else is an integral of an analytic function against the Chebyshev weight, for which Gauss–Chebyshev quadrature converges exponentially. `bloch_solver.py` collocates at the $n$ Chebyshev nodes per arc. With 18 slits and $n=96$ this is a dense $1728\times1728$ system with condition number about $700$.

## 4. The inradius

For straight slits (`bloch_design.py`): the distance from a point to the boundary is the minimum of its distances to the circle and to each segment, all given by exact formulas. A polar grid over the sector $0\le\arg z\le60^\circ$ (enough by symmetry) locates the local maxima; because the distance function is 1-Lipschitz, every maximum within the grid spacing of the best one is kept. Each is refined by Newton's method on the three equations "distance to feature $k$ equals $r$".

For curved slits (`bloch_curved.py`): the same, with the nearest point of each slit found by Newton's method on the curve parameter.

## 5. Checks, and how to run them

All from this folder, with the Python environment in `..\constants\.venv`.

| Command | What it shows |
|---|---|
| `python bloch_solver.py` | the solver against two exact formulas (one slit; six slits via $z^3$ and Koebe maps): agreement to $10^{-16}$ |
| `python designs.py` | both ratios recomputed, and the figure |
| `python bloch_verify.py 5.23368 2.0 0.80248 3.74887 0.32992 2.83043` | straight domain: inradius by three routines and 30 million random points; conformal radius at five resolutions and after a rotation |
| `python bloch_residual.py 5.23368 2.0 0.80248 3.74887 0.32992 2.83043` | the computed Green function evaluated on the slits at points not used in the computation |
| `python bloch_verify_curved.py <R t2 then angle tip bend per class>` | the same for a curved domain (slow: about 45 minutes) |

Why the residual check matters: the function $U_h$ produced by the solver is harmonic in $D\setminus\{0\}$, has the right pole and vanishes on the circle, by construction. So $g-U_h$ is harmonic in $D$ and equals $-U_h$ on the slits. By the maximum principle

$$|\log\rho-\log\rho_h|\le\max_{\text{slits}}|U_h| .$$

The measured maximum is below $10^{-12}$.

Earlier records reproduced by the same code: six slits give $0.65779$ (Robinson: $0.658$); twelve straight slits give $0.6564163$ (Beller–Hummel: $0.6564155$).

## 6. What is not proved

Everything above is floating-point computation. A proof needs:

1. *Inradius of the straight domain.* Finitely many pockets, each bounded by segments, tips and a circular arc; exact algebra.
2. *A rigorous bound on $\max|U_h|$ on the slits.* $U_h$ is an explicit function once the Chebyshev coefficients are fixed (as rational numbers). What is needed is a validated evaluation of the integrals defining it and a validated bound of its maximum on each slit, by interval arithmetic. Because the claimed improvement is $3\times10^{-5}$ to $9\times10^{-5}$ in the ratio, a bound of $10^{-6}$ would be ample.

## 7. Questions you should be able to answer

* *Why does one domain give a bound?* Section 1.
* *Why is a slit disc allowed?* It is simply connected; no regularity of the boundary is required. The 2009 domain is also a slit disc.
* *How do you know the inradius is not larger?* Three routines and a brute-force search agree; for the straight domain an independent computation also found 1.
* *How accurate is the conformal radius?* 12 digits by self-consistency; an independent Schwarz–Christoffel computation of the first straight domain agreed to 9 digits.
* *Is it optimal?* No. The domains are local optima of a small family. Carroll and Ortega-Cerdà's necessary condition (harmonic symmetry of free arcs) is not imposed.
* *What did the AI do?* See the disclosure in `NOTE.md`.
