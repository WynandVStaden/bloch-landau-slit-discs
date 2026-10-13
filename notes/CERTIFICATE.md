# The rigorous certificate

Goal: turn "the computer says the ratio is 0.656285866" into a proof of an inequality
$$\mathcal U\le\frac{R(D)}{\rho(D)}\le\frac{1+\varepsilon}{\rho_{\rm low}} .$$
Two independent statements are proved, each by a program in `certify/`:

1. `inradius_certificate.py`: every disc inside $D$ has radius at most $1+\varepsilon$.
2. `conformal_certificate.py`: the conformal radius of $D$ is at least $\rho_{\rm low}$.

Both use *ball arithmetic* (the Arb library, through python-flint): every number is stored as a midpoint and a guaranteed error radius, and every operation enlarges the radius to cover its own rounding error. Ordinary floating point is used only to make choices (where to subdivide, how many nodes); a bad choice can make a certificate fail but cannot make a false one pass.

The domain is defined by the decimal parameters printed in `NOTE.md`, read as exact rational numbers, with $\pi/3$ exact.

## 1. Inradius

*Fact used:* if $q$ is a boundary point then $\operatorname{dist}(p,\partial D)\le|p-q|$.

By symmetry it is enough to treat $p$ in the sector $0\le\arg p\le60^\circ$, $|p|\le R$. The square $[0,8]^2$ is cut by a quadtree into cells, each of one of three kinds:

* outside the disc: checked as $x_{\min}^2+y_{\min}^2\ge R^2$;
* outside the sector: checked as $y_{\min}^2>3x_{\max}^2$;
* covered: a boundary point $q$ (the "witness") is named, and all four corners of the cell are within $1+\varepsilon$ of $q$. Since $p\mapsto|p-q|$ is convex, every point of the cell is then within $1+\varepsilon$ of the boundary.

The cells' corners are dyadic rationals; the witness is enclosed to $2^{-69}$ by ball arithmetic and the comparisons are done in exact integer arithmetic. The program also checks that the cells tile the square (no cell contains another, areas add up to that of the square).

Result (2026-10-13): $R(D)\le1+10^{-9}$ for G and for B, with 4993 cells for G.

An independent audit of both programs and of the argument below (2026-10-13, scripts in `referee/certificate_audit/`) found no fatal defect; the repairs it asked for (exact box endpoints in the branch-and-bound, a range check on the cells, rounding the printed bound downwards, and the clarifications marked in the text) have been made.

## 2. Conformal radius

### 2.1 The comparison principle

Let $g$ be the Green function of $D$ with pole at $0$: harmonic in $D\setminus\{0\}$, zero on the boundary, $g(z)=-\log|z|+\log\rho(D)+o(1)$.

**Lemma.** Let $v$ be harmonic in $D\setminus\{0\}$, continuous up to the boundary, with $v(z)=-\log|z|+c+o(1)$ at $0$ and $v\le0$ on $\partial D$. Then $\log\rho(D)\ge c$.

*Proof.* $g-v$ is harmonic in $D$ (the singularities cancel) and $\ge0$ on the boundary, hence $\ge0$ in $D$; at $0$ it equals $\log\rho(D)-c$. ∎

So a lower bound for $\rho$ needs no knowledge of the true Green function: any function with the right pole that is $\le0$ on the boundary will do. The computed solution is used as that function.

### 2.2 The function

In the coordinate $w=\log z$, $L=\log R$, slit $j$ and its mirror image in the line $\operatorname{Re}w=L$ form the arc
$$w_j(x)=L+a_jx+i\big(\theta_j+b_{1,j}(a_jx)^2+b_{2,j}(a_jx)^4\big),\qquad-1\le x\le1,$$
the slit itself being $-1\le x\le0$. Let $\psi_j$ be *odd* polynomials with rational Chebyshev coefficients $c_{j,n}$, equal for slits related by the symmetry group, and
$$U(w)=(L-\operatorname{Re}w)+\sum_j\int_{-1}^1\Phi\big(w-w_j(y)\big)\psi_j(y)\frac{dy}{\sqrt{1-y^2}},\qquad\Phi(w)=-\log\big|2\sinh\tfrac w2\big| .$$
The following hold exactly, whatever the coefficients are:

* $U$ is $2\pi i$-periodic, so it is a function of $z$; it is harmonic in $D\setminus\{0\}$ (the mirror halves of the arcs lie outside the disc) and continuous, including at the slit tips (the logarithmic potential of $T_n(y)/\sqrt{1-y^2}$ is continuous across the endpoints).
* $U=0$ on $|z|=R$: the reflection $w\mapsto2L-\bar w$ maps $w_j(x)$ to $w_j(-x)$, and $\psi_j$ is odd, so $U$ is odd under the reflection.
* As $\operatorname{Re}w\to-\infty$, $\Phi(w-w')=\tfrac12\operatorname{Re}(w-w')+o(1)$, and $\int\psi_j\,dy/\sqrt{1-y^2}=0$, $\int y\,\psi_j\,dy/\sqrt{1-y^2}=\tfrac\pi2c_{j,1}$, so
$$U=-\log|z|+c+o(1),\qquad c=L-\frac\pi4\sum_ja_jc_{j,1}.$$

Hence, by the Lemma applied to $v=U-\delta$:
$$\boxed{\ \log\rho(D)\ \ge\ L-\frac\pi4\sum_ja_jc_{j,1}-\delta\qquad\text{whenever }\delta\ge0\text{ and }U\le\delta\text{ on every slit.}\ }$$
The coefficients come from the floating-point solver, rounded to doubles; if they were poor the only consequence would be a large $\delta$.

Note that the proof does not need the slits to be disjoint: $D$ can be read as the component of $0$ in the disc minus the slits, which is simply connected in any case because every slit is attached to the circle. (An independent check found the slits pairwise disjoint, at distance at least $0.1484$.)

### 2.3 Bounding $U$ on a slit

By symmetry one slit per class is enough. Put $f(x)=U(w_i(x))$ for $-1\le x\le1$; it is odd.

**Evaluating $f$ at a point.** The self-interaction splits as
$$\Phi\big(w_i(x)-w_i(y)\big)=-\log|x-y|-\log|K_i(x,y)|,\qquad K_i(x,y)=\frac{2\sinh\frac{w_i(x)-w_i(y)}2}{x-y},$$
and $\int_{-1}^1-\log|x-y|\,T_n(y)\,dy/\sqrt{1-y^2}=\frac\pi nT_n(x)$ exactly ($n\ge1$). All remaining integrals have the form $\int\operatorname{Re}h(y)\,\psi(y)\,dy/\sqrt{1-y^2}$ with $h=\log(\text{kernel})$ analytic near $[-1,1]$, and are computed by the $N$-point Gauss–Chebyshev rule with the remainder bound below.

**Lemma (coefficients from the real part).** Let $h$ be analytic in a neighbourhood of the closed Bernstein ellipse $E_r$ (foci $\pm1$, semi-axes $\frac12(r\pm r^{-1})$) and $|\operatorname{Re}h|\le A$ on its boundary. Then the Chebyshev coefficients $\alpha_k$ of the real function $\operatorname{Re}h$ on $[-1,1]$ satisfy $|\alpha_k|\le4A/(r^k+r^{-k})$.

*Proof.* With $x=\frac12(u+u^{-1})$, $h$ becomes a Laurent series $\sum b_ku^k$ with $b_{-k}=b_k$, and $\alpha_k=2\operatorname{Re}b_k$. The $k$-th Fourier coefficient of $\operatorname{Re}h$ on $|u|=r$ is $\frac12(b_kr^k+\bar b_kr^{-k})$, of modulus at most $A$, and its real part is $\frac12\operatorname{Re}b_k\,(r^k+r^{-k})$. ∎

Only $|\operatorname{Re}h|=\big|\log|\text{kernel}|\big|$ is needed, so it suffices to bound the modulus of the kernel above and below (away from $0$) on a rectangle containing the ellipse. This is done by branch-and-bound in ball arithmetic; the lower bound also proves that the kernel does not vanish on the simply connected set (rectangle) $\times[-1,1]$, so it has a single-valued logarithm there, continuous in both variables and analytic in the complex one.

**Quadrature remainder.** The rule is exact for polynomials of degree $<2N$. Splitting $\operatorname{Re}h$ at degree $K=2N-1-\deg\psi$,
$$\Big|\int-\text{rule}\Big|\le2\pi\,\|\psi\|_\infty\sum_{k>K}|\alpha_k|\le\frac{8\pi A\|\psi\|_\infty}{r-1}\,r^{-K}.$$

**From points to the whole slit.** Let $p$ interpolate $f$ at the $M$ zeros of $T_M$. Then
$$\sup_{[-1,1]}|f|\le\sum_q|\text{Chebyshev coefficients of }p|+2\sum_{k\ge M}|a_k|,$$
where $a_k$ are the Chebyshev coefficients of $f$. By the Lemma again, now in the variable $x$, $|a_k|\le4B\rho^{-k}$, where $B$ bounds, on an ellipse $E_\rho$ in the complex $x$-plane, the sum of the explicit polynomial terms and of $\pi\|\psi_j\|_\infty\sup\big|\log|\text{kernel}_{ij}|\big|$. So the second term is at most $8B\rho^{1-M}/(\rho-1)$.

The size of the ellipses is limited by the nearest neighbouring slit (the kernel vanishes when the complexified arc reaches another arc). For G the closest slits are $0.03$ apart in the $w$-plane, which gives $\rho\approx1.01$ and a few thousand nodes on the longest slits.

### 2.4 Tests

The program is first run on the two cases with closed-form answers:

| domain | exact | certified lower bound | gap |
|---|---|---|---|
| unit disc minus $[\tfrac12,1]$ | $8/9=0.888888888889$ | $0.888888888328$ | $5.6\times10^{-10}$ |
| radius 3, slits $[1,3]$ and $[2,3]$ on alternate axes | $1.520743036313$ | $1.520743034966$ | $1.3\times10^{-9}$ |

A certified lower bound must never exceed the exact value; both pass.

## 3. What a referee would still want

* The Arb library and python-flint are trusted, as is the Python integer arithmetic.
* The two lemmas above are standard but should be written out with references (Trefethen, *Approximation Theory and Approximation Practice*, ch. 8, for the classical version with $|h|$ in place of $|\operatorname{Re}h|$).
* The programs should be read by a second person; the argument is short enough to check by hand, and the tests in 2.4 guard against sign and normalisation errors.
* The comparison with $0.656319277$ is with a value that is itself not rigorous; the comparison with the proved bound $0.6563937$ is the one that makes a theorem.
