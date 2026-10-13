# Lessons: the mathematics behind the result, from the basics

Notes from the walk-through sessions. More formal detail is in `WALKTHROUGH.md`.

## Lesson 1: complex numbers, maps and the problem

**1. A complex number is a point.** $z=x+iy$ is the point $(x,y)$. It has a length $|z|$ (distance from the origin) and an angle (anticlockwise from the positive $x$-axis). Example: $1+i$ has length $\sqrt2$ and angle $45^\circ$.

**2. Multiplication rotates and scales.** The only rule is $i^2=-1$. To multiply two complex numbers, multiply their lengths and add their angles. Example: $i(3+2i)=-2+3i$, a quarter turn. So multiplying by a fixed $c$ rotates by the angle of $c$ and scales by $|c|$.

**3. A function moves points.** $f(z)=z^2$ squares lengths and doubles angles: $(1+i)^2=2i$.

**4. The derivative describes a tiny piece.** $f(z_0+h)\approx f(z_0)+f'(z_0)h$. The small step $h$ is multiplied by the number $f'(z_0)$, so it is rotated by the angle of $f'(z_0)$ and scaled by $|f'(z_0)|$. For $f(z)=z^2$ at $z_0=1+i$: $f'(z_0)=2+2i$, scale $2\sqrt2$, rotation $45^\circ$.

**5. Why angles survive.** Two small arrows at $z_0$ are multiplied by the same number, so the angle between them is unchanged. A function with a complex derivative (non-zero) is *conformal*. If it is also one-to-one it is *univalent*.

**6. Every shape is a stretched disc (Riemann mapping theorem).** If $D$ is a region with no holes, not the whole plane, containing $0$, there is a univalent $f$ from the unit disc onto $D$ with $f(0)=0$, unique up to rotating the disc. The number $|f'(0)|$ therefore belongs to $D$: it is the *conformal radius* $\rho(D)$.
* Disc of radius 5: $f(z)=5z$, $\rho=5$.
* Half-plane whose edge is 1 away from the origin: $f(z)=2z/(1-z)$, $\rho=2$, although the shape is infinite.

**7. The problem.** For a shape $D$ let $R(D)$ be the radius of the largest disc ("coin") inside it. The constant is
$$\mathcal U=\min_D\frac{R(D)}{\rho(D)} .$$
* A *lower* bound is a theorem about every shape.
* An *upper* bound is one shape with its two numbers. That is what we work on.

**8. Why cuts.** A thin cut stops coins from sitting across it, so $R$ drops a lot; but it is a small target as seen from the centre, so $\rho$ drops little.

Worked example: the disc of radius 2 has $R=\rho=2$. Cut it along the real axis from $x=1$ to $x=2$. The coin slides away from the cut: centre $(-\tfrac12,0)$, radius $1.5$. The conformal radius becomes $16/9\approx1.778$. Ratio $0.84$. One cut is not enough; the coin just moves, which is why cuts are needed all round at several depths.

**9. What $\rho$ measures exactly.** Release a random walker (Brownian motion) at $0$ and stop it when it first hits the boundary. Then
$$\log\rho(D)=\text{average of }\log(\text{distance of the landing point from }0).$$
So $\rho$ is the geometric mean of the landing distances: average the logs, then exponentiate. The weights are the landing probabilities (*harmonic measure*), not equal weights by direction.
* Geometric mean of 1 and 4: logs $0$ and $1.386$, average $0.693$, $e^{0.693}=2$.
* Made-up example: land at distance 2 with probability $0.8$ and at distance 1 with probability $0.2$: $\log\rho=0.8\log2=0.555$, $\rho=1.74$.

This explains the strategy: to lower $\rho$ the walker must actually land on a cut, near its inner end. A thin cut is hard to hit there.

### Questions answered so far
* Disc of radius 2: $R=2$, $\rho=2$, ratio 1.
* One cut from 1 to 2: $R=1.5$ (not 1), $\rho=16/9$.

### Open question
Add a second cut opposite the first, from $x=-2$ to $x=-1$. Where does the largest coin go, and roughly what is $R$?

## Next: Lesson 2

Green's functions: the function that encodes "where the walker lands", why its constant term at the pole is $\log\rho$, and how that turns the problem into solving Laplace's equation on the shape.
