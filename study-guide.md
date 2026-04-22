# Continuous Mathematics — Numerical Methods Study Guide

Covers **Integration**, **Root-Finding**, **Optimization**. For each method: formula, applicable dimension, when to use, conditions, error/convergence. Required proofs vs bonus/non-syllabus material are listed at the end.

---

## Foundational Background (Chapters 1–2, 4)

### Derivatives (nD)

- **Gradient** (column vector of first partial derivatives — points in the direction of steepest increase):

$$\nabla f = \left(\frac{\partial f}{\partial x_1}, \ldots, \frac{\partial f}{\partial x_d}\right)^{T}$$

- **Hessian** — the $d \times d$ matrix of second partial derivatives. It is symmetric when second partials are continuous (Schwarz's theorem):

$$H(f)\_{ij} = \frac{\partial^{2} f}{\partial x_i \partial x_j}$$

- **Jacobian** of $f : \mathbb{R}^n \to \mathbb{R}^m$ — the $m \times n$ matrix of first partial derivatives of each output with respect to each input:

$$J\_{ij} = \frac{\partial f_i}{\partial x_j}$$

Note: The Hessian of a scalar function $f$ equals the Jacobian of $\nabla f$. This is why Newton's method for optimization (finding zeros of $\nabla f$) uses the Hessian as its "Jacobian".

### Taylor's Theorem

Taylor's theorem approximates a function near a point using its derivatives. The **remainder term** quantifies the approximation error.

**1D, order $n$ (Lagrange remainder):**

$$f(x+h) = \sum_{k=0}^{n} \frac{h^{k}}{k!} f^{(k)}(x) + \frac{h^{n+1}}{(n+1)!} f^{(n+1)}(\xi)$$

for some $\xi$ between $x$ and $x+h$. The last term is the remainder — it is **exact**, not an approximation, but $\xi$ is unknown.

**nD, second order:**

$$f(\mathbf{x}+\mathbf{h}) = f(\mathbf{x}) + \mathbf{h}^{T}\nabla f(\mathbf{x}) + \tfrac{1}{2} \mathbf{h}^{T} H(f)(\mathbf{x}+\theta\mathbf{h}) \mathbf{h}$$

for some $\theta \in (0,1)$. The three terms represent: the function value, the linear (gradient) correction, and the quadratic (curvature) correction. This expansion is fundamental to deriving Newton's method and analysing convergence.

### Convexity

A function $f$ is **convex** on a convex set $S$ if, for all $\mathbf{x}, \mathbf{y} \in S$ and $\lambda \in [0,1]$:

$$f(\lambda\mathbf{x} + (1-\lambda)\mathbf{y}) \leq \lambda f(\mathbf{x}) + (1-\lambda)f(\mathbf{y})$$

Geometrically, this means the line segment between any two points on the graph of $f$ lies on or above the graph itself. Equivalently:

- Twice-differentiable $f$ is convex $\iff$ the Hessian is positive semidefinite everywhere.
- Strictly convex $\Leftarrow$ Hessian is positive definite (the converse is not true in general).
- **Key property:** for a convex function, every local minimum is automatically a global minimum, and any stationary point ($\nabla f = \mathbf{0}$) is a global minimum.

### Jensen's Inequality

For convex $f$ and probability weights $\lambda_i \geq 0$ with $\sum \lambda_i = 1$:

$$f\left(\sum_i \lambda_i \mathbf{x}_i\right) \leq \sum_i \lambda_i f(\mathbf{x}_i)$$

In words: the function of a weighted average is at most the weighted average of the function values. This generalises the definition of convexity from two points to any finite collection.

### Stationary Point Classification (nD)

At a stationary point $\mathbf{x}^{\ast}$ (where $\nabla f(\mathbf{x}^{\ast}) = \mathbf{0}$), the nature of the extremum is determined by the Hessian $H(f)(\mathbf{x}^{\ast})$:

| Hessian property | Conclusion |
|---|---|
| Positive definite (all eigenvalues $> 0$) | Strict local **minimum** |
| Negative definite (all eigenvalues $< 0$) | Strict local **maximum** |
| Indefinite (eigenvalues of both signs) | **Saddle point** |
| Semidefinite (some eigenvalue $= 0$) | **Inconclusive** — need higher-order terms |

**Why this works:** From the second-order Taylor expansion, $f(\mathbf{x}^{\ast} + \mathbf{h}) - f(\mathbf{x}^{\ast}) \approx \tfrac{1}{2} \mathbf{h}^T H \mathbf{h}$, so the sign of $\mathbf{h}^T H \mathbf{h}$ for all directions $\mathbf{h}$ determines whether $f$ increases or decreases around $\mathbf{x}^{\ast}$.

**Tests for positive definiteness:** (i) all eigenvalues positive, (ii) all leading principal minors positive, or (iii) Cholesky factorisation succeeds without breakdown.

### Accuracy (Chapter 4)

- **Error**: $\tilde{u} - u$. **Absolute error**: $\lvert \tilde{u} - u \rvert$. **Relative error**: $\lvert \tilde{u} - u \rvert / \lvert u \rvert$.
- For vectors, replace $\lvert\cdot\rvert$ by the Euclidean norm $\lVert\cdot\rVert$.
- **Truncation error**: error from approximating an infinite procedure by a finite one (e.g. Taylor cutoff, $n$-strip integration). This is inherent to the mathematical method.
- **Roundoff error**: error from finite floating-point representation. This is inherent to the computer hardware. *Not* the same as truncation error.
- **Machine epsilon**: $\varepsilon \approx 1.1 \times 10^{-16}$ (double precision), $\approx 6 \times 10^{-8}$ (single precision). It is the smallest $\varepsilon$ such that $1 + \varepsilon \neq 1$ in floating point — equivalently, the maximum relative rounding error when storing a real number.
- **Catastrophic cancellation**: subtracting two nearly equal floating-point numbers can produce a result with huge relative error, because the leading significant digits cancel and you are left with the (inaccurate) trailing digits. Often avoidable by algebraic reformulation, e.g. rewrite $1-\cos x$ as $2\sin^{2}(x/2)$.
- **Forward error**: $\tilde{f}(x) - f(x)$ — how far the computed answer is from the true answer. **Backward error**: the smallest perturbation $\delta x$ such that $\tilde{f}(x) = f(x + \delta x)$ — how much the input would need to change to make the computed answer exact.

### Rates of Convergence (Ch 4)

Let $\epsilon_n = x_n - x^{\ast}$ be the error at step $n$. As $n \to \infty$:

| Type | Condition | Intuition |
|---|---|---|
| Linear (order 1) | $\lvert \epsilon_{n+1}\rvert / \lvert \epsilon_n \rvert \to a$, with $0 < a < 1$ | Error shrinks by a constant factor each step |
| Sublinear | $\lvert \epsilon_{n+1}\rvert / \lvert \epsilon_n \rvert \to 1$ | Convergence slows down over time |
| Logarithmic | Sublinear, plus $\lvert \epsilon_{n+2}-\epsilon_{n+1}\rvert / \lvert \epsilon_{n+1}-\epsilon_n\rvert \to 1$ | Error steps shrink very slowly |
| Superlinear | $\lvert \epsilon_{n+1}\rvert / \lvert \epsilon_n \rvert \to 0$ | Faster than any linear rate |
| Order $q>1$ | $\lvert \epsilon_{n+1}\rvert / \lvert \epsilon_n\rvert^{q} \to a$, with $a > 0$ | Error raised to power $q$ each step |
| Quadratic | Order $q = 2$ | Correct decimal digits roughly double each step |

**Lemma 4.2 (REQUIRED):** If the sequence $A_n$ converges with order $q$, then the subsequence $A_{2n}$ (every other term) converges with order $q^{2}$. Practical use: if you compare two methods — one with order $q$ costing 1 unit per step, vs. one with order $p$ costing 2 units per step — the second method is worth it only if $p > q^2$.

**Rule of thumb:** Linear convergence adds a constant number of correct decimal digits per step. Quadratic convergence doubles the number of correct digits per step (once the error is small enough).

### Lagrange Multipliers (equality constraints)

To minimise $f(\mathbf{x})$ subject to the constraints $g_i(\mathbf{x}) = 0$ for $i = 1, \ldots, m$:

1. Form the **Lagrangian** by introducing one multiplier $\lambda_i$ per constraint:

$$L(\mathbf{x}, \boldsymbol{\lambda}) = f(\mathbf{x}) - \sum_{i=1}^{m} \lambda_i g_i(\mathbf{x})$$

2. Solve the system $\nabla_{\mathbf{x}} L = \mathbf{0}$ (stationarity) together with $\nabla_{\boldsymbol{\lambda}} L = \mathbf{0}$ (which simply reproduces the original constraints $g_i = 0$).

**Geometric intuition:** at a constrained optimum, the gradient of $f$ must be a linear combination of the constraint gradients $\nabla g_i$, because otherwise you could move along the constraint surface and still decrease $f$. The multipliers $\lambda_i$ are the coefficients of this linear combination.

**KKT conditions (for inequality constraints) are NOT on the syllabus this year.**

---

## Part I — Numerical Integration (Chapter 3)

**Goal:** approximate $\displaystyle \int_a^b f(x) dx$.

Most classical rules are **Newton–Cotes**: replace $f$ by a polynomial interpolant on equally-spaced nodes, then integrate the polynomial exactly.

### 1D Rules (single strip of width $h = b-a$)

**Midpoint rule.** One node at $m = (a+b)/2$. Approximate the area by a rectangle of height $f(m)$:

$$\int_a^b f(x) dx \approx (b-a) f(m)$$

Error: $\displaystyle -\frac{(b-a)^{3}}{24} f''(\xi)$ for some $\xi \in (a,b)$. Exact for polynomials of degree $\leq 1$.

**Trapezium rule.** Two nodes at $a, b$. Approximate the area by a trapezoid:

$$\int_a^b f(x) dx \approx \frac{b-a}{2}\bigl(f(a) + f(b)\bigr)$$

Error: $\displaystyle \frac{(b-a)^{3}}{12} f''(\xi)$. Exact for polynomials of degree $\leq 1$.

**Simpson's rule.** Three nodes at $a, m, b$. Fit a quadratic through the three points and integrate it exactly:

$$\int_a^b f(x) dx \approx \frac{b-a}{6}\bigl(f(a) + 4f(m) + f(b)\bigr)$$

Error: $\displaystyle -\frac{(b-a)^{5}}{2880} f^{(4)}(\xi)$. Exact for polynomials of degree $\leq 3$ (one degree better than the expected $\leq 2$, thanks to the symmetry of the nodes — the cubic error term vanishes).

**Boole's rule.** Five equally-spaced nodes $f_0, f_1, f_2, f_3, f_4$:

$$\int_a^b f(x) dx \approx \frac{b-a}{90}\bigl(7 f_0 + 32 f_1 + 12 f_2 + 32 f_3 + 7 f_4\bigr)$$

Error: $O\bigl((b-a)^{7} f^{(6)}\bigr)$. Exact for polynomials of degree $\leq 5$.

### Composite Rules ($n$ strips of width $h = (b-a)/n$)

Sum the single-strip rule over $n$ sub-intervals $[x_{i-1}, x_i]$. The error improves because $h$ is now $(b-a)/n$ instead of $b-a$.

**Composite midpoint:**

$$\int_a^b f(x) dx \approx h\sum_{i=1}^{n} f(m_i), \qquad \lvert E \rvert \leq \frac{(b-a)^{3}}{24 n^{2}} \max_{x \in [a,b]} \lvert f''(x) \rvert$$

where $m_i$ is the midpoint of the $i$-th sub-interval.

**Composite trapezium:**

$$\int_a^b f(x) dx \approx h\left(\tfrac{1}{2}f(a) + \sum_{i=1}^{n-1} f(x_i) + \tfrac{1}{2}f(b)\right), \qquad \lvert E \rvert \leq \frac{(b-a)^{3}}{12 n^{2}} \max \lvert f'' \rvert$$

Note the endpoints are weighted by $1/2$. Interior nodes shared between adjacent strips are counted once.

**Composite Simpson's** (requires $n$ even, i.e. an even number of strips):

$$\int_a^b f(x) dx \approx \frac{h}{3}\left( f_0 + 4(f_1 + f_3 + \cdots + f_{n-1}) + 2(f_2 + f_4 + \cdots + f_{n-2}) + f_n \right)$$

$$\lvert E \rvert \leq \frac{(b-a)^{5}}{180 n^{4}} \max \lvert f^{(4)} \rvert$$

The alternating weights 1–4–2–4–2–…–4–1 arise from applying Simpson's rule to each pair of strips.

**Convergence rates:** Midpoint and trapezium are $O(n^{-2})$; Simpson's is $O(n^{-4})$; Boole's is $O(n^{-6})$. Doubling $n$ reduces midpoint/trapezium error by $\approx 4\times$ and Simpson's error by $\approx 16\times$.

### When to Use Which

- $f$ only $C^{2}$ smooth $\to$ midpoint or trapezium (higher-order rules won't achieve their theoretical rate because $f^{(4)}$ doesn't exist).
- $f$ is $C^{4}$ smooth $\to$ Simpson's (the default workhorse — best effort-to-accuracy tradeoff for smooth functions).
- Very smooth $f$ $\to$ Boole's or higher-order rules (caveat: high-order Newton–Cotes rules develop negative weights and become numerically unstable; Gaussian quadrature is preferred but is not in this syllabus).
- Endpoint singularities $\to$ midpoint (avoids evaluating at endpoints where $f$ may blow up).
- Non-smooth regions (e.g. kinks) $\to$ split the interval at the non-smooth point, or the convergence rate will degrade to $O(n^{-1})$ or worse.

### Comparison

| Feature | Midpoint | Trapezium | Simpson's |
|---|---|---|---|
| $f$ evaluations (composite, shared) | $n$ | $n+1$ | $2n+1$ |
| Convergence rate | $O(n^{-2})$ | $O(n^{-2})$ | $O(n^{-4})$ |
| Evaluates $f$ at endpoints? | No | Yes | Yes |
| Exact for polynomial degree | $\leq 1$ | $\leq 1$ | $\leq 3$ |

**Midpoint vs Trapezium:** The midpoint error is roughly half the magnitude of the trapezium error and opposite in sign. This is why Simpson's rule — which is the weighted combination $\frac{2 \cdot \text{midpoint} + \text{trapezium}}{3}$ — achieves a much higher-order cancellation.

### Multi-Dimensional Integration

- **Tensor-product Newton–Cotes** (e.g. Simpson's in $d$ dimensions) costs $O(n^{d})$ evaluations: this is the **curse of dimensionality**. The convergence rate per axis is unchanged (e.g. $O(n^{-4})$ for Simpson's), but the total number of function evaluations grows exponentially in $d$.
- For large $d$ or irregular regions, switch to **Monte Carlo** integration.

### Monte Carlo Integration

Draw $X_1, \ldots, X_N$ i.i.d. uniformly on region $R \subset \mathbb{R}^{d}$ with volume $A(R)$. The estimator is:

$$\mathrm{MC}\_{N} = A(R) \cdot \frac{1}{N} \sum_{i=1}^{N} f(X_i)$$

This works because $\mathbb{E}[f(X)] = \frac{1}{A(R)}\int_{R} f d\mathbf{x}$, so multiplying by $A(R)$ recovers the integral.

**Unbiased:**

$$\mathbb{E}[\mathrm{MC}\_{N}] = \int_{R} f d\mathbf{x}$$

**Variance:**

$$\mathrm{Var}(\mathrm{MC}\_{N}) = \frac{A(R)^{2} \mathrm{Var}(f(X))}{N}$$

The variance decreases as $1/N$, regardless of the dimension $d$.

**RMS error:** $O(N^{-1/2})$ — **independent of dimension $d$**. This is the key advantage of MC: for $d \geq 5$, MC typically beats deterministic quadrature despite its slow convergence rate.

**Chebyshev's inequality** gives a distribution-free confidence bound:

$$\mathbb{P}\left(\lvert \mathrm{MC}\_{N} - I \rvert \geq \frac{k \sigma}{\sqrt{N}}\right) \leq \frac{1}{k^{2}}$$

where $\sigma^2 = A(R)^2 \mathrm{Var}(f(X))$ and $I$ is the true integral. For example, with $k = 10$ the probability of exceeding $10\sigma/\sqrt{N}$ error is at most $1\%$.

**Central Limit Theorem:** for large $N$, $\mathrm{MC}\_{N}$ is approximately normally distributed, giving tighter Gaussian confidence intervals (e.g. $\pm 1.96\sigma/\sqrt{N}$ for 95% confidence).

#### Variance Reduction

These techniques reduce the variance $\mathrm{Var}(f(X))$, making MC converge faster without increasing $N$:

- **Stratified sampling**: partition $R$ into sub-regions, sample independently inside each. By sampling each sub-region separately, you eliminate the between-region component of variance. Variance is never worse than plain MC.
- **Importance sampling**: instead of sampling uniformly, draw from a density $p(\mathbf{x})$ chosen to be roughly proportional to $\lvert f(\mathbf{x})\rvert$. Weight each sample by $f(X)/p(X)$ to keep the estimator unbiased. Huge gains when $p$ closely matches the shape of $\lvert f \rvert$.
- **Control variates**: find a function $\tilde{f}$ whose integral is known analytically. Estimate $\int (f - \tilde{f}) d\mathbf{x}$ by MC and add the known $\int \tilde{f} d\mathbf{x}$. This works well when $f - \tilde{f}$ has much smaller variance than $f$ alone.

### When to Use Monte Carlo

- $d \geq 5$ approximately (tensor grids become infeasible).
- Irregular or complicated region $R$ (hard to set up tensor-product grids).
- $f$ only available as a simulation output (black-box function).
- Low-to-moderate accuracy is acceptable ($\sqrt{N}$ convergence is slow — to halve the error, you need $4\times$ as many samples).

---

## Part II — Root-Finding (Chapter 5)

**Goal:** find $\mathbf{x}^{\ast}$ with $\mathbf{f}(\mathbf{x}^{\ast}) = \mathbf{0}$.

### 1D — Bracketing Method

#### Interval Bisection

Requires a **bracket**: an interval $[a, b]$ with $f(a) f(b) < 0$ and $f$ continuous on $[a,b]$. By the Intermediate Value Theorem, there is at least one root inside.

**Algorithm:** Set $m = (a+b)/2$. If $f(a)f(m) < 0$, set $b = m$; otherwise set $a = m$. Repeat.

- **Convergence:** linear, with rate $1/2$: $\lvert \epsilon_{n+1}\rvert / \lvert \epsilon_n \rvert \leq 1/2$ always.
- **Error bound:** $\lvert \epsilon_n \rvert \leq (b_0 - a_0)/2^{n+1}$ after $n$ iterations.
- **When to use:** guaranteed convergence; no derivatives needed; very robust. The downside is slow convergence — it gains only about 1 binary digit of accuracy per step.

### 1D — Open Methods (no bracket required)

These methods converge faster but have no convergence guarantee — they can diverge if the initial guess is poor.

#### Newton's Method (1D)

$$x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}$$

**Derivation:** Expand $f$ in a first-order Taylor series around $x_n$: $f(x) \approx f(x_n) + f'(x_n)(x - x_n)$. Setting this linear approximation to zero and solving for $x$ gives the formula above. Geometrically, $x_{n+1}$ is where the tangent line at $(x_n, f(x_n))$ crosses the $x$-axis.

- **Convergence:** quadratic if $f'(x^{\ast}) \neq 0$ (simple root), $f''$ is continuous near $x^{\ast}$, and $x_0$ is sufficiently close to $x^{\ast}$.
- At a **multiple root** ($f'(x^{\ast}) = 0$), convergence degrades to linear. Remedy: use $x_{n+1} = x_n - m \cdot f(x_n)/f'(x_n)$ if the multiplicity $m$ is known.
- **Pitfalls:** can diverge, cycle, or overshoot (e.g. if $f'(x_n) \approx 0$ at an iterate); requires computing $f'$.

#### Secant Method

Replace $f'(x_n)$ with a finite-difference approximation using the two most recent iterates:

$$x_{n+1} = x_n - f(x_n) \cdot \frac{x_n - x_{n-1}}{f(x_n) - f(x_{n-1})}$$

Geometrically, draw the secant line through $(x_{n-1}, f(x_{n-1}))$ and $(x_n, f(x_n))$, and take $x_{n+1}$ as its $x$-intercept.

- **Convergence:** superlinear, order $\varphi = (1+\sqrt{5})/2 \approx 1.618$ (the golden ratio).
- Only **one** new $f$ evaluation per step (cache the value from the previous step), compared to one $f$ + one $f'$ for Newton.
- **Pitfall:** the denominator $f(x_n) - f(x_{n-1})$ becomes tiny near a root $\to$ catastrophic cancellation; no bracket guarantee.

#### Halley's Method

$$x_{n+1} = x_n - \frac{2 f(x_n) f'(x_n)}{2 f'(x_n)^{2} - f(x_n) f''(x_n)}$$

- **Convergence:** cubic (order 3) — correct digits roughly triple each step.
- **Pitfall:** needs both $f'$ and $f''$; expensive per step, so only worthwhile when $f', f''$ are cheap.

#### Muller's Method

Fit a **quadratic** through the three most recent points $(x_{n-2}, f(x_{n-2})), (x_{n-1}, f(x_{n-1})), (x_n, f(x_n))$. Take $x_{n+1}$ = root of that quadratic closest to $x_n$.

- **Convergence:** superlinear, order $\approx 1.839$.
- **Special feature:** can find **complex roots** even starting from real iterates (because the quadratic formula can produce complex values).

#### Inverse Quadratic Interpolation (IQI)

Fit a quadratic to the points $(f(x_i), x_i)$ — swap the roles of $x$ and $y$ — and evaluate at $y = 0$ to get $x_{n+1}$.

- **Order:** $\approx 1.839$ (same as Muller).
- **Pitfall:** fails if any two of the three $f$-values coincide (the interpolating polynomial is not unique).

#### Brent's Method

A hybrid of bisection, secant, and IQI. It maintains a bracket at all times for safety. At each step, it attempts a superlinear step (secant or IQI); if that step falls outside the current bracket or is too large, it falls back to bisection.

- **Convergence:** superlinear in typical cases; guaranteed linear convergence in the worst case (because of the bracket).
- **Use:** the **robust default** 1D root-finder. The classic `zeroin` implementation guarantees a solution to within $4\varepsilon$ relative tolerance.

### Comparison (1D)

| Method | Order | Needs derivatives | Bracket | Guaranteed? |
|---|---|---|---|---|
| Bisection | 1 (linear) | No | Yes | Yes |
| Newton | 2 | $f'$ | No | No |
| Secant | $\approx 1.618$ | No | No | No |
| Halley | 3 | $f', f''$ | No | No |
| Muller / IQI | $\approx 1.839$ | No | No | No |
| Brent | $\geq 1$ (typically $\approx 1.839$) | No | Yes | Yes |

### Error Analysis — Newton's Method (Lemma 5.3, REQUIRED)

**Setup:** Let $x^{\ast}$ be a simple root ($f(x^{\ast}) = 0$, $f'(x^{\ast}) \neq 0$) with $f''$ continuous near $x^{\ast}$.

**Proof sketch:** By Taylor expansion of $f(x^{\ast})$ around $x_n$:

$$f(x^{\ast}) = f(x_n) + f'(x_n)(x^{\ast} - x_n) + \tfrac{1}{2} f''(\xi_n)(x^{\ast} - x_n)^{2} = 0$$

Dividing by $f'(x_n)$ and rearranging, using $\epsilon_n = x_n - x^{\ast}$:

$$\epsilon_{n+1} = -\frac{f''(\xi_n)}{2 f'(x_n)} \epsilon_n^{2}$$

for some $\xi_n$ between $x_n$ and $x^{\ast}$. As $n \to \infty$, $x_n \to x^{\ast}$, so $\xi_n \to x^{\ast}$ and $f'(x_n) \to f'(x^{\ast})$:

$$\frac{\lvert \epsilon_{n+1} \rvert}{\lvert \epsilon_n \rvert^{2}} \to \frac{\lvert f''(x^{\ast}) \rvert}{2 \lvert f'(x^{\ast}) \rvert}$$

This establishes **quadratic convergence** with asymptotic error constant $\lvert f''(x^{\ast}) \rvert / (2\lvert f'(x^{\ast}) \rvert)$.

### A Posteriori Error Estimate

Given a candidate root $\tilde{x}$, estimate the distance to the true root using a first-order Taylor argument:

$$\lvert \tilde{x} - x^{\ast} \rvert \approx \frac{\lvert f(\tilde{x}) \rvert}{\lvert f'(\tilde{x}) \rvert}$$

This is exact when $f$ is linear. For nonlinear $f$ it is a good estimate when $\tilde{x}$ is close to $x^{\ast}$. It's essentially one Newton step's size — if the correction is small, you're close to the root.

### Multi-Dimensional Root-Finding

$\mathbf{f} : \mathbb{R}^{d} \to \mathbb{R}^{d}$, find $\mathbf{x}^{\ast}$ with $\mathbf{f}(\mathbf{x}^{\ast}) = \mathbf{0}$. Note: need as many equations as unknowns ($d$ equations for $d$ unknowns).

#### Newton's Method (nD)

Linearise $\mathbf{f}$ at $\mathbf{x}_n$: $\mathbf{f}(\mathbf{x}) \approx \mathbf{f}(\mathbf{x}_n) + J(\mathbf{x}_n)(\mathbf{x} - \mathbf{x}_n)$, where $J$ is the Jacobian.

Setting the linear approximation to zero and solving:

$$J(\mathbf{x}_n) \Delta\mathbf{x} = -\mathbf{f}(\mathbf{x}_n), \qquad \mathbf{x}_{n+1} = \mathbf{x}_n + \Delta\mathbf{x}$$

Do **not** form $J^{-1}$ explicitly; solve the linear system directly (LU, Cholesky, etc.) — this is both more numerically stable and more efficient.

- **Cost:** $O(d^{2})$ partial-derivative evaluations (to fill the Jacobian) plus $O(d^{3})$ to solve the linear system per step.
- **Convergence:** quadratic if $J(\mathbf{x}^{\ast})$ is nonsingular and $\mathbf{x}_0$ is sufficiently close.
- **Fails** if $J(\mathbf{x}_n)$ is singular or nearly singular.

#### Broyden's Method (Quasi-Newton)

Avoid computing the full $d \times d$ Jacobian at each step by building up an approximation $\hat{J}_n$ via rank-1 updates.

Let $\Delta\mathbf{x} = \mathbf{x}_n - \mathbf{x}\_{n-1}$ and $\Delta\mathbf{f} = \mathbf{f}(\mathbf{x}_n) - \mathbf{f}(\mathbf{x}\_{n-1})$:

$$\hat{J}\_n = \hat{J}\_{n-1} + \frac{\Delta\mathbf{f} - \hat{J}\_{n-1} \Delta\mathbf{x}}{\lVert \Delta\mathbf{x} \rVert^{2}} \Delta\mathbf{x}^{T}$$

This is the **minimal change** to $\hat{J}\_{n-1}$ that satisfies the **secant equation** $\hat{J}\_n \Delta\mathbf{x} = \Delta\mathbf{f}$ (the multi-dimensional analogue of the secant condition).

- **Convergence:** superlinear (not quadratic), but much cheaper per step than Newton.
- Use the **Sherman–Morrison formula** to update $\hat{J}^{-1}$ directly in $O(d^{2})$ — avoids the $O(d^{3})$ linear system solve each step.
- **Initialisation:** $\hat{J}_0 = I$ (identity — first step is then just $-\mathbf{f}$), or use the true $J(\mathbf{x}_0)$ if it is cheap to compute.

### Termination Criteria (Root-Finding)

Stop iterating when any of the following hold:

1. **Step is small:** $\lvert x_n - x_{n-1} \rvert < \text{tol} \cdot (1 + \lvert x_n \rvert)$ — the relative change in $x$ is below tolerance. The $(1 + \lvert x_n \rvert)$ factor makes this a relative test for large $x$ and an absolute test for small $x$.
2. **Residual is small:** $\lvert f(x_n) \rvert < \text{tol}$ — the function value is close to zero (backward error is small).
3. **Budget exhausted:** $n = N$ — maximum number of iterations reached.
4. **Step ill-defined:** singular Jacobian, zero denominator, etc.
5. **Tolerance:** can be as small as $O(\varepsilon)$ for root-finding (unlike optimisation, where $\sqrt{\varepsilon}$ is the practical limit).

---

## Part III — Optimization (Chapters 2, 6)

**Goal:** find a local minimum of $f : \mathbb{R}^{d} \to \mathbb{R}$. The global minimum is generally unattainable for non-convex $f$ (the problem is NP-hard in general).

### Tolerance Reality Check

Near a local minimum, the second-order Taylor expansion gives:

$$f(\mathbf{x}) - f(\mathbf{x}^{\ast}) \approx \tfrac{1}{2}(\mathbf{x} - \mathbf{x}^{\ast})^{T} H (\mathbf{x} - \mathbf{x}^{\ast})$$

Because $f$ is quadratic near the minimum, a perturbation of size $\delta$ in $\mathbf{x}$ causes a perturbation of size $\sim \delta^2$ in $f$. When $\delta^2 \sim \varepsilon$ (machine epsilon), then $f(\mathbf{x})$ and $f(\mathbf{x}^{\ast})$ are indistinguishable in floating point. This means:

$$\frac{\lVert \mathbf{x} - \mathbf{x}^{\ast} \rVert}{\lVert \mathbf{x}^{\ast} \rVert} = O(\sqrt{\varepsilon})$$

**Sensible default:** $\text{tol} = \sqrt{\varepsilon} \approx 10^{-8}$ in double precision. You **cannot** locate a minimum to full machine precision — this is a fundamental difference from root-finding, which can achieve $\text{tol} = O(\varepsilon)$.

### Termination Criteria

Stop iterating when any of the following hold:

1. **Step is small:** $\lVert \mathbf{x}_n - \mathbf{x}\_{n-1} \rVert < \text{tol}(1 + \lVert \mathbf{x}_n \rVert)$.
2. **Value change is small:** $\lvert f(\mathbf{x}_n) - f(\mathbf{x}\_{n-1}) \rvert < \text{tol}(1 + \lvert f(\mathbf{x}_n) \rvert)$ — but be cautious: $f$ could be creeping towards $-\infty$ with tiny steps.
3. **Gradient is small:** $\lVert \nabla f(\mathbf{x}_n) \rVert < \text{tol}(1 + \lVert \nabla f(\mathbf{x}_0) \rVert)$ — this is the best default when the gradient is available, since it directly measures proximity to a stationary point.
4. **Budget exhausted** or step ill-defined (singular Hessian, etc.).

### 1D Optimization

#### Golden Section Search

Requires a **bracket** $(a, b, c)$ with $a < b < c$, $f(b) < f(a)$, and $f(b) < f(c)$ — this guarantees a local minimum lies in $(a, c)$.

Let $\phi = (\sqrt{5} - 1)/2 \approx 0.618$. Each step selects a new test point $z$ in the **larger** of the two sub-intervals:

$$z = a + \phi(b-a) \qquad \text{or} \qquad z = c - \phi(c-b)$$

Compare $f(z)$ with $f(b)$ to form a new, smaller bracket from three of the four points $\{a, b, z, c\}$.

**Derivation of $\phi$:** we want the ratio of bracket reduction to be the same regardless of which sub-interval the new bracket falls in, *and* we want points from one iteration to be reusable in the next. This self-similarity requirement gives:

$$\phi^{2} + \phi = 1 \qquad\Longrightarrow\qquad \phi = \frac{\sqrt{5} - 1}{2}$$

- **Convergence:** linear with rate $\phi \approx 0.618$ per step (slightly slower than bisection's $0.5$ for root-finding — but optimal for derivative-free minimisation with one function evaluation per step).
- **Derivatives not required;** $f$ need only be continuous and unimodal on the bracket.
- **Use when:** $f$ is non-differentiable or derivatives are expensive.

#### Successive Parabolic Interpolation

Fit a parabola through the three bracket points $(a, f(a)), (b, f(b)), (c, f(c))$; take $z$ = vertex (minimum) of the parabola.

- **Convergence:** superlinear ($\approx$ order 1.325).
- **Pitfall:** can fail if the three points are collinear (parabola degenerates to a line), or the fitted parabola opens downward (vertex is a maximum).

#### Brent's Method (for Minimisation)

Combines parabolic interpolation (fast when applicable) with golden section search (safe fallback). Accepts the parabolic step only if it falls inside the bracket and is sufficiently small; otherwise defaults to golden section. Robust default 1D minimiser.

### Multi-Dimensional Optimization — Line Search Framework

All of the following methods share a common structure. At each step:

$$\mathbf{x}\_{n+1} = \mathbf{x}\_n + \alpha_n \mathbf{d}\_n$$

choosing:

1. **Direction** $\mathbf{d}_n$ — must be a **descent direction**, meaning the directional derivative is negative:

$$\mathbf{g}_n^{T} \mathbf{d}_n < 0$$

where $\mathbf{g}_n = \nabla f(\mathbf{x}_n)$. This guarantees that $f$ decreases along $\mathbf{d}_n$ for sufficiently small steps.

2. **Step length** $\alpha_n > 0$, chosen to approximately minimise $f(\mathbf{x}_n + \alpha \mathbf{d}_n)$ as a function of $\alpha$ (an "inexact line search").

3. **Infinite travel condition** to prevent stalling on a flat region:

$$\sum_{n=1}^{\infty} \alpha_n = \infty$$

This ensures the iterates can travel arbitrarily far from the starting point if needed.

4. Successive directions should not oscillate wildly (each undoing the progress of the previous step).

#### Step Length — Backtracking

A practical way to choose $\alpha_n$. Start with an initial guess:

$$\alpha' = \alpha_{n-1} \frac{\mathbf{g}\_{n-1}^{T} \mathbf{d}\_{n-1}}{\mathbf{g}\_n^{T} \mathbf{d}\_n}$$

(for $n = 0$, the user supplies $\alpha_0'$). This scales the previous step length by the ratio of directional derivatives. Then repeatedly shrink $\alpha \leftarrow \rho\alpha$ (typically $\rho = 1/2$) until the **Armijo (sufficient decrease) condition** holds:

$$f(\mathbf{x}\_n + \alpha \mathbf{d}\_n) < f(\mathbf{x}\_n) + \sigma \alpha \mathbf{g}\_n^{T} \mathbf{d}\_n, \qquad \sigma \in [10^{-4}, 10^{-1}]$$

The right-hand side is the function value predicted by a linear model, scaled down by $\sigma$. The condition says the actual decrease must be at least a fraction $\sigma$ of the predicted decrease.

Backtracking is guaranteed to terminate whenever $\mathbf{d}_n$ is a descent direction (because for sufficiently small $\alpha$, the Taylor expansion ensures the Armijo condition holds).

**Wolfe conditions** (not implemented in detail here): Armijo plus a **curvature condition** that prevents the step from being too short. Together they ensure both sufficient decrease and sufficient progress.

### Direction Choices

#### 1. Coordinate Descent

$\mathbf{d}_n$ cycles through the coordinate directions $\pm\mathbf{e}_i$. No derivatives are needed; each step minimises along one axis.

- Each step is very cheap; but convergence is slow, especially when variables are correlated.
- The method can get stuck on non-axis-aligned valleys.

#### 2. Gradient (Steepest) Descent

$$\mathbf{d}_n = -\mathbf{g}_n$$

This is automatically a descent direction (since $\mathbf{g}_n^T(-\mathbf{g}_n) = -\lVert\mathbf{g}_n\rVert^2 < 0$). It is the direction that gives the fastest instantaneous rate of decrease of $f$: the argmin of $\mathbf{g}^{T}\mathbf{d}$ over $\lVert \mathbf{d} \rVert \leq 1$.

- Tends to **zigzag** in narrow valleys (when the Hessian is ill-conditioned — the eigenvalues span a wide range). Each step heads straight towards the valley floor but overshoots sideways.
- Converges linearly for well-behaved convex problems; the rate depends on the condition number of the Hessian.
- Simple and robust, but slow for ill-conditioned problems.

#### 3. Newton's Method

$$\mathbf{d}_n = -H(f)(\mathbf{x}_n)^{-1} \mathbf{g}_n, \qquad \alpha_n = 1 \text{ by default}$$

**Derivation:** Approximate $f$ by its second-order Taylor expansion at $\mathbf{x}_n$:

$$f(\mathbf{x}_n + \mathbf{d}) \approx f(\mathbf{x}_n) + \mathbf{g}_n^T\mathbf{d} + \tfrac{1}{2}\mathbf{d}^T H \mathbf{d}$$

Minimising this quadratic in $\mathbf{d}$ by setting the gradient to zero gives $H\mathbf{d} = -\mathbf{g}_n$. In practice, **solve** $H \mathbf{d}_n = -\mathbf{g}_n$ rather than forming $H^{-1}$.

- **Convergence:** quadratic near a local minimum with positive-definite Hessian.
- **Cost:** $O(d^{2})$ second-derivative evaluations plus $O(d^{3})$ linear solve per step.
- **Failure modes (important!):**
  - $H$ not positive definite $\to$ $\mathbf{d}_n$ may not be a descent direction; the method can race towards a maximum or saddle point.
  - $H$ singular $\to$ step undefined.
  - Non-convex regions $\to$ divergence is likely without safeguards.

**Relation to Chapter 5:** Newton's method for optimisation is exactly Newton's method for root-finding applied to $\nabla f = \mathbf{0}$. However, because a zero gradient doesn't distinguish minima from maxima or saddle points, the optimisation version needs extra safeguards (descent check, Hessian positive-definiteness).

**Why not just minimise $\lVert \mathbf{f}(\mathbf{x})\rVert^{2}$ to solve a root-finding problem?** Because $\lVert \mathbf{f}(\mathbf{x})\rVert^{2}$ is rarely convex; it typically has spurious local minima where $\mathbf{f} \neq \mathbf{0}$. Optimisation methods can get trapped at these, failing to find real roots. Useful only to generate a starting point for the proper root-finding methods of Chapter 5. **Exception:** when $\mathbf{f}$ is affine (i.e. $\mathbf{f}(\mathbf{x}) = A\mathbf{x} - \mathbf{b}$), then $\lVert \mathbf{f} \rVert^{2}$ is convex, and gradient/conjugate-gradient methods solve the system cleanly.

#### 4. BFGS (Quasi-Newton)

Build up a **symmetric positive-definite** approximation $\hat{H}_n$ to the Hessian via rank-2 updates. Let $\Delta\mathbf{x} = \mathbf{x}_n - \mathbf{x}\_{n-1}$ and $\Delta\mathbf{g} = \mathbf{g}_n - \mathbf{g}\_{n-1}$:

$$\hat{H}\_n = \hat{H}\_{n-1} + \frac{\Delta\mathbf{g} \Delta\mathbf{g}^{T}}{\Delta\mathbf{g}^{T} \Delta\mathbf{x}} - \frac{\hat{H}\_{n-1} \Delta\mathbf{x} \Delta\mathbf{x}^{T} \hat{H}\_{n-1}}{\Delta\mathbf{x}^{T} \hat{H}\_{n-1} \Delta\mathbf{x}}$$

This update adds information from the latest step while preserving symmetry.

**Curvature condition:** The update preserves positive-definiteness only if:

$$\Delta\mathbf{x}^{T} \Delta\mathbf{g} > 0$$

For convex $f$ this always holds (since the gradient is monotone). For non-convex $f$ it can fail for small step sizes. Simple recovery: **reset $\hat{H}_n = I$** — the next step becomes gradient descent, and BFGS can rebuild a good approximation from there.

Direction: $\mathbf{d}_n = -\hat{H}_n^{-1} \mathbf{g}_n$. Use the Sherman–Morrison–Woodbury formula to update $\hat{H}^{-1}$ directly in $O(d^2)$, avoiding the $O(d^{3})$ linear system solve at each step.

- **Initialisation:** $\hat{H}_0 = I$ makes the first step equivalent to gradient descent. As iterations proceed, $\hat{H}_n$ accumulates curvature information and approaches the true Hessian on well-behaved problems.
- **Convergence:** superlinear (not quadratic in general, because the Hessian approximation lags behind); globally convergent under mild conditions with a line search.
- **No second derivatives required** — huge savings vs Newton, especially in high dimensions.
- **L-BFGS** is the limited-memory variant that stores only the last $m$ update vectors instead of the full $d \times d$ matrix. Standard in large-scale machine learning.

#### 5. Stochastic / Minibatch Gradient Descent

When $f$ decomposes as a sum (typical in machine learning):

$$l(\mathbf{w}) = \frac{1}{N}\sum_{i=1}^{N} l_i(\mathbf{w})$$

computing the full gradient $\nabla l$ requires processing all $N$ data points — expensive when $N$ is large.

**Minibatch gradient:** compute the gradient over a random subset $B \subset \{1, \ldots, N\}$:

$$\mathbf{d}_n = -\frac{1}{\lvert B \rvert}\sum_{i \in B} \nabla l_i(\mathbf{w}_n)$$

This is an unbiased estimate of $-\nabla l(\mathbf{w}_n)$ with variance proportional to $1/\lvert B \rvert$.

**Stochastic gradient descent (SGD):** the extreme case $\lvert B \rvert = 1$ — one random data point per step. Much noisier but vastly cheaper per iteration. The noise can sometimes help escape shallow local minima (an advantage in non-convex optimisation).

### Comparison (Multi-Dimensional)

| Method | Needs | Per-step cost | Convergence | Notes |
|---|---|---|---|---|
| Coordinate descent | $f$ only | $O(d)$ | Slow | Naive; no derivative info |
| Gradient descent | $\nabla f$ | $O(d)$ | Linear | Zigzags; tuning-sensitive |
| Newton | $\nabla f, H$ | $O(d^{3})$ | Quadratic near min | Fragile in non-convex regions |
| BFGS | $\nabla f$ | $O(d^{2})$ | Superlinear | Robust, industry default |

### When to Use What — Optimization

- $d = 1$, non-smooth $\to$ golden section or Brent.
- $d = 1$, smooth $\to$ 1D Newton on $f'$ (i.e. Newton-opt).
- Convex, small $d$, cheap Hessian $\to$ Newton with backtracking.
- Non-convex or larger $d$ $\to$ BFGS (or L-BFGS for very large $d$).
- Huge $d$, ML-style sum structure $\to$ SGD or minibatch.
- Large sparse linear system $A\mathbf{x} = \mathbf{b}$ (an implicitly convex quadratic $\frac{1}{2}\mathbf{x}^T A \mathbf{x} - \mathbf{b}^T\mathbf{x}$) $\to$ gradient or conjugate-gradient methods.

### Optimization vs Root-Finding — Summary

| | Optimization | Root-Finding |
|---|---|---|
| Best attainable tolerance | $\sqrt{\varepsilon} \approx 10^{-8}$ | $\varepsilon \approx 10^{-16}$ |
| "Better" iterate means | $f(\mathbf{x}\_{n+1}) < f(\mathbf{x}\_n)$ | Smaller residual $\lVert\mathbf{f}(\mathbf{x}\_n)\rVert$ (not always monotone) |
| Newton variant | Quadratic near min; fails if $H$ not PD | Quadratic near simple root; fails if $J$ singular |
| Second-order structure | Hessian (symmetric) | Jacobian (generally not symmetric) |
| Key difference | Must ensure **descent** (not just stationarity) | Must ensure **residual** shrinks |

---

## Required Proofs (ON Syllabus — be able to reproduce)

1. **Composite integration error bounds** (§3.4): derive the composite midpoint, trapezium, and Simpson's bounds

$$\frac{(b-a)^{3}}{24 n^{2}} \max \lvert f'' \rvert, \qquad \frac{(b-a)^{3}}{12 n^{2}} \max \lvert f'' \rvert, \qquad \frac{(b-a)^{5}}{180 n^{4}} \max \lvert f^{(4)} \rvert$$

from the single-strip error formulas, summing over strips and using continuity of derivatives / the intermediate-value theorem.

2. **Monte Carlo convergence**: unbiasedness ($\mathbb{E}[\mathrm{MC}\_N] = I$); variance formula $\mathrm{Var}(\mathrm{MC}\_N) = A(R)^{2} \mathrm{Var}(f)/N$; error rate $O(N^{-1/2})$ via Chebyshev's inequality.

3. **Lemma 4.2**: if $A_n$ has order-$q$ convergence then $A_{2n}$ has order $q^{2}$.

4. **Lemma 5.3 (1D Newton quadratic convergence)**: the identity

$$\epsilon_{n+1} = -\frac{f''(\xi_n)}{2 f'(x_n)} \epsilon_n^{2}$$

and the consequent quadratic rate.

5. **Stationary-point classification** via Hessian definiteness (from the second-order Taylor expansion).

6. **Convexity results**: local min equals global min for convex $f$; twice-differentiable $f$ is convex iff Hessian is PSD; Jensen's inequality.

7. **Golden-ratio derivation** for golden section search: $\phi^{2} + \phi = 1$.

8. **Lagrange multiplier method** — statement and usage (the proof is bonus; see below).

9. **Convergence of Newton's method for 1D optimisation** (quadratic, derived by applying 1D Newton root-finding to $f'$).

---

## Bonus / NOT Part of the Syllabus (for reference only)

The notes explicitly flag the following as outside the syllabus. Do not study the content; just know it exists so you can skip it.

1. **§1.8** — Proof of Taylor's theorem (nD via 1D construction). Use the statement; skip the proof.
2. **§2.8** — Proof of the Lagrange multiplier method. Use the method; the proof (IFT / geometric argument) is optional.
3. **KKT conditions** (inequality-constrained optimisation) — "we will not use this method this year".
4. **§3.9** — Single-strip trapezium and Simpson's error-bound *derivations* via Rolle-type arguments. (Composite bounds ARE required; the single-strip proofs are bonus.)
5. **§4.5** — Acceleration methods: Aitken's $\Delta^{2}$ process, Steffensen's method, Theorem 4.4.
6. **§5.8** — Convergence of Newton's method in $d$ dimensions (nD quadratic-convergence proof). The 1D version IS required.
7. **§6.6** — Linear-convergence proof for gradient descent on strongly-convex Lipschitz functions (Lemma 6.3, Theorem 6.4, condition-number argument).
8. Derivation of conjugate gradient descent (Polak–Ribière).
9. Derivation of the Sherman–Morrison update for BFGS's inverse Hessian.
10. Hessian-modification schemes for Newton's method (adding $\lambda I$, Cholesky-based $\lambda$ search).
11. The distinction between forward and backward error is mentioned but "outside our syllabus".
12. Details of the Wolfe conditions beyond their statement.
13. The Nelder–Mead / simplex method.
14. Internals of Brent's hybrid-switching logic.

---

## Quick Reference — "Which Method?"

**Integration**

- Smooth $f$, 1D, default $\to$ composite Simpson's.
- Rough $f$ or quick estimate $\to$ composite trapezium.
- Endpoint singularity $\to$ composite midpoint.
- $d \geq 5$ or irregular region $\to$ Monte Carlo (+ variance reduction if needed).

**Root-Finding (1D)**

- Have a bracket, want guarantee $\to$ bisection or Brent.
- Smooth, good starting guess, cheap derivative $\to$ Newton.
- Smooth, good guess, no derivative $\to$ secant.
- Production default $\to$ Brent.

**Root-Finding (nD)**

- Jacobian cheap $\to$ Newton.
- Jacobian expensive $\to$ Broyden.

**Optimization (1D)**

- Non-smooth or derivative-free $\to$ golden section or Brent.
- Smooth $\to$ 1D Newton on $f'$.

**Optimization (nD)**

- Small $d$, convex, cheap Hessian $\to$ Newton + backtracking.
- Default robust workhorse $\to$ BFGS.
- Huge $d$ with sum structure (ML) $\to$ SGD / minibatch.
- Large linear system $A\mathbf{x} = \mathbf{b}$ $\to$ gradient or conjugate-gradient methods.
