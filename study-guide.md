# Continuous Mathematics — Numerical Methods Study Guide

Covers **Integration**, **Root-Finding**, **Optimization**. For each method: formula, applicable dimension, when to use, conditions, error/convergence. Required proofs vs bonus/non-syllabus material are listed at the end.

---

## Foundational Background (Chapters 1–2, 4)

### Derivatives (nD)

- **Gradient** (column vector):

$$\nabla f = \left(\frac{\partial f}{\partial x_1}, \ldots, \frac{\partial f}{\partial x_d}\right)^{T}$$

- **Hessian** — symmetric when second partials continuous:

$$H(f)_{ij} = \frac{\partial^{2} f}{\partial x_i\, \partial x_j}$$

- **Jacobian** of $f : \mathbb{R}^n \to \mathbb{R}^m$:

$$J_{ij} = \frac{\partial f_i}{\partial x_j}$$

Note: Hessian of scalar $f$ equals the Jacobian of $\nabla f$.

### Taylor's Theorem

**1D, order $n$ (Lagrange remainder):**

$$f(x+h) \;=\; \sum_{k=0}^{n} \frac{h^{k}}{k!}\, f^{(k)}(x) \;+\; \frac{h^{n+1}}{(n+1)!}\, f^{(n+1)}(\xi)$$

for some $\xi$ between $x$ and $x+h$.

**nD, second order:**

$$f(\mathbf{x}+\mathbf{h}) \;=\; f(\mathbf{x}) + \mathbf{h}^{T}\nabla f(\mathbf{x}) + \tfrac{1}{2}\,\mathbf{h}^{T}\,H(f)(\mathbf{x}+\theta\mathbf{h})\,\mathbf{h}$$

for some $\theta \in (0,1)$.

### Convexity

- $f$ convex on convex set $S$ iff for all $\mathbf{x}, \mathbf{y} \in S$, $\lambda \in [0,1]$:

$$f(\lambda\mathbf{x} + (1-\lambda)\mathbf{y}) \;\leq\; \lambda f(\mathbf{x}) + (1-\lambda)f(\mathbf{y})$$

- Twice-differentiable $f$ convex $\iff$ Hessian positive semidefinite everywhere. Strictly convex $\impliedby$ Hessian positive definite.
- Convex $f$: every local min is global; stationary point $\Rightarrow$ global min.

### Jensen's Inequality

For convex $f$ and probability weights $\lambda_i \geq 0$ with $\sum \lambda_i = 1$:

$$f\!\left(\sum_i \lambda_i \mathbf{x}_i\right) \;\leq\; \sum_i \lambda_i f(\mathbf{x}_i)$$

### Stationary Point Classification (nD)

At stationary point $\mathbf{x}^{\ast}$ (gradient zero):

- $H(f)(\mathbf{x}^{\ast})$ positive definite $\to$ strict local **minimum**.
- Negative definite $\to$ strict local **maximum**.
- Indefinite $\to$ **saddle point**.
- Semidefinite (singular) $\to$ inconclusive; need higher-order test.

**Tests for definiteness:** all eigenvalues positive / all leading principal minors (pivots) positive / Cholesky factorisation succeeds.

### Accuracy (Chapter 4)

- **Error**: $\tilde{u} - u$. **Absolute error**: $\lvert \tilde{u} - u \rvert$. **Relative error**: $\lvert \tilde{u} - u \rvert / \lvert u \rvert$.
- For vectors, replace $\lvert\cdot\rvert$ by Euclidean norm $\lVert\cdot\rVert$.
- **Truncation error**: error from approximating an infinite procedure by a finite one (e.g. Taylor cutoff, $n$-strip integration).
- **Roundoff error**: error from finite floating-point representation. *Not* the same as truncation error.
- **Machine epsilon**: $\varepsilon \approx 1.1 \times 10^{-16}$ (double), $\approx 6 \times 10^{-8}$ (single). Smallest relative rounding error on storing a real.
- **Catastrophic cancellation**: subtracting two nearly equal floating-point numbers produces huge relative error. Often avoidable by algebraic reformulation (e.g. $1-\cos x \to 2\sin^{2}(x/2)$).
- **Forward error**: $\tilde{f}(x) - f(x)$. **Backward error**: $\tilde{x} - x$ where $\tilde{f}(x) = f(\tilde{x})$.

### Rates of Convergence (Ch 4)

Let $\epsilon_n = x_n - x^{\ast}$. As $n \to \infty$:

| Type | Condition |
|---|---|
| Linear (order 1) | $\lvert \epsilon_{n+1}\rvert / \lvert \epsilon_n \rvert \to a$, with $0 < a < 1$ |
| Sublinear | $\lvert \epsilon_{n+1}\rvert / \lvert \epsilon_n \rvert \to 1$ |
| Logarithmic | sublinear plus $\lvert \epsilon_{n+2}-\epsilon_{n+1}\rvert / \lvert \epsilon_{n+1}-\epsilon_n\rvert \to 1$ |
| Superlinear | $\lvert \epsilon_{n+1}\rvert / \lvert \epsilon_n \rvert \to 0$ |
| Order $q>1$ | $\lvert \epsilon_{n+1}\rvert / \lvert \epsilon_n\rvert^{q} \to a$, with $a > 0$ |
| Quadratic | order $2$ |

**Lemma 4.2 (REQUIRED):** If $A_n$ converges with order $q$, then $A_{2n}$ converges with order $q^{2}$. Lets you compare a fast-convergent slow method vs a slow-convergent fast method.

Rule of thumb: linear $\Rightarrow$ adds a constant number of correct decimals per step; quadratic $\Rightarrow$ doubles correct decimals per step (once error is small).

### Lagrange Multipliers (equality constraints)

To minimize $f(\mathbf{x})$ subject to $g_i(\mathbf{x}) = 0$ for $i = 1, \ldots, m$: introduce the Lagrangian

$$L(\mathbf{x}, \boldsymbol{\lambda}) \;=\; f(\mathbf{x}) - \sum_{i=1}^{m} \lambda_i\, g_i(\mathbf{x})$$

and solve $\nabla_{\mathbf{x}} L = \mathbf{0}$ together with $\nabla_{\boldsymbol{\lambda}} L = \mathbf{0}$ (the latter reproduces the constraints).

**KKT conditions (for inequality constraints) are NOT on the syllabus this year.**

---

## Part I — Numerical Integration (Chapter 3)

**Goal:** approximate $\displaystyle \int_a^b f(x)\,dx$.

Most classical rules are **Newton–Cotes**: replace $f$ by a polynomial interpolant on equally-spaced nodes, integrate the polynomial exactly.

### 1D Rules (single strip of width $h = b-a$)

**Midpoint rule.** One node $m = (a+b)/2$:

$$\int_a^b f(x)\,dx \;\approx\; (b-a)\, f(m)$$

Error: $\displaystyle -\frac{(b-a)^{3}}{24} f''(\xi)$. Exact for polynomials of degree $\leq 1$.

**Trapezium rule.** Two nodes $a, b$:

$$\int_a^b f(x)\,dx \;\approx\; \frac{b-a}{2}\bigl(f(a) + f(b)\bigr)$$

Error: $\displaystyle \frac{(b-a)^{3}}{12} f''(\xi)$. Exact for polynomials of degree $\leq 1$.

**Simpson's rule.** Three nodes $a, m, b$:

$$\int_a^b f(x)\,dx \;\approx\; \frac{b-a}{6}\bigl(f(a) + 4f(m) + f(b)\bigr)$$

Error: $\displaystyle -\frac{(b-a)^{5}}{2880} f^{(4)}(\xi)$. Exact for polynomials of degree $\leq 3$ (one degree better than expected thanks to symmetry).

**Boole's rule.** Five equally-spaced nodes:

$$\int_a^b f(x)\,dx \;\approx\; \frac{b-a}{90}\bigl(7 f_0 + 32 f_1 + 12 f_2 + 32 f_3 + 7 f_4\bigr)$$

Error: $O\bigl((b-a)^{7} f^{(6)}\bigr)$. Exact for polynomials of degree $\leq 5$.

### Composite Rules ($n$ strips of width $h = (b-a)/n$)

Sum single-strip rules over $n$ sub-intervals.

**Composite midpoint:**

$$\int_a^b f(x)\,dx \;\approx\; h\sum_{i=1}^{n} f(m_i), \qquad \lvert E \rvert \;\leq\; \frac{(b-a)^{3}}{24\, n^{2}} \max_{x \in [a,b]} \lvert f''(x) \rvert$$

**Composite trapezium:**

$$\int_a^b f(x)\,dx \;\approx\; h\left(\tfrac{1}{2}f(a) + \sum_{i=1}^{n-1} f(x_i) + \tfrac{1}{2}f(b)\right), \qquad \lvert E \rvert \;\leq\; \frac{(b-a)^{3}}{12\, n^{2}} \max \lvert f'' \rvert$$

**Composite Simpson's** ($n$ even):

$$\int_a^b f(x)\,dx \;\approx\; \frac{h}{3}\!\left( f_0 + 4(f_1 + f_3 + \cdots + f_{n-1}) + 2(f_2 + f_4 + \cdots + f_{n-2}) + f_n \right)$$

$$\lvert E \rvert \;\leq\; \frac{(b-a)^{5}}{180\, n^{4}} \max \lvert f^{(4)} \rvert$$

**Rates:** midpoint and trapezium are $O(n^{-2})$; Simpson's is $O(n^{-4})$; Boole's is $O(n^{-6})$.

### When to Use Which

- $f$ only $C^{2}$ smooth $\to$ midpoint or trapezium.
- $f$ is $C^{4}$ smooth $\to$ Simpson's (the default workhorse).
- Very smooth $f$ $\to$ Boole's or higher (caveat: high-order Newton–Cotes has negative weights and becomes unstable; Gaussian quadrature is preferred but is not in this syllabus).
- Endpoint singularities $\to$ midpoint (avoids evaluating at endpoints).
- Non-smooth regions $\to$ split the interval, or rate will degrade.

### Comparison

| Feature | Midpoint | Trapezium | Simpson's |
|---|---|---|---|
| $f$ evaluations (composite, shared) | $n$ | $n+1$ | $2n+1$ |
| Convergence rate | $O(n^{-2})$ | $O(n^{-2})$ | $O(n^{-4})$ |
| Evaluates $f$ at endpoints? | No | Yes | Yes |
| Exact for polynomial degree | $\leq 1$ | $\leq 1$ | $\leq 3$ |

Midpoint error is roughly half the magnitude of trapezium and opposite in sign; Simpson's is the linear combination $(2\cdot\text{midpoint} + \text{trapezium})/3$.

### Multi-Dimensional Integration

- **Tensor-product Newton–Cotes** (e.g. Simpson's in $d$ dimensions) costs $O(n^{d})$ evaluations: the **curse of dimensionality**. Rate is still $O(n^{-4})$ in the per-axis resolution, but the evaluation cost scales exponentially in $d$.
- For large $d$ or irregular regions, switch to **Monte Carlo**.

### Monte Carlo Integration

Estimator: draw $X_1, \ldots, X_N$ i.i.d. uniform on region $R \subset \mathbb{R}^{d}$ with volume $A(R)$:

$$\mathrm{MC}_{N} \;=\; A(R) \cdot \frac{1}{N} \sum_{i=1}^{N} f(X_i)$$

**Unbiased:**

$$\mathbb{E}[\mathrm{MC}_{N}] \;=\; \int_{R} f\,d\mathbf{x}$$

**Variance:**

$$\mathrm{Var}(\mathrm{MC}_{N}) \;=\; \frac{A(R)^{2}\, \mathrm{Var}(f(X))}{N}$$

**RMS error:** $O(N^{-1/2})$ — **independent of dimension $d$**. This is the reason MC wins in high dimension.

**Chebyshev's inequality** gives a cheap confidence bound:

$$\mathbb{P}\!\left(\lvert \mathrm{MC}_{N} - I \rvert \geq \frac{k\,\sigma}{\sqrt{N}}\right) \;\leq\; \frac{1}{k^{2}}$$

**Central Limit Theorem:** for large $N$, $\mathrm{MC}_{N}$ is approximately normal, giving tighter Gaussian confidence intervals.

#### Variance Reduction

- **Stratified sampling**: partition $R$ into sub-regions, sample independently inside each. Variance is never worse than plain MC.
- **Importance sampling**: draw from density $p(\mathbf{x})$ proportional to $\lvert f(\mathbf{x})\rvert$, weight by $f(X)/p(X)$. Huge gains when $p$ resembles $\lvert f \rvert$.
- **Control variates**: find $\tilde{f}$ with known integral; estimate $\int (f - \tilde{f})$ by MC and add the known $\int \tilde{f}$. Works when $f - \tilde{f}$ has smaller variance than $f$.

### When to Use Monte Carlo

- $d \geq \sim 5$ (tensor grids infeasible).
- Irregular region $R$.
- $f$ only available as a simulation output.
- Low-to-moderate accuracy acceptable ($\sqrt{N}$ is slow).

---

## Part II — Root-Finding (Chapter 5)

**Goal:** find $\mathbf{x}^{\ast}$ with $\mathbf{f}(\mathbf{x}^{\ast}) = \mathbf{0}$.

### 1D — Bracketing Method

#### Interval Bisection

Requires a **bracket**: $a < b$ with $f(a) f(b) < 0$ and $f$ continuous on $[a,b]$ — guarantees at least one root inside.

Iteration: set $m = (a+b)/2$; replace the endpoint that preserves the sign change.

- **Convergence:** linear, $\lvert \epsilon_{n+1}\rvert / \lvert \epsilon_n \rvert \leq 1/2$ always.
- **Error bound:** $\lvert \epsilon_n \rvert \leq (b_0 - a_0)/2^{\,n+1}$.
- **When to use:** guaranteed convergence; no derivatives needed; very robust.

### 1D — Open Methods (no bracket)

#### Newton's Method (1D)

$$x_{n+1} \;=\; x_n - \frac{f(x_n)}{f'(x_n)}$$

Derivation: first-order Taylor expansion of $f$ at $x_n$, solve linear approximation for zero.

- **Convergence:** quadratic if $f'(x^{\ast}) \neq 0$, $f''$ continuous near $x^{\ast}$, and $x_0$ sufficiently close.
- At a multiple root ($f'(x^{\ast}) = 0$) convergence degrades to linear.
- **Pitfalls:** can diverge or overshoot; requires $f'$.

#### Secant Method

Replace $f'(x_n)$ with a finite-difference estimate from the previous two iterates:

$$x_{n+1} \;=\; x_n \;-\; f(x_n) \cdot \frac{x_n - x_{n-1}}{f(x_n) - f(x_{n-1})}$$

- **Convergence:** superlinear, order $\varphi = (1+\sqrt{5})/2 \approx 1.618$ (golden ratio).
- One new $f$ evaluation per step (cache the previous value).
- **Pitfall:** denominator becomes tiny near a root $\to$ catastrophic cancellation; no bracket.

#### Halley's Method

$$x_{n+1} \;=\; x_n \;-\; \frac{2\, f(x_n)\, f'(x_n)}{2\, f'(x_n)^{2} \;-\; f(x_n)\, f''(x_n)}$$

- **Convergence:** cubic (order 3).
- **Pitfall:** needs both $f'$ and $f''$; expensive per step.

#### Muller's Method

Fit a **quadratic** through $(x_n, x_{n-1}, x_{n-2})$; take $x_{n+1}$ = root of that quadratic closest to $x_n$.

- **Convergence:** superlinear, order $\approx 1.839$.
- **Use:** can find complex roots even from real iterates (quadratic formula can produce complex values).

#### Inverse Quadratic Interpolation (IQI)

Fit a quadratic to the points $(f(x_i), x_i)$ — swap roles of $x$ and $y$ — and take $x_{n+1}$ as the value where $y = 0$.

- **Order:** $\approx 1.839$ (same as Muller).
- **Pitfall:** fails if two $f$ values coincide.

#### Brent's Method

Hybrid of bisection, secant, and IQI. Maintains a bracket at all times. Accepts a superlinear (secant/IQI) step only if it lies inside a "safe" region; otherwise falls back to bisection.

- **Convergence:** superlinear in typical cases; linear worst-case guarantee because of the bracket.
- **Use:** robust default 1D root-finder (the classic `zeroin` implementation guarantees relative error $\leq 4\varepsilon$).

### Comparison (1D)

| Method | Order | Needs derivatives | Bracket | Guaranteed? |
|---|---|---|---|---|
| Bisection | 1 (linear) | No | Yes | Yes |
| Newton | 2 | $f'$ | No | No |
| Secant | $\approx 1.618$ | No | No | No |
| Halley | 3 | $f', f''$ | No | No |
| Muller / IQI | $\approx 1.839$ | No | No | No |
| Brent | $\geq 1$ (typ. $\approx 1.839$) | No | Yes | Yes |

### Error Analysis — Newton's Method (Lemma 5.3, REQUIRED)

For 1D Newton near a simple root with $f''$ continuous and $f'(x^{\ast}) \neq 0$:

$$\epsilon_{n+1} \;=\; -\frac{f''(\xi_n)}{2\, f'(x_n)}\, \epsilon_n^{2}$$

for some $\xi_n$ between $x_n$ and $x^{\ast}$. Hence

$$\frac{\lvert \epsilon_{n+1} \rvert}{\lvert \epsilon_n \rvert^{2}} \;\to\; \frac{\lvert f''(x^{\ast}) \rvert}{2\,\lvert f'(x^{\ast}) \rvert}$$

establishing **quadratic convergence**.

### A Posteriori Error Estimate

Given a candidate $\tilde{x}$, estimate the distance to a root via a Taylor expansion:

$$\lvert \tilde{x} - x^{\ast} \rvert \;\approx\; \frac{\lvert f(\tilde{x}) \rvert}{\lvert f'(\tilde{x}) \rvert}$$

(exact when $f$ is linear).

### Multi-Dimensional Root-Finding

$\mathbf{f} : \mathbb{R}^{d} \to \mathbb{R}^{d}$, find $\mathbf{x}^{\ast}$ with $\mathbf{f}(\mathbf{x}^{\ast}) = \mathbf{0}$.

#### Newton's Method (nD)

Linearise: $\mathbf{f}(\mathbf{x}) \approx \mathbf{f}(\mathbf{x}_n) + J(\mathbf{x}_n)(\mathbf{x} - \mathbf{x}_n)$.

Iteration: **solve** the linear system

$$J(\mathbf{x}_n)\, \Delta\mathbf{x} \;=\; -\mathbf{f}(\mathbf{x}_n), \qquad \mathbf{x}_{n+1} \;=\; \mathbf{x}_n + \Delta\mathbf{x}$$

Do **not** form $J^{-1}$; solve directly (LU, Cholesky, etc.).

- **Cost:** $O(d^{3})$ per step plus $d^{2}$ partial-derivative evaluations.
- **Convergence:** quadratic if $J(\mathbf{x}^{\ast})$ is nonsingular and $\mathbf{x}_0$ is close.
- **Fails** if $J(\mathbf{x}_n)$ is singular.

#### Broyden's Method (Quasi-Newton)

Avoid the $d^{2}$ derivatives by building up an approximate Jacobian $\hat{J}_n$ via rank-1 updates. Let $\Delta\mathbf{x} = \mathbf{x}_n - \mathbf{x}_{n-1}$ and $\Delta\mathbf{f} = \mathbf{f}(\mathbf{x}_n) - \mathbf{f}(\mathbf{x}_{n-1})$:

$$\hat{J}_n \;=\; \hat{J}_{n-1} \;+\; \frac{\Delta\mathbf{f} - \hat{J}_{n-1}\, \Delta\mathbf{x}}{\lVert \Delta\mathbf{x} \rVert^{2}}\, \Delta\mathbf{x}^{T}$$

Satisfies the **secant equation** $\hat{J}_n\, \Delta\mathbf{x} = \Delta\mathbf{f}$.

- **Convergence:** superlinear (not quadratic).
- Use the **Sherman–Morrison formula** to update $\hat{J}^{-1}$ directly in $O(d^{2})$ — avoids the $O(d^{3})$ solve.
- Initial $\hat{J}_0$: identity, or the true $J(\mathbf{x}_0)$ if cheap.

### Termination Criteria (Root-Finding)

- $\lvert x_n - x_{n-1} \rvert < \text{tol} \cdot (1 + \lvert x_n \rvert)$ (step small).
- $\lvert f(x_n) \rvert < \text{tol}$ (residual / backward error).
- $n = N$ (iteration budget exhausted).
- Step ill-defined (singular Jacobian, zero denominator).
- $\text{tol}$ can be as small as $O(\varepsilon)$ for root-finding (unlike optimization).

---

## Part III — Optimization (Chapters 2, 6)

**Goal:** find a local minimum of $f : \mathbb{R}^{d} \to \mathbb{R}$. The global minimum is generally unattainable for non-convex $f$.

### Tolerance Reality Check

Near a local minimum, a second-order Taylor expansion gives

$$f(\mathbf{x}) - f(\mathbf{x}^{\ast}) \;\approx\; \tfrac{1}{2}(\mathbf{x} - \mathbf{x}^{\ast})^{T}\, H\, (\mathbf{x} - \mathbf{x}^{\ast})$$

So $f(\mathbf{x})$ and $f(\mathbf{x}^{\ast})$ become indistinguishable in floating point when $\lVert \mathbf{x} - \mathbf{x}^{\ast} \rVert / \lVert \mathbf{x}^{\ast} \rVert = O(\sqrt{\varepsilon})$.

**Sensible default:** $\text{tol} = \sqrt{\varepsilon} \approx 10^{-8}$ in double precision. You **cannot** locate a minimum to full machine precision.

### Termination Criteria

- $\lVert \mathbf{x}_n - \mathbf{x}_{n-1} \rVert < \text{tol}\,(1 + \lVert \mathbf{x}_n \rVert)$ — step small.
- $\lvert f(\mathbf{x}_n) - f(\mathbf{x}_{n-1}) \rvert < \text{tol}\,(1 + \lvert f(\mathbf{x}_n) \rvert)$ — value change small (risky: $f$ could be creeping to $-\infty$).
- $\lVert \nabla f(\mathbf{x}_n) \rVert < \text{tol}\,(1 + \lVert \nabla f(\mathbf{x}_0) \rVert)$ — gradient small (good default when gradient is available).
- $n = N$ (budget) or the iteration step was ill-defined.

### 1D Optimization

#### Golden Section Search

Requires a **bracket** $(a, b, c)$ with $a < b < c$, $f(b) < f(a)$, and $f(b) < f(c)$ — guarantees a local minimum lies in $(a, c)$.

Let $\phi = (\sqrt{5} - 1)/2 \approx 0.618$. Each step picks a new point $z$ in the larger sub-interval via

$$z = a + \phi(b-a) \quad \text{or} \quad z = c - \phi(c-b)$$

and keeps the three of $\{a, b, z, c\}$ that remain a bracket.

**Derivation of $\phi$:** requiring both (i) equal candidate-bracket lengths and (ii) the ratio to persist across iterations gives

$$\phi^{2} + \phi = 1 \quad\Longrightarrow\quad \phi = \frac{\sqrt{5} - 1}{2}$$

- **Convergence:** linear with rate $\phi \approx 0.618$ per step (slightly slower than bisection's $0.5$ for root-finding — but optimal for derivative-free minimization).
- **Derivatives not required;** $f$ need only be continuous and unimodal on the bracket.
- **Use when:** $f$ is non-differentiable or derivatives are expensive.

#### Successive Parabolic Interpolation

Fit a parabola through $(a, f(a)), (b, f(b)), (c, f(c))$; new $z$ = vertex of the parabola.

- **Convergence:** superlinear ($\approx$ order 1.325).
- **Pitfall:** can fail if the three points are collinear, or the vertex is a maximum rather than a minimum.

#### Brent's Method (for Minimization)

Combines parabolic interpolation (fast when applicable) with golden section search (safe fallback). Robust default 1D minimizer.

### Multi-Dimensional Optimization — Line Search Framework

Iterate

$$\mathbf{x}_{n+1} \;=\; \mathbf{x}_n + \alpha_n\, \mathbf{d}_n$$

choosing at each step:

1. **Direction** $\mathbf{d}_n$. Must be a **descent direction**: with $\mathbf{g}_n = \nabla f(\mathbf{x}_n)$, require

$$\mathbf{g}_n^{T}\, \mathbf{d}_n \;<\; 0$$

2. **Step length** $\alpha_n > 0$, approximately minimising $f(\mathbf{x}_n + \alpha\, \mathbf{d}_n)$.

3. **Infinite travel condition** to prevent stalling:

$$\sum_{n=1}^{\infty} \alpha_n \;=\; \infty$$

4. Successive directions should not veer wildly (each undoing the previous).

#### Step Length — Backtracking

Start with

$$\alpha' \;=\; \alpha_{n-1}\, \frac{\mathbf{g}_{n-1}^{T}\, \mathbf{d}_{n-1}}{\mathbf{g}_{n}^{T}\, \mathbf{d}_{n}}$$

(for $n = 0$, the user supplies $\alpha_0'$). Then repeatedly set $\alpha \leftarrow \rho\alpha$ (typically $\rho = 1/2$) until the **Armijo rule** holds:

$$f(\mathbf{x}_n + \alpha\, \mathbf{d}_n) \;<\; f(\mathbf{x}_n) + \sigma\,\alpha\, \mathbf{g}_n^{T}\, \mathbf{d}_n, \qquad \sigma \in [10^{-4},\, 10^{-1}]$$

Backtracking is guaranteed to terminate whenever $\mathbf{d}_n$ is a descent direction (by the second-order Taylor remainder).

**Wolfe conditions** (not implemented here): Armijo plus a curvature condition preventing the step from being too short. Ensure both sufficient decrease and non-trivial progress.

### Direction Choices

#### 1. Coordinate Descent

$\mathbf{d}_n$ cycles through $\pm\mathbf{e}_i$. No derivatives needed; $\mathbf{d}_n$ is independent of $f$. Each step cheap; convergence slow.

#### 2. Gradient (Steepest) Descent

$$\mathbf{d}_n \;=\; -\mathbf{g}_n$$

Automatically a descent direction. Argmin of $\mathbf{g}^{T}\mathbf{d}$ over $\lVert \mathbf{d} \rVert \leq 1$ — locally the fastest-decreasing direction.

- Tends to **zigzag** in narrow valleys (when the Hessian is ill-conditioned).
- Converges linearly for well-behaved convex problems, slowly otherwise.

#### 3. Newton's Method

$$\mathbf{d}_n \;=\; -H(f)(\mathbf{x}_n)^{-1}\, \mathbf{g}_n, \qquad \alpha_n = 1 \text{ by default}$$

Derived by minimising the local quadratic approximation of $f$ at $\mathbf{x}_n$. **Solve** $H\, \mathbf{d}_n = -\mathbf{g}_n$ rather than inverting.

- **Convergence:** quadratic near a local minimum with positive-definite Hessian.
- **Cost:** $O(d^{2})$ second-derivative evaluations plus $O(d^{3})$ linear solve per step.
- **Failure modes (important!):**
  - $H$ not positive definite $\to$ $\mathbf{d}_n$ may not be a descent direction; can race towards a maximum or saddle point.
  - $H$ singular $\to$ step undefined.
  - Non-convex regions $\to$ divergence likely.

**Relation to Chapter 5:** Newton-optimisation is exactly Newton-root-finding applied to $\nabla f = \mathbf{0}$. But because zero-gradient alone doesn't distinguish minima from maxima/saddles, optimisation needs extra care (descent check, positive-definite Hessian).

**Why not just minimise $\lVert \mathbf{f}(\mathbf{x})\rVert^{2}$ to solve a root-finding problem?** Because $\lVert \mathbf{f}(\mathbf{x})\rVert^{2}$ is rarely convex; it may have spurious local minima that are not zeros of $\mathbf{f}$. Useful only to generate a starting point for Chapter-5 methods. Exception: when $\mathbf{f}$ is linear (i.e. solving $A\mathbf{x} = \mathbf{b}$), $\lVert \mathbf{f} \rVert^{2}$ is convex and optimisation methods apply cleanly — this is how gradient and conjugate-gradient methods handle large sparse linear systems.

#### 4. BFGS (Quasi-Newton)

Build up a **symmetric positive-definite** approximation $\hat{H}_n$ to the Hessian via rank-2 updates. Let $\Delta\mathbf{x} = \mathbf{x}_n - \mathbf{x}_{n-1}$ and $\Delta\mathbf{g} = \mathbf{g}_n - \mathbf{g}_{n-1}$:

$$\hat{H}_n \;=\; \hat{H}_{n-1} \;+\; \frac{\Delta\mathbf{g}\, \Delta\mathbf{g}^{T}}{\Delta\mathbf{g}^{T}\, \Delta\mathbf{x}} \;-\; \frac{\hat{H}_{n-1}\, \Delta\mathbf{x}\, \Delta\mathbf{x}^{T}\, \hat{H}_{n-1}}{\Delta\mathbf{x}^{T}\, \hat{H}_{n-1}\, \Delta\mathbf{x}}$$

**Curvature condition:** positive-definiteness is preserved only if

$$\Delta\mathbf{x}^{T}\, \Delta\mathbf{g} \;>\; 0$$

For convex $f$ this always holds; for non-convex $f$ it can fail for small step sizes. Simple recovery: **reset $\hat{H}_n = I$** — the next step is then gradient descent, and BFGS can rebuild a good Hessian approximation.

Direction: $\mathbf{d}_n = -\hat{H}_n^{-1}\, \mathbf{g}_n$. Use Sherman–Morrison to track $\hat{H}^{-1}$ directly, avoiding $O(d^{3})$ solves.

- Initial $\hat{H}_0 = I$: the first step is gradient descent. As iterations proceed, $\hat{H}_n$ approaches the true Hessian on well-behaved $f$.
- **Convergence:** superlinear (not quadratic in general); globally convergent under mild conditions.
- **No second derivatives required** — huge savings vs Newton.
- **L-BFGS** is the limited-memory variant, standard in machine learning applications of large dimension.

#### 5. Stochastic / Minibatch Gradient Descent

When $f$ decomposes as a sum $l(\mathbf{w}) = \frac{1}{N}\sum_{i=1}^{N} l_i(\mathbf{w})$ (typical in machine learning):

**Minibatch:** compute the gradient only over a subset $B$:

$$\mathbf{d}_n \;=\; -\frac{1}{\lvert B \rvert}\sum_{i \in B} \nabla l_i(\mathbf{w}_n)$$

**Stochastic (SGD):** one random index per step. Much noisier but vastly cheaper per iteration; can sometimes escape shallow local minima.

### Comparison (Multi-Dimensional)

| Method | Needs | Per-step cost | Convergence | Notes |
|---|---|---|---|---|
| Coordinate descent | $f$ only | $O(d)$ | — | Naive |
| Gradient descent | $\nabla f$ | $O(d)$ | linear | Zigzags; tuning-sensitive |
| Newton | $\nabla f, H$ | $O(d^{3})$ | quadratic near min | Fragile in non-convex |
| BFGS | $\nabla f$ | $O(d^{2})$ | superlinear | Robust, industry default |

### When to Use What — Optimization

- $d = 1$, non-smooth $\to$ golden section or Brent.
- $d = 1$, smooth $\to$ 1D Newton on $f'$ (i.e. Newton-opt).
- Convex, small $d$, cheap Hessian $\to$ Newton with backtracking.
- Non-convex or larger $d$ $\to$ BFGS (or L-BFGS).
- Huge $d$, ML-style sum structure $\to$ SGD or minibatch.
- Large linear system $A\mathbf{x} = \mathbf{b}$ (an implicitly convex quadratic) $\to$ gradient / conjugate-gradient methods.

### Optimization vs Root-Finding — Summary

| | Optimization | Root-Finding |
|---|---|---|
| Best attainable tol | $\sqrt{\varepsilon}$ | $\varepsilon$ |
| "Better" iterate | $f(\mathbf{x}_{n+1}) < f(\mathbf{x}_n)$ | residual smaller (not always monotone) |
| Newton variant | Quadratic near min; fails if $H$ not PD | Quadratic near simple root |
| Second-order structure | Hessian symmetric | Jacobian generally not symmetric |

---

## Required Proofs (ON Syllabus — be able to reproduce)

1. **Composite integration error bounds** (§3.4): derive the composite midpoint, trapezium, and Simpson's bounds

$$\frac{(b-a)^{3}}{24 n^{2}} \max \lvert f'' \rvert, \qquad \frac{(b-a)^{3}}{12 n^{2}} \max \lvert f'' \rvert, \qquad \frac{(b-a)^{5}}{180 n^{4}} \max \lvert f^{(4)} \rvert$$

from the single-strip error formulas, summing over strips and using continuity of derivatives / the intermediate-value theorem.

2. **Monte Carlo convergence**: unbiasedness; $\mathrm{Var}(\mathrm{MC}_{N}) = A(R)^{2}\, \mathrm{Var}(f)/N$; $O(N^{-1/2})$ error via Chebyshev's inequality.

3. **Lemma 4.2**: if $A_n$ has order-$q$ convergence then $A_{2n}$ has order $q^{2}$.

4. **Lemma 5.3 (1D Newton quadratic convergence)**: the identity

$$\epsilon_{n+1} \;=\; -\frac{f''(\xi_n)}{2\, f'(x_n)}\, \epsilon_n^{2}$$

and the consequent quadratic rate.

5. **Stationary-point classification** via Hessian definiteness (from the second-order Taylor expansion).

6. **Convexity results**: local min equals global min for convex $f$; twice-differentiable $f$ is convex iff Hessian is PSD; Jensen's inequality.

7. **Golden-ratio derivation** for golden section search: $\phi^{2} + \phi = 1$.

8. **Lagrange multiplier method** — statement and usage (the proof is bonus; see below).

9. **Convergence of Newton's method for 1D optimization** (quadratic, derived by applying 1D Newton root-finding to $f'$).

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
