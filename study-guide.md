# Continuous Mathematics — Numerical Methods Study Guide

Covers **Integration**, **Root-Finding**, **Optimization**. For each method: formula, applicable dimension, when to use, conditions, error/convergence. Required proofs vs bonus/non-syllabus material are listed at the end.

---

## Foundational Background (Chapters 1–2, 4)

### Derivatives (nD)
- **Gradient**: $\nabla f = \left(\frac{\partial f}{\partial x_1}, \ldots, \frac{\partial f}{\partial x_d}\right)^T$ (column vector).
- **Hessian**: $H(f)_{ij} = \frac{\partial^2 f}{\partial x_i \partial x_j}$ — symmetric when second partials continuous.
- **Jacobian** of $f : \mathbb{R}^n \to \mathbb{R}^m$: $J_{ij} = \frac{\partial f_i}{\partial x_j}$. Note: Hessian of scalar $f$ = Jacobian of $\nabla f$.

### Taylor's Theorem
- **1D, order $n$ (Lagrange remainder)**:

$$f(x+h) = \sum_{k=0}^{n} \frac{h^k}{k!} f^{(k)}(x) + \frac{h^{n+1}}{(n+1)!} f^{(n+1)}(\xi)$$

for some $\xi$ between $x$ and $x+h$.

- **nD, second order**:

$$f(\mathbf{x}+\mathbf{h}) = f(\mathbf{x}) + \mathbf{h}^T \nabla f(\mathbf{x}) + \tfrac{1}{2} \mathbf{h}^T H(f)(\mathbf{x}+\theta \mathbf{h}) \mathbf{h}$$

for some $\theta \in (0,1)$.

### Convexity
- $f$ convex on convex set $S$ iff for all $\mathbf{x}, \mathbf{y} \in S$, $\lambda \in [0,1]$:

$$f(\lambda \mathbf{x} + (1-\lambda)\mathbf{y}) \leq \lambda f(\mathbf{x}) + (1-\lambda) f(\mathbf{y}).$$

- Twice-differentiable $f$ convex $\Leftrightarrow$ Hessian positive semidefinite everywhere. Strictly convex $\Leftarrow$ Hessian positive definite.
- Convex $f$: every local min is global; stationary point $\Rightarrow$ global min.

### Stationary Point Classification (nD)
At stationary point $\mathbf{x}^*$ ($\nabla f = 0$):
- $H(f)(\mathbf{x}^*)$ positive definite $\to$ strict local **minimum**.
- Negative definite $\to$ strict local **maximum**.
- Indefinite $\to$ **saddle point**.
- Semidefinite $\to$ inconclusive (need higher-order test).

Tests for definiteness: all eigenvalues positive / all leading principal minors (pivots) positive / Cholesky succeeds.

### Accuracy (Chapter 4)
- **Error**: $\tilde{u} - u$; **absolute**: $|\tilde{u} - u|$; **relative**: $|\tilde{u} - u|/|u|$.
- **Truncation error**: approximating infinite procedure by finite (e.g., Taylor cutoff, $n$-strip integration).
- **Roundoff error**: finite floating-point representation. NOT the same as truncation.
- **Machine epsilon**: $\varepsilon \approx 1.1 \cdot 10^{-16}$ (double), $\approx 6 \cdot 10^{-8}$ (single).
- **Catastrophic cancellation**: subtracting nearly equal numbers $\to$ huge relative error.
- **Forward error**: $\tilde{f}(x) - f(x)$. **Backward error**: $\tilde{x} - x$ where $\tilde{f}(x) = f(\tilde{x})$.

### Rates of Convergence (Ch 4)
Let $\epsilon_n = x_n - x^*$. As $n \to \infty$:

| Type | Condition |
|---|---|
| Linear (order 1) | $\|\epsilon_{n+1}\|/\|\epsilon_n\| \to a$, $0 < a < 1$ |
| Sublinear | $\|\epsilon_{n+1}\|/\|\epsilon_n\| \to 1$ |
| Logarithmic | sublinear + $\|\epsilon_{n+2} - \epsilon_{n+1}\|/\|\epsilon_{n+1} - \epsilon_n\| \to 1$ |
| Superlinear | $\|\epsilon_{n+1}\|/\|\epsilon_n\| \to 0$ |
| Order-$q$ ($q>1$) | $\|\epsilon_{n+1}\|/\|\epsilon_n\|^q \to a$, $a > 0$ |
| Quadratic | order-2 |

**Lemma 4.2**: If $A_n$ converges order $q$, then $A_{2n}$ converges order $q^2$. (Required — comparing half-speed quadratic vs order-1.5.)

Rule of thumb: linear $\Rightarrow$ adds a constant number of decimals per step; quadratic $\Rightarrow$ doubles correct decimals per step.

### Lagrange Multipliers (equality constraints)
To minimize $f(\mathbf{x})$ subject to $g_i(\mathbf{x}) = 0$ ($i = 1, \ldots, m$): introduce Lagrangian

$$L(\mathbf{x}, \boldsymbol{\lambda}) = f(\mathbf{x}) - \sum_i \lambda_i g_i(\mathbf{x}),$$

solve $\nabla_{\mathbf{x}} L = 0$, $\nabla_{\boldsymbol{\lambda}} L = 0$.

**KKT conditions (inequality constraints): NOT on syllabus this year.**

---

## Part I — Numerical Integration (Chapter 3)

Goal: approximate $\int_a^b f(x)\,dx$. Most rules are **Newton–Cotes**: replace $f$ by an interpolating polynomial on equally-spaced nodes, integrate exactly.

### 1D Methods (single strip width $h = b-a$)

| Method | Nodes | Formula | Exact for polys of deg | Error (single strip) |
|---|---|---|---|---|
| **Midpoint** | 1 ($m = (a+b)/2$) | $(b-a) f(m)$ | $\leq 1$ | $-\dfrac{(b-a)^3}{24} f''(\xi)$ |
| **Trapezium** | 2 ($a, b$) | $\dfrac{b-a}{2}(f(a) + f(b))$ | $\leq 1$ | $\dfrac{(b-a)^3}{12} f''(\xi)$ |
| **Simpson's** | 3 ($a, m, b$) | $\dfrac{b-a}{6}(f(a) + 4f(m) + f(b))$ | $\leq 3$ | $-\dfrac{(b-a)^5}{2880} f^{(4)}(\xi)$ |
| **Boole's** | 5 equally spaced | $\dfrac{b-a}{90}(7f_0 + 32f_1 + 12f_2 + 32f_3 + 7f_4)$ | $\leq 5$ | $O((b-a)^7 f^{(6)})$ |

### Composite Rules ($n$ strips of width $h = (b-a)/n$)

- **Composite Midpoint**: $h \sum_i f(m_i)$; abs error bound $\dfrac{(b-a)^3}{24 n^2} \max |f''|$.
- **Composite Trapezium**: $h \left( \tfrac{1}{2} f(a) + \sum_{i=1}^{n-1} f(x_i) + \tfrac{1}{2} f(b) \right)$; error bound $\dfrac{(b-a)^3}{12 n^2} \max |f''|$.
- **Composite Simpson's** ($n$ even): $\dfrac{h}{3}\bigl(f_0 + 4(f_1 + f_3 + \cdots) + 2(f_2 + f_4 + \cdots) + f_n\bigr)$; error bound $\dfrac{(b-a)^5}{180 n^4} \max |f^{(4)}|$.

**Convergence:** composite midpoint/trapezium are $O(n^{-2})$; Simpson's is $O(n^{-4})$; Boole's $O(n^{-6})$.

### When to Use
- $f$ has limited smoothness (only $C^2$) $\to$ midpoint/trapezium.
- $f$ smooth ($C^4$) $\to$ Simpson's (standard default; best accuracy-per-work ratio).
- Very smooth $\to$ Boole's or higher Newton–Cotes (beware: weights become negative and unstable for very high order; use Gaussian quadrature instead — not in syllabus).
- Endpoint singularities $\to$ midpoint (doesn't evaluate at endpoints).
- Non-smooth / discontinuous derivative $\to$ degrades rates; adapt or split the interval.

### Comparison

| Feature | Midpoint | Trapezium | Simpson's |
|---|---|---|---|
| $f$ evals per strip (composite, shared) | $n$ | $n+1$ | $2n+1$ |
| Rate | $O(n^{-2})$ | $O(n^{-2})$ | $O(n^{-4})$ |
| Needs $f$ at endpoints | No | Yes | Yes |
| Exact on polys deg | $\leq 1$ | $\leq 1$ | $\leq 3$ |

Midpoint error has $\sim\tfrac{1}{2}$ the magnitude of trapezium (and opposite sign) — Simpson's is $(2 \cdot \text{midpoint} + \text{trapezium})/3$ combination.

### Multidimensional Integration

- **Tensor product** (e.g. Simpson's in $d$ dims): cost $O(n^d)$ — **curse of dimensionality**. Error still $O(n^{-4})$ but $n$ points per axis means $n^d$ total.
- **Monte Carlo** is the fallback for large $d$.

### Monte Carlo Integration

- Estimator: $\mathrm{MC}_N = A(R) \cdot \dfrac{1}{N} \sum_{i=1}^N f(X_i)$, where $X_i$ i.i.d. uniform on region $R \subset \mathbb{R}^d$, $A(R) =$ volume.
- Unbiased: $\mathbb{E}[\mathrm{MC}_N] = \int_R f$.
- **Variance**: $\mathrm{Var}(\mathrm{MC}_N) = A(R)^2 \mathrm{Var}(f(X)) / N$.
- **RMS error**: $O(N^{-1/2})$ — **independent of dimension $d$**.
- **Chebyshev's inequality**: $\mathbb{P}\!\left(|\mathrm{MC}_N - I| \geq k \sigma / \sqrt{N}\right) \leq 1/k^2$.
- **Central Limit Theorem**: $\mathrm{MC}_N$ approximately normal for large $N$ $\to$ tighter Gaussian-based confidence intervals.

#### Variance Reduction
- **Stratified sampling**: subdivide $R$, sample independently in each subregion. Variance $\leq$ plain MC (never worse if strata chosen sensibly).
- **Importance sampling**: sample from $p(\mathbf{x})$ proportional to $|f(\mathbf{x})|$, weight by $f(X)/p(X)$. Reduces variance when $p$ mimics $|f|$.
- **Control variates**: use $\tilde{f}$ with known integral, estimate $\int (f - \tilde{f}) + \int \tilde{f}$ via MC on the difference.

### When to Use MC
- $d$ large ($\geq \sim 5$): deterministic grids infeasible.
- Irregular region $R$.
- $f$ only available via simulation.
- Don't need high accuracy (only $\sqrt{N}$ convergence).

---

## Part II — Root-Finding (Chapter 5)

Goal: find $\mathbf{x}^*$ with $f(\mathbf{x}^*) = 0$.

### 1D: Bracketing Method

#### Interval Bisection
- Requires **bracket**: $a < b$ with $f(a) f(b) < 0$ and $f$ continuous.
- Iterate: $m = (a+b)/2$; replace one endpoint so sign condition preserved.
- **Convergence**: linear, $|\epsilon_{n+1}| / |\epsilon_n| \leq 1/2$ (always). Always converges.
- **Error bound**: $|\epsilon_n| \leq (b_0 - a_0) / 2^{n+1}$.
- **When to use**: guaranteed method, no derivative info, robust.

### 1D: Open Methods (no bracket)

#### Newton's Method (1D)
- Formula: $x_{n+1} = x_n - \dfrac{f(x_n)}{f'(x_n)}$.
- **Convergence**: quadratic *if* $f'(x^*) \neq 0$, $f''$ continuous, and $x_0$ sufficiently close.
- Linear only (order 1) at multiple roots ($f'(x^*) = 0$).
- **Pitfalls**: may diverge; can overshoot; needs $f'$.
- Derivation: first-order Taylor expansion of $f$ at $x_n$.

#### Secant Method
- Replace $f'(x_n)$ by finite difference from previous iterate:

$$x_{n+1} = x_n - f(x_n) \cdot \frac{x_n - x_{n-1}}{f(x_n) - f(x_{n-1})}.$$

- **Convergence**: superlinear, order $\varphi = (1+\sqrt{5})/2 \approx 1.618$ (golden ratio).
- **When to use**: $f'$ unavailable or expensive. Costs 1 function eval per step.
- **Pitfall**: denominator can blow up (catastrophic cancellation); no bracket.

#### Halley's Method
- $x_{n+1} = x_n - \dfrac{2 f(x_n) f'(x_n)}{2 f'(x_n)^2 - f(x_n) f''(x_n)}$.
- **Convergence**: cubic (order 3).
- **Pitfall**: needs both $f'$ and $f''$; more work per step.

#### Muller's Method
- Fits a **quadratic** through three iterates ($x_n, x_{n-1}, x_{n-2}$); root of quadratic nearest $x_n$ is $x_{n+1}$.
- **Convergence**: superlinear, order $\approx 1.839$.
- **Use**: complex roots (quadratic formula can return complex values).

#### Inverse Quadratic Interpolation (IQI)
- Fits quadratic to $(f(x_i), x_i)$ — swaps roles of $x$ and $y$ — sets $x_{n+1}$ at the interpolated value where $y = 0$.
- **Order** $\approx 1.839$ (same as Muller).
- Fails if two $f$ values coincide.

#### Brent's Method
- **Hybrid** of bisection, secant, and IQI. Maintains a bracket.
- Accepts superlinear step only if it lies within a "safe" region; else falls back to bisection.
- **Convergence**: superlinear in most cases, guaranteed linear worst-case.
- **When to use**: default robust 1D root-finder (e.g. `zeroin`). Guarantees relative error $\leq 4\varepsilon$.

### Comparison (1D)

| Method | Order | Needs derivatives | Bracket | Guaranteed |
|---|---|---|---|---|
| Bisection | 1 (linear) | No | Yes | Yes |
| Newton | 2 | $f'$ | No | No |
| Secant | 1.618 | No | No | No |
| Halley | 3 | $f', f''$ | No | No |
| Muller / IQI | 1.839 | No | No | No |
| Brent | $\geq 1$ ($\approx 1.839$) | No | Yes | Yes |

### Multi-Dimensional Root-Finding

$\mathbf{f} : \mathbb{R}^d \to \mathbb{R}^d$, want $\mathbf{x}^*$ with $\mathbf{f}(\mathbf{x}^*) = \mathbf{0}$.

#### Newton's Method (nD)
- Linearize: $\mathbf{f}(\mathbf{x}) \approx \mathbf{f}(\mathbf{x}_n) + J(\mathbf{x}_n)(\mathbf{x} - \mathbf{x}_n)$.
- Iteration: **solve** $J(\mathbf{x}_n) \Delta \mathbf{x} = -\mathbf{f}(\mathbf{x}_n)$; then $\mathbf{x}_{n+1} = \mathbf{x}_n + \Delta \mathbf{x}$.
- Do **not** compute $J^{-1}$ explicitly — solve the system (LU, Cholesky, …).
- **Cost**: $O(d^3)$ per step; $d^2$ derivative evaluations.
- **Convergence**: quadratic if $J(\mathbf{x}^*)$ nonsingular and $\mathbf{x}_0$ close.
- Fails if $J(\mathbf{x}_n)$ singular.

#### Broyden's Method (Quasi-Newton)
- Approximates Jacobian with rank-1 updates (avoid computing $d^2$ derivatives).
- Update:

$$\hat{J}_n = \hat{J}_{n-1} + \frac{\Delta \mathbf{f} - \hat{J}_{n-1} \Delta \mathbf{x}}{\|\Delta \mathbf{x}\|^2} \Delta \mathbf{x}^T$$

where $\Delta \mathbf{x} = \mathbf{x}_n - \mathbf{x}_{n-1}$, $\Delta \mathbf{f} = \mathbf{f}(\mathbf{x}_n) - \mathbf{f}(\mathbf{x}_{n-1})$.

- Satisfies **secant equation**: $\hat{J}_n \Delta \mathbf{x} = \Delta \mathbf{f}$.
- **Convergence**: superlinear (not quadratic).
- Sherman–Morrison formula updates $\hat{J}^{-1}$ in $O(d^2)$ $\to$ avoid full linear solve.
- Initial $\hat{J}_0$ typically identity or true $J(\mathbf{x}_0)$.

### Error Analysis — Newton's Method (Lemma 5.3) — REQUIRED
For 1D Newton, assuming $f''$ continuous and $f'(x^*) \neq 0$:

$$\epsilon_{n+1} = -\frac{f''(\xi_n)}{2 f'(x_n)} \epsilon_n^2$$

$\to |\epsilon_{n+1}| / |\epsilon_n|^2 \to |f''(x^*)| / (2|f'(x^*)|)$, establishing **quadratic convergence**.

### A Posteriori Error Estimate
Given $f$ and a candidate $\tilde{x}$, bound how far $\tilde{x}$ is from a root:

$$|\tilde{x} - x^*| \approx \frac{|f(\tilde{x})|}{|f'(\tilde{x})|}$$

(from Taylor; exact for linear $f$).

### Termination Criteria (Root-Finding)
- $|x_n - x_{n-1}| < \text{tol} \cdot (1 + |x_n|)$ (step small).
- $|f(x_n)| < \text{tol}$ (residual / backward error).
- $n = N$ (iteration budget).
- Step ill-defined (singular Jacobian, zero denominator).
- $\text{tol}$ can sensibly be $O(\varepsilon)$ for root-finding (unlike optimization).

---

## Part III — Optimization (Chapters 2, 6)

Goal: find local min of $f : \mathbb{R}^d \to \mathbb{R}$. Global min is hard for non-convex.

### Tolerance Reality Check (Important!)
Near minimum, $f(\mathbf{x}) - f(\mathbf{x}^*) \approx \tfrac{1}{2}(\mathbf{x} - \mathbf{x}^*)^T H (\mathbf{x} - \mathbf{x}^*)$. Thus $|f(\mathbf{x}) - f(\mathbf{x}^*)| / |f(\mathbf{x}^*)| < \varepsilon$ when $|\mathbf{x} - \mathbf{x}^*| / |\mathbf{x}^*| = O(\sqrt{\varepsilon})$.

$\to$ **Use $\text{tol} = \sqrt{\varepsilon} \approx 10^{-8}$ (double precision) — cannot locate a minimum to full machine precision.**

### Termination Criteria
- $\|\mathbf{x}_n - \mathbf{x}_{n-1}\| < \text{tol}(1 + \|\mathbf{x}_n\|)$ (step).
- $|f(\mathbf{x}_n) - f(\mathbf{x}_{n-1})| < \text{tol}(1 + |f(\mathbf{x}_n)|)$ (value; risky — $f$ may creep to $-\infty$).
- $\|\nabla f(\mathbf{x}_n)\| < \text{tol} \cdot (1 + \|\nabla f(\mathbf{x}_0)\|)$ (gradient — good default if $\nabla f$ available).
- $n = N$, or iteration ill-defined.

### 1D Optimization

#### Golden Section Search
- Needs a **bracket** $(a, b, c)$ with $a < b < c$, $f(b) < f(a)$ and $f(b) < f(c)$.
- Pick new point $z$ in larger subinterval at golden ratio $\phi = (\sqrt{5} - 1)/2 \approx 0.618$:
  - $z = a + \phi(b - a)$ or $z = c - \phi(c - b)$.
- Keep smaller bracket (3 of $\{a, b, z, c\}$).
- **Convergence**: linear, ratio $\phi \approx 0.618$ per step (worse than bisection's 0.5 for root-finding, but best possible here).
- **Derivatives-free**; $f$ need only be continuous (and unimodal on bracket).
- **Use when**: $f$ non-differentiable or derivative expensive.

Derivation of $\phi$ comes from requiring: (i) new point splits the larger interval so either new bracket has same length, (ii) the ratio persists across iterations $\to \phi^2 + \phi = 1$.

#### Successive Parabolic Interpolation
- Fit quadratic through $(a, f(a)), (b, f(b)), (c, f(c))$; new $z$ = vertex.
- **Convergence**: superlinear ($\sim$order 1.325).
- **Pitfall**: can fail if points collinear, or vertex is a maximum.

#### Brent's Method (for Minimization)
- Combines parabolic interpolation (when it makes enough progress) and golden section (fallback).
- Default robust 1D minimizer.

### Multi-Dimensional Optimization — Line Search Framework

Iterate: $\mathbf{x}_{n+1} = \mathbf{x}_n + \alpha_n \mathbf{d}_n$, choosing:
- **Direction** $\mathbf{d}_n$: must be **descent direction** $\mathbf{g}_n^T \mathbf{d}_n < 0$ (where $\mathbf{g}_n = \nabla f(\mathbf{x}_n)$).
- **Step length** $\alpha_n > 0$.
- **Infinite travel condition**: $\sum \alpha_n = \infty$ — prevents stalling.

#### Choosing Step Length — Backtracking
Start with $\alpha' = \alpha_{n-1} \cdot (\mathbf{g}_{n-1}^T \mathbf{d}_{n-1}) / (\mathbf{g}_n^T \mathbf{d}_n)$ (or user-supplied for first step).

Repeatedly $\alpha \leftarrow \rho \alpha$ (e.g. $\rho = \tfrac{1}{2}$) until **Armijo rule** holds:

$$f(\mathbf{x}_n + \alpha \mathbf{d}_n) < f(\mathbf{x}_n) + \sigma \alpha \mathbf{g}_n^T \mathbf{d}_n,$$

$\sigma \in [10^{-4}, 10^{-1}]$.

**Wolfe conditions** (mentioned, not implemented): Armijo + curvature condition $|\mathbf{d}^T \nabla f(\mathbf{x} + \alpha \mathbf{d})| \leq c_2 |\mathbf{d}^T \mathbf{g}|$; ensure step neither too short nor too long.

### Direction Choices

#### 1. Coordinate Gradient Descent
- $\mathbf{d}_n$ cycles through $\pm \mathbf{e}_i$. No derivatives needed.
- Naive; independent of $f$. Each step fast.

#### 2. Gradient (Steepest) Descent
- $\mathbf{d}_n = -\mathbf{g}_n$. Guaranteed descent direction.
- Is argmin of $\mathbf{g}^T \mathbf{d}$ subject to $\|\mathbf{d}\| \leq 1$.
- **Convergence**: linear for strongly convex Lipschitz $f$ with step $\alpha = 1/\bar{\lambda}$ (see bonus). Rate $\sqrt{1 - \alpha \underline{\lambda}}$.
- Zigzags on narrow valleys (ill-conditioned Hessian).

#### 3. Newton's Method
- $\mathbf{d}_n = -H(f)(\mathbf{x}_n)^{-1} \mathbf{g}_n$, $\alpha_n = 1$ by default.
- Minimizes local quadratic approximation.
- Solve $H \mathbf{d}_n = -\mathbf{g}_n$ (don't invert).
- **Convergence**: quadratic near a minimum with positive-definite Hessian.
- **Cost**: $O(d^2)$ derivatives + $O(d^3)$ linear solve per step.
- **Fails** when:
  - $H$ not positive definite $\to \mathbf{d}_n$ not a descent direction; may converge to maximum/saddle!
  - $H$ singular.
  - Non-convex regions $\to$ can diverge or race to a maximum.
- Fix (non-syllabus tangent): add $\lambda I$ to $H$ until positive definite (Levenberg / damped Newton).

**Relation to Ch 5**: Newton-opt = Newton-root-finding applied to $\nabla f = 0$. But zero-gradient doesn't distinguish min/max/saddle — motivation for more careful methods here.

#### 4. Conjugate Gradient Descent (mentioned, not examined in detail)
- Polak–Ribière:

$$\mathbf{d}_n = -\mathbf{g}_n + \mathbf{d}_{n-1} \cdot \frac{(\mathbf{g}_n - \mathbf{g}_{n-1})^T \mathbf{g}_n}{\mathbf{g}_{n-1}^T \mathbf{g}_{n-1}}.$$

- Avoids zigzag; uses previous direction info without storing Hessian.

#### 5. BFGS (Quasi-Newton)
Approximate Hessian $\hat{H}_n$ with a rank-2 update preserving symmetry and positive definiteness. Let $\Delta \mathbf{x} = \mathbf{x}_n - \mathbf{x}_{n-1}$, $\Delta \mathbf{g} = \mathbf{g}_n - \mathbf{g}_{n-1}$:

$$\hat{H}_n = \hat{H}_{n-1} + \frac{\Delta \mathbf{g} \Delta \mathbf{g}^T}{\Delta \mathbf{g}^T \Delta \mathbf{x}} - \frac{\hat{H}_{n-1} \Delta \mathbf{x} \Delta \mathbf{x}^T \hat{H}_{n-1}}{\Delta \mathbf{x}^T \hat{H}_{n-1} \Delta \mathbf{x}}.$$

- **Curvature condition**: $\Delta \mathbf{x}^T \Delta \mathbf{g} > 0$. If fails (can happen in non-convex region), **reset $\hat{H}_n = I$** ($\to$ next step is gradient descent).
- $\hat{H}_0 = I$ typically.
- Direction: $\mathbf{d}_n = -\hat{H}_n^{-1} \mathbf{g}_n$. Use Sherman–Morrison to track $\hat{H}^{-1}$ directly (avoid $O(d^3)$).
- **Convergence**: superlinear (not necessarily quadratic). Globally convergent under mild conditions.
- **L-BFGS**: limited-memory variant, popular in ML.
- No $d^2$ derivatives required — huge savings.

### Stochastic / Minibatch Gradient Descent
For $l(\mathbf{w}) = (1/N) \sum_i l_i(\mathbf{w})$:
- **Minibatch**: $\mathbf{d}_n = -\dfrac{1}{|B|} \sum_{i \in B} \nabla l_i(\mathbf{w}_n)$.
- **Stochastic (SGD)**: one random $i$ per step. Many more iterations but each much cheaper. Noisy gradient; can escape shallow local minima.

### Comparison

| Method | Needs | Per-step cost | Convergence | Notes |
|---|---|---|---|---|
| Coord. descent | $f$ only | $O(d)$ | — | Naive |
| Gradient descent | $\nabla f$ | $O(d)$ | linear | Zigzags; needs step-size tuning |
| Newton | $\nabla f, H$ | $O(d^3)$ | quadratic near min | Fragile in non-convex |
| BFGS | $\nabla f$ | $O(d^2)$ | superlinear | Robust, industry standard |
| CG (Polak-Ribière) | $\nabla f$ | $O(d)$ | — | Avoids zigzag |

### When to Use What (Optimization)
- $d = 1$, non-smooth $\to$ Golden section / Brent.
- Convex, small $d$, smooth, have Hessian $\to$ Newton.
- Non-convex or large $d$ $\to$ BFGS / L-BFGS.
- Very large $d$, huge data (ML) $\to$ SGD / minibatch / Adam (Adam not in syllabus).
- Linear least-squares ($\|A\mathbf{x} - \mathbf{b}\|^2$): gradient or CG — always convex.

### Optimization vs Root-Finding

| | Optimization | Root-Finding |
|---|---|---|
| Best accuracy | $\text{tol} = \sqrt{\varepsilon}$ | $\text{tol} = \varepsilon$ |
| "Better" iterate | $f(\mathbf{x}_{n+1}) < f(\mathbf{x}_n)$ | $\|f(\mathbf{x}_{n+1})\| < \|f(\mathbf{x}_n)\|$ (not always) |
| Newton variant | Quadratic near min, fails if $H$ not PD | Quadratic near root |
| Structure | Hessian symmetric | Jacobian not necessarily |

Converting between: $\mathbf{f}(\mathbf{x}) = \mathbf{0}$ $\Rightarrow$ $\min \|\mathbf{f}(\mathbf{x})\|^2$, but $\|\mathbf{f}(\mathbf{x})\|^2$ rarely convex; may have spurious minima. Good for generating starting points only.

---

## Required Proofs (On Syllabus)

Be ready to reproduce these.

1. **Error bounds for composite numerical integration** (Chapter 3, §3.4):
   - Composite midpoint: $|E| \leq \dfrac{(b-a)^3}{24 n^2} \max |f''|$.
   - Composite trapezium: $|E| \leq \dfrac{(b-a)^3}{12 n^2} \max |f''|$.
   - Composite Simpson's: $|E| \leq \dfrac{(b-a)^5}{180 n^4} \max |f^{(4)}|$.

   (Derive from single-strip bounds + sum over strips. Uses Taylor + intermediate value theorem on continuous derivatives.)

2. **Monte Carlo convergence** (Ch 3):
   - Unbiasedness of MC estimator.
   - Variance $\mathrm{Var}(\mathrm{MC}_N) = A(R)^2 \cdot \mathrm{Var}(f) / N$.
   - Error $O(N^{-1/2})$ via Chebyshev's inequality.

3. **Lemma 4.2**: Order-$q$ iteration taken in pairs has order $q^2$. Proof in notes (Ch 4).

4. **Lemma 5.3**: Newton's method 1D quadratic convergence error formula $\epsilon_{n+1} = -\dfrac{f''(\xi)}{2 f'(x_n)} \epsilon_n^2$.

5. **Classification of stationary points** via Hessian definiteness (derived from Taylor 2nd order).

6. **Convexity results**:
   - Convex $f$: local min = global min.
   - Twice-diff $f$ convex iff Hessian PSD.
   - Jensen's inequality (basic).

7. **Golden ratio derivation** for golden section search ($\phi^2 + \phi = 1$).

8. **Lagrange multipliers method** (statement + usage; derivation/proof is bonus — see below).

9. **Convergence of Newton's method for optimization** in 1D (quadratic, via reduction to Newton root-finding on $f'$).

---

## Bonus / NOT Part of Syllabus

Explicitly flagged "not part of the syllabus" in the notes:

1. **§1.8 — Proof of Taylor's Theorem** (nD via 1D construction). Use the statement only.

2. **§2.8 — Proof of the Lagrange Multiplier Method**. Use the method; the proof (geometric/IFT argument) is optional reading.

3. **KKT Conditions** (inequality-constrained optimization). Described for context; "we will not use this method this year".

4. **§3.9 — Some More Error Bounds**: single-strip Trapezium and Simpson's error bound derivations via Rolle/IVT-style arguments. (Composite bounds ARE required; single-strip proofs are bonus.)

5. **§4.5 — Acceleration Methods**:
   - Aitken extrapolation / $\Delta^2$ method.
   - Steffensen's method (quadratic convergence from linear via Aitken).
   - Theorem 4.4 ($\Delta^2$ accelerates linear convergence strictly).

6. **§5.8 — Convergence of Newton's Method in $d$ dimensions** (quadratic convergence proof using Taylor + operator norms). 1D version IS required.

7. **§6.6 — Linear Convergence of Gradient Descent** for strongly convex Lipschitz $f$ (Lemma 6.3 & Theorem 6.4 with condition number argument). Bonus material.

8. **Conjugate Gradient Descent (Polak-Ribière) derivation**. Formula given for context; full derivation beyond syllabus.

9. **Sherman-Morrison formula** update for BFGS inverse Hessian — given as efficient implementation; derivation optional.

10. **Hessian modification schemes** for Newton's method (add $\lambda I$ until PD; Cholesky-based $\lambda$ search).

11. **Forward vs backward error distinction** (mentioned in Ch 4; distinction is "outside our syllabus").

12. **Details of Wolfe conditions** beyond stating them.

13. **Nelder-Mead / simplex method** (mentioned in passing for nD derivative-free opt; not examined).

14. **Brent's method internals** (hybrid logic — concept examined, internal switching rules not).

---

## Quick Reference — "Which Method?"

**Integration**
- Smooth $f$, 1D, default $\to$ Composite Simpson's.
- Rough $f$ or 1D quick $\to$ Trapezium (composite).
- Endpoint singularity $\to$ Midpoint.
- $d \geq 5$ or irregular region $\to$ Monte Carlo (+ variance reduction).

**Root-Finding (1D)**
- Have a bracket + want guarantee $\to$ Bisection or Brent.
- Smooth + good guess + cheap derivative $\to$ Newton.
- Smooth + good guess + no derivative $\to$ Secant.
- Production code $\to$ Brent.

**Root-Finding (nD)**
- Can afford full Jacobian $\to$ Newton.
- Can't $\to$ Broyden.

**Optimization (1D)**
- Non-smooth / derivative-free $\to$ Golden section or Brent.
- Smooth $\to$ 1D Newton on $f'$ (i.e., Newton-opt).

**Optimization (nD)**
- Small $d$, Hessian cheap, convex $\to$ Newton + backtracking.
- Default robust workhorse $\to$ BFGS.
- Huge $d$ + data $\to$ SGD / minibatch.
- Large linear system $A\mathbf{x} = \mathbf{b}$ (implicitly convex quadratic) $\to$ CG.
