# Continuous Mathematics — Summary Tables

*Oxford HT 2026 · Andrew D. Ker · Based on lecture notes + study guide*

---

## Table 1: Root Finding — 1D Methods

**Goal:** Find $x^*$ such that $f(x^*) = 0$, where $f : \mathbb{R} \to \mathbb{R}$.

### Variable Key

| Symbol | Meaning |
|--------|---------|
| $x^*$ | True root: $f(x^*) = 0$ |
| $x_n$ | Current iterate (approximation to $x^*$) |
| $\epsilon_n = x_n - x^*$ | Error at step $n$ |
| $[a, b]$ | Bracket: interval known to contain a root ($f(a)f(b)<0$) |
| $m = (a+b)/2$ | Midpoint of current bracket |
| $f'(x_n)$, $f''(x_n)$ | First and second derivatives of $f$ at $x_n$ |
| tol | Stopping tolerance |

### Method Comparison Table

| Property | **Bisection** | **Newton's** | **Secant** | **Halley's** | **Muller's / IQI** | **Brent's** |
|---|---|---|---|---|---|---|
| **Next iterate formula** | $m = \frac{a+b}{2}$; if $f(a)f(m)<0$ set $b=m$, else set $a=m$ | $x_{n+1} = x_n - \dfrac{f(x_n)}{f'(x_n)}$ | $x_{n+1} = x_n - f(x_n)\cdot\dfrac{x_n - x_{n-1}}{f(x_n)-f(x_{n-1})}$ | $x_{n+1} = x_n - \dfrac{2f(x_n)f'(x_n)}{2[f'(x_n)]^2 - f(x_n)f''(x_n)}$ | **Muller:** fit quadratic through last 3 points, take root closest to $x_n$. **IQI:** fit quadratic to swapped $(f_i, x_i)$ points, evaluate at $f=0$ | Hybrid: attempt IQI/secant step; accept only if inside bracket; else bisect |
| **Requires bracket?** | Yes | No | No | No | No | Yes |
| **Derivatives needed** | None | $f'$ | None | $f'$ and $f''$ | None | None |
| **Convergence order** | Linear $q=1$ | Quadratic $q=2$ | Superlinear $q \approx 1.618$ | Cubic $q=3$ | $q \approx 1.839$ | $q \geq 1$ (typically $\approx 1.839$) |
| **Convergence constant** | Rate $= \frac{1}{2}$ always | $C = \frac{ \vert f''(x^*) \vert }{2 \vert f'(x^*) \vert }$ | — | — | — | — |
| **Error bound** | $ \vert \epsilon_n \vert  \leq \frac{b_0 - a_0}{2^{n+1}}$ | $ \vert \epsilon_{n+1} \vert  \approx C \vert \epsilon_n \vert ^2$ (Lemma 5.3) | — | — | — | Worst case same as bisection |
| **Cost per step** | 1 $f$ eval | 1 $f$ + 1 $f'$ | 1 $f$ (cache previous) | 1 $f$ + 1 $f'$ + 1 $f''$ | 1 $f$ (cache 2 previous) | 1 $f$ |
| **Conditions for convergence** | $f$ continuous; $f(a)f(b)<0$ | $f'(x^*)\neq 0$ (simple root); $f''$ continuous near $x^*$; $x_0$ sufficiently close | $f$ smooth; $x_0, x_1$ close to $x^*$; $f(x_0)\neq f(x_1)$ | $f'(x^*)\neq 0$; $f, f', f''$ available near $x^*$ | $f$ smooth; 3 starting points with distinct $f$-values | $f(a)f(b)<0$; $f$ continuous |
| **Limitations** | Slow (1 binary digit/step); must have bracket | Diverges if $f'(x_n)\approx 0$; fails at multiple roots (degrades to linear); no bracket | Catastrophic cancellation near root; no bracket guarantee | Needs $f''$; expensive per step; can diverge | IQI fails if two $f$-values equal; Muller can produce complex iterates | Switching logic opaque |
| **When to use** | Need guaranteed convergence; non-smooth $f$; robust fallback | Smooth $f$; $f'$ cheap; good initial guess | No $f'$ available; smooth $f$; good initial guess | $f', f''$ both cheap; need maximum speed per step | Complex roots desired; when IQI/Muller's order beats Newton | **Production default**: need robustness + speed |

### Key Required Formula (Lemma 5.3)

Expand $f(x^*) = 0$ in Taylor series around $x_n$:

$$f(x^*) = f(x_n) + f'(x_n)(x^* - x_n) + \tfrac{1}{2}f''(\xi_n)(x^* - x_n)^2 = 0$$

Dividing by $f'(x_n)$ and using $\epsilon_n = x_n - x^*$:

$$\boxed{\epsilon_{n+1} = -\frac{f''(\xi_n)}{2f'(x_n)}\,\epsilon_n^2} \qquad \Longrightarrow \qquad \frac{|\epsilon_{n+1}|}{|\epsilon_n|^2} \to \frac{|f''(x^*)|}{2|f'(x^*)|}$$

This proves **quadratic convergence** of Newton's method at simple roots. $\xi_n$ lies between $x_n$ and $x^*$ (Lagrange remainder).

### A Posteriori Error Estimate

$$|\tilde{x} - x^*| \approx \frac{|f(\tilde{x})|}{|f'(\tilde{x})|}$$

Exact when $f$ is linear. Essentially: the size of one Newton step from $\tilde{x}$.

---

### Worked Examples

**Bisection** — find $\sqrt{2}$, i.e. root of $f(x) = x^2 - 2$. Initial bracket $[1, 2]$ ($f(1)=-1<0$, $f(2)=2>0$).

| $n$ | $a$ | $b$ | $m = (a+b)/2$ | $f(m)$ | New bracket |
|-----|-----|-----|---------------|--------|-------------|
| 0 | 1.000 | 2.000 | 1.500 | $+0.250$ | $[1.000, 1.500]$ (since $f(1)<0$, $f(1.5)>0$) |
| 1 | 1.000 | 1.500 | 1.250 | $-0.438$ | $[1.250, 1.500]$ |
| 2 | 1.250 | 1.500 | 1.375 | $-0.109$ | $[1.375, 1.500]$ |
| 3 | 1.375 | 1.500 | 1.438 | $+0.066$ | $[1.375, 1.438]$ |

Error bound after step $n$: $|\epsilon_n| \leq 1/2^{n+1}$. After 10 steps: error $\leq 0.001$. After 53 steps: double-precision accuracy.

---

**Newton's Method** — root of $f(x) = x^2 - 2$, $f'(x) = 2x$. Start $x_0 = 1$.

$$x_{n+1} = x_n - \frac{x_n^2 - 2}{2x_n} = \frac{x_n}{2} + \frac{1}{x_n}$$

| $n$ | $x_n$ | $f(x_n)$ | $x_{n+1} = x_n - f/f'$ | $ \vert \epsilon_n \vert  =  \vert x_n - \sqrt{2} \vert $ |
|----|-------|---------|----------------------|----------------------------------|
| 0 | 1.0000 | $-1.0000$ | 1.5000 | 0.4142 |
| 1 | 1.5000 | $+0.2500$ | 1.4167 | 0.0858 |
| 2 | 1.4167 | $+0.0069$ | 1.4142 | 0.0025 |
| 3 | 1.41422 | $\approx 5\times10^{-7}$ | 1.41421 | $\approx 6\times10^{-7}$ |

Note errors roughly **square**: $0.0858 \approx 0.4142^2/0.5$, $0.0025 \approx 0.0858^2/3$. Quadratic convergence — digits **double** each step.

---

**Secant Method** — root of $f(x) = x^2 - 2$. Start $x_0 = 1$, $x_1 = 2$; $f_0 = -1$, $f_1 = 2$.

$$x_2 = 2 - 2 \cdot \frac{2-1}{2-(-1)} = 2 - \frac{2}{3} = \tfrac{4}{3} \approx 1.333, \quad f(x_2) \approx -0.222$$

$$x_3 = \tfrac{4}{3} - (-0.222)\cdot\frac{\tfrac{4}{3}-2}{-0.222-2} \approx 1.400, \quad f(x_3) \approx -0.040$$

$$x_4 \approx 1.4138, \quad |\epsilon_4| \approx 4\times10^{-4}$$

**Key:** only **one** $f$-evaluation per step (cache $f(x_1)$ from previous step). Convergence faster than linear but slower than quadratic — order $\varphi = (1+\sqrt{5})/2 \approx 1.618$.

---

**Halley's Method** — root of $f(x) = x^2 - 2$; $f'(x) = 2x$; $f''(x) = 2$. Start $x_0 = 1$.

$$x_1 = 1 - \frac{2(-1)(2)}{2(4) - (-1)(2)} = 1 - \frac{-4}{8+2} = 1 + 0.4 = 1.4$$

$$x_2 = 1.4 - \frac{2(-0.04)(2.8)}{2(7.84)-(-0.04)(2)} = 1.4 - \frac{-0.224}{15.68+0.08} = 1.4 + 0.01417 \approx 1.41417$$

$x_3 \approx 1.41421$ (correct to 5 d.p.) — **cubic convergence**, digits triple each step.

---

**Muller's Method** — root of $f(x) = x^2 - 2$. Three starting points: $x_0=0$, $x_1=1$, $x_2=2$.

Fit quadratic $p(x) = ax^2 + bx + c$ through $(0,-2),(1,-1),(2,2)$:

$c = -2$; $a+b-2 = -1 \Rightarrow a+b=1$; $4a+2b-2=2 \Rightarrow 2a+b=2$. So $a=1, b=0, c=-2$.

$p(x) = x^2 - 2$. Roots: $x = \pm\sqrt{2}$. Take positive: $x_3 = \sqrt{2} \approx 1.414$. **Exact in one step** here because $f$ is already quadratic!

**Key feature:** Muller's quadratic formula can give complex roots — if the fit gives $b^2-4ac < 0$, the method finds **complex roots** even starting from real iterates. No other method in this syllabus can do this.

---

**IQI (Inverse Quadratic Interpolation)** — same points: $x_0=0, f_0=-2$; $x_1=1, f_1=-1$; $x_2=2, f_2=2$.

Fit $x = A f^2 + Bf + C$ (quadratic in $f$), evaluate at $f=0$ using Lagrange form:

$$x_3 = x_0\frac{(0-f_1)(0-f_2)}{(f_0-f_1)(f_0-f_2)} + x_1\frac{(0-f_0)(0-f_2)}{(f_1-f_0)(f_1-f_2)} + x_2\frac{(0-f_0)(0-f_1)}{(f_2-f_0)(f_2-f_1)}$$

$$= 0\cdot\frac{(1)(-2)}{(-1)(-4)} + 1\cdot\frac{(2)(-2)}{(1)(-3)} + 2\cdot\frac{(2)(1)}{(4)(3)} = 0 + \frac{-4}{-3} + \frac{4}{12} = \frac{4}{3} + \frac{1}{3} = \frac{5}{3} \approx 1.667$$

**Pitfall:** if $f_0 = f_1$ (or any two $f$-values equal), the denominator vanishes — IQI undefined.

---

**Brent's Method** — maintains bracket $[a,b]$ (with $f(a)f(b)<0$) and tracks best point $s$. At each step: attempt IQI (or secant) step; accept it only if (i) it falls strictly inside $[a,b]$ and (ii) it moves by at least half the previous step — otherwise fall back to bisection. Guarantees convergence; typically superlinear. **Always the default 1D root-finder in practice** (implements `zeroin`).

---

## Table 2: Root Finding — $n$D Methods

**Goal:** Find $\mathbf{x}^* \in \mathbb{R}^d$ such that $\mathbf{f}(\mathbf{x}^*) = \mathbf{0}$, where $\mathbf{f} : \mathbb{R}^d \to \mathbb{R}^d$ ($d$ equations, $d$ unknowns).

### Variable Key

| Symbol | Meaning |
|--------|---------|
| $\mathbf{x}_n \in \mathbb{R}^d$ | Current iterate (vector) |
| $\mathbf{f}(\mathbf{x}_n) \in \mathbb{R}^d$ | Residual vector at current iterate |
| $J(\mathbf{x}_n) \in \mathbb{R}^{d\times d}$ | Jacobian matrix: $J_{ij} = \partial f_i / \partial x_j$ |
| $\hat{J}_n \in \mathbb{R}^{d\times d}$ | Broyden's Jacobian approximation |
| $\Delta\mathbf{x} = \mathbf{x}_n - \mathbf{x}_{n-1}$ | Step vector taken in last iteration |
| $\Delta\mathbf{f} = \mathbf{f}(\mathbf{x}_n) - \mathbf{f}(\mathbf{x}_{n-1})$ | Change in residual over last step |
| $\Delta\mathbf{x}$ in Broyden update | Same as step vector |

### Method Comparison Table

| Property | **Newton's Method ($n$D)** | **Broyden's Method (Quasi-Newton)** |
|---|---|---|
| **Step formula** | Solve $J(\mathbf{x}_n)\,\Delta\mathbf{x} = -\mathbf{f}(\mathbf{x}_n)$; then $\mathbf{x}_{n+1} = \mathbf{x}_n + \Delta\mathbf{x}$ | Solve $\hat{J}_n\,\Delta\mathbf{x} = -\mathbf{f}(\mathbf{x}_n)$; then $\mathbf{x}_{n+1} = \mathbf{x}_n + \Delta\mathbf{x}$ |
| **Jacobian update** | Recompute full $J(\mathbf{x}_n)$ at every step ($d^2$ partial derivatives) | Rank-1 update: $\hat{J}_n = \hat{J}_{n-1} + \dfrac{(\Delta\mathbf{f} - \hat{J}_{n-1}\Delta\mathbf{x})\,\Delta\mathbf{x}^T}{ \Vert \Delta\mathbf{x} \Vert ^2}$ |
| **Initialisation** | Exact Jacobian $J(\mathbf{x}_0)$ | $\hat{J}_0 = I$ (identity); first step is then $-\mathbf{f}(\mathbf{x}_0)$ |
| **Cost per step** | $O(d^2)$ partial derivatives + $O(d^3)$ LU solve | $O(d^2)$ rank-1 update + $O(d^2)$ Sherman–Morrison inverse update (avoids $O(d^3)$ solve) |
| **Convergence** | Quadratic (near $\mathbf{x}^*$ if $J(\mathbf{x}^*)$ nonsingular, $\mathbf{x}_0$ close) | Superlinear (not quadratic); curvature accumulates over iterations |
| **Conditions** | $J(\mathbf{x}^*)$ nonsingular; $\mathbf{f}$ smooth; $\mathbf{x}_0$ sufficiently close | $\mathbf{f}$ smooth; $\hat{J}_0$ invertible; $\Delta\mathbf{x} \neq \mathbf{0}$ each step |
| **Limitations** | $J$ singular or ill-conditioned $\Rightarrow$ step undefined or unstable; $O(d^3)$ solve expensive for large $d$ | Approximate Jacobian may degrade; requires reset ($\hat{J}=I$) periodically; slower per-iterate convergence |
| **When to use** | Jacobian is cheap to compute; need fast (quadratic) convergence | Jacobian expensive or unavailable; large $d$; willing to trade convergence rate for cheaper steps |

> **Implementation note:** Never form $J^{-1}$ explicitly. Always **solve** $J\,\Delta\mathbf{x} = -\mathbf{f}$ using LU decomposition. Computing $J^{-1}$ costs the same $O(d^3)$ but with larger constants and worse numerical stability.

> **Secant equation:** The Broyden update satisfies $\hat{J}_n\,\Delta\mathbf{x} = \Delta\mathbf{f}$ — the minimal-change update to $\hat{J}_{n-1}$ consistent with the observed $({\Delta\mathbf{x}}, {\Delta\mathbf{f}})$ pair. This is the $n$D analogue of the secant condition.

---

### Worked Examples

**Newton's Method ($n$D)** — solve $\mathbf{f}(x,y) = \begin{pmatrix} x^2 - y \\ y^2 - x \end{pmatrix} = \mathbf{0}$. Solutions: $(0,0)$ and $(1,1)$. Start $\mathbf{x}_0 = (2,2)^T$.

**Jacobian:** $J(x,y) = \begin{pmatrix}2x & -1 \\ -1 & 2y\end{pmatrix}$

**Step 0:** $J(2,2) = \begin{pmatrix}4&-1\\-1&4\end{pmatrix}$, $\mathbf{f}(2,2) = \begin{pmatrix}2\\2\end{pmatrix}$.

Solve $\begin{pmatrix}4&-1\\-1&4\end{pmatrix}\Delta\mathbf{x} = \begin{pmatrix}-2\\-2\end{pmatrix}$:

Row equations: $4\Delta x - \Delta y = -2$ and $-\Delta x + 4\Delta y = -2$. Subtract: $5(\Delta x - \Delta y)=0$, so $\Delta x = \Delta y$. Then $3\Delta x = -2 \Rightarrow \Delta x = -\tfrac{2}{3}$.

$$\mathbf{x}_1 = (2,2)^T + (-\tfrac{2}{3}, -\tfrac{2}{3})^T = (\tfrac{4}{3},\, \tfrac{4}{3})^T \approx (1.333, 1.333)^T$$

**Step 1:** $\mathbf{f}(4/3,4/3)=\begin{pmatrix}16/9-4/3\\16/9-4/3\end{pmatrix}=\begin{pmatrix}4/9\\4/9\end{pmatrix}$. Similar calculation gives $\mathbf{x}_2 \approx (1.067,1.067)^T$, then $\mathbf{x}_3 \approx (1.002,1.002)^T$ — **quadratic convergence** toward $(1,1)$.

---

**Broyden's Method ($n$D)** — same system, $\mathbf{x}_0 = (2,2)^T$, $\hat{J}_0 = I$.

**Step 0:** $\Delta\mathbf{x} = -\hat{J}_0^{-1}\mathbf{f}(\mathbf{x}_0) = -(2,2)^T = (-2,-2)^T$. So $\mathbf{x}_1 = (0,0)^T$, which happens to be a root!

**Demonstrating the update** with $\mathbf{x}_0=(1.5,1.5)^T$, $\hat{J}_0=I$, $\mathbf{f}(1.5,1.5)=(0.75,0.75)^T$:

$\Delta\mathbf{x} = -(0.75,0.75)^T = (-0.75,-0.75)^T$. $\mathbf{x}_1 = (0.75, 0.75)^T$. $\mathbf{f}(0.75,0.75)=(-0.1875,-0.1875)^T$.

$\Delta\mathbf{f} = (-0.1875,\!-0.1875)^T - (0.75,0.75)^T = (-0.9375,-0.9375)^T$.

$\hat{J}_0\Delta\mathbf{x} = I(-0.75,-0.75)^T = (-0.75,-0.75)^T$.

Numerator: $\Delta\mathbf{f} - \hat{J}_0\Delta\mathbf{x} = (-0.1875,-0.1875)^T$. $\|\Delta\mathbf{x}\|^2 = 1.125$.

$$\hat{J}_1 = I + \frac{(-0.1875,-0.1875)^T (-0.75,-0.75)}{1.125} = I + \frac{1}{1.125}\begin{pmatrix}0.1406 & 0.1406 \\ 0.1406 & 0.1406\end{pmatrix}$$

$\hat{J}_1 \approx \begin{pmatrix}1.125 & 0.125 \\ 0.125 & 1.125\end{pmatrix}$ — a rank-1 update that incorporates information about $\mathbf{f}$'s change without computing any new partial derivatives.

---

## Table 3: Optimization — 1D Methods

**Goal:** Find local minimum $x^* = \arg\min_{x \in \mathbb{R}} f(x)$ where $f'(x^*)=0$, $f''(x^*)>0$.

### Variable Key

| Symbol | Meaning |
|--------|---------|
| $(a, b, c)$ | Bracket: $a < b < c$; $f(b) < f(a)$ and $f(b) < f(c)$; minimum lies somewhere in $(a,c)$ |
| $b$ | Current best point (lowest known $f$ value) |
| $z$ | New evaluation point at each step |
| $\phi = \frac{\sqrt{5}-1}{2} \approx 0.618$ | Golden ratio (conjugate); satisfies $\phi^2 + \phi = 1$ |
| $f'(x)$, $f''(x)$ | First and second derivatives |
| $x^*$ | True local minimiser |

### Method Comparison Table

| Property | **Golden Section Search** | **Successive Parabolic Interpolation (SPI)** | **Brent's (1D min)** | **Newton on $f'$** |
|---|---|---|---|---|
| **Core idea** | Probe the larger sub-interval; shrink bracket by factor $\phi$ | Fit quadratic through 3 bracket points; jump to its vertex | Hybrid: SPI if safe (inside bracket, not too large); else golden section | Apply 1D Newton root-finding to $f'(x)=0$ |
| **New point formula** | $z = a + \phi(b-a)$ or $z = c - \phi(c-b)$ (in larger of $[a,b]$ or $[b,c]$) | $z = $ vertex of parabola through $(a,f(a)),(b,f(b)),(c,f(c))$ | SPI step accepted if $z \in (a,c)$ and step not too large; else golden section | $x_{n+1} = x_n - \dfrac{f'(x_n)}{f''(x_n)}$ |
| **Requires bracket?** | Yes: $(a,b,c)$ | Yes: $(a,b,c)$ | Yes | No |
| **Derivatives needed** | None | None | None | $f'$ and $f''$ |
| **Convergence order** | Linear, rate $\phi \approx 0.618$ | Superlinear, $q \approx 1.325$ | Superlinear typical; linear worst case | Quadratic (applying Lemma 5.3 to $f'$) |
| **Cost per step** | 1 $f$ eval (reuse 1 point from previous step) | 1 $f$ eval (reuse 2 points) | 1 $f$ eval | 1 $f'$ + 1 $f''$ |
| **Conditions** | $f$ continuous and unimodal on bracket | 3 non-collinear points; fitted parabola must open upward | Same as golden section; $f$ unimodal | $f''(x^*) > 0$; $f'''$ continuous near $x^*$ |
| **Limitations** | Slightly slower than bisection (rate 0.618 vs 0.5) | Fails if 3 points collinear or parabola opens down | Complex switching logic | Needs $f''$; fails at inflections; can converge to max not min |
| **When to use** | $f$ non-differentiable or derivatives expensive; need guaranteed shrinkage | Fast convergence when parabola is a good local model | **Default 1D minimiser** | $f'', f'$ cheap and $f$ strongly convex locally |

**Derivation of $\phi$:** We need the same proportional bracket reduction whether the new bracket is the left or right sub-interval, **and** the evaluation point from one step should be reusable as an interior point in the next. This self-similarity gives:

$$\phi = \frac{\text{smaller sub-interval}}{\text{total interval}}, \quad (1-\phi) = \phi^2 \quad \Longrightarrow \quad \phi^2 + \phi = 1 \quad \Longrightarrow \quad \phi = \frac{\sqrt{5}-1}{2}$$

---

### Worked Examples

**Golden Section Search** — minimise $f(x) = (x-1)^2$. Bracket: $(a,b,c) = (0,1,3)$.

$f(0)=1$, $f(1)=0$, $f(3)=4$. Sub-intervals: $[0,1]$ width 1, $[1,3]$ width 2. Larger: $[1,3]$.

$$z = 3 - \phi(3-1) = 3 - 0.618\times 2 = 1.764, \quad f(1.764) = (0.764)^2 = 0.584$$

Since $f(z)=0.584 > f(b)=0$: $z$ becomes the new right boundary. New bracket: $(0,\, 1,\, 1.764)$.

Bracket width: $1.764 - 0 = 1.764 = \phi \times 3$. ✓ Old point $b=1$ is **reused** in next step.

| Step | $a$ | $b$ | $c$ | Width | New $z$ from larger sub-interval |
|------|-----|-----|-----|-------|----------------------------------|
| 0 | 0 | 1 | 3 | 3 | 1.764 (right sub-interval is larger) |
| 1 | 0 | 1 | 1.764 | 1.764 | 0.764 (now left sub-interval is larger) |
| 2 | 0.764 | 1 | 1.764 | 1.000 | 1.236 (both equal; pick either) |

After $n$ steps: bracket width $= 3 \times \phi^n$. To reach width $< 0.01$: need $n > \log(300)/\log(1/\phi) \approx 11$ steps.

---

**SPI (Successive Parabolic Interpolation)** — minimise $f(x) = (x-2)^2$. Points: $(0,4),(1,1),(3,1)$.

Fit quadratic $p(x) = Ax^2 + Bx + C$:

$C=4$; $A+B+C=1 \Rightarrow A+B=-3$; $9A+3B+C=1 \Rightarrow 3A+B=-1$. Subtract: $2A=2$, $A=1$, $B=-4$, $C=4$.

$p(x) = x^2 - 4x + 4 = (x-2)^2$. Vertex at $z = -B/(2A) = 4/2 = \mathbf{2}$. **Exact minimum in one step** (SPI is exact when $f$ is itself a quadratic).

For non-quadratic $f$: SPI converges superlinearly ($q\approx1.325$) because each step reuses two previous evaluations and the parabola is a better local model than a straight line.

---

**Brent's Method (1D min)** — hybrid. At each step:
1. Fit parabola to the 3 current best points, compute vertex $z$.
2. **Accept** $z$ if: $z$ is inside $(a,c)$ **and** the step is not too large (no more than half the step from 2 iterations ago).
3. Otherwise: **golden section step** (safe shrinkage).

This prevents the parabola from jumping outside the bracket or stagnating. **Default 1D minimiser** in most software (e.g. `scipy.optimize.minimize_scalar`).

---

**Newton on $f'$** — minimise $f(x) = x^4 - 4x^2$. $f'(x) = 4x^3 - 8x$, $f''(x) = 12x^2 - 8$. Start $x_0 = 2.5$.

$$x_1 = 2.5 - \frac{f'(2.5)}{f''(2.5)} = 2.5 - \frac{4(15.625)-20}{12(6.25)-8} = 2.5 - \frac{62.5-20}{75-8} = 2.5 - \frac{42.5}{67} \approx 1.866$$

$$x_2 = 1.866 - \frac{4(1.866)^3 - 8(1.866)}{12(1.866)^2-8} \approx 1.866 - \frac{26.05-14.93}{41.79-8} = 1.866 - \frac{11.12}{33.79} \approx 1.537$$

$$x_3 \approx 1.426, \quad x_4 \approx 1.415, \quad x_5 \approx 1.4142 \approx \sqrt{2}$$

Converging to local min at $x^* = \sqrt{2}$ (where $f'(\sqrt{2})=0$, $f''(\sqrt{2})=24-8=16>0$). Convergence is quadratic once close.

---

## Table 4: Optimization — $n$D Methods

**Goal:** Find local minimum $\mathbf{x}^* = \arg\min_{\mathbf{x}\in\mathbb{R}^d} f(\mathbf{x})$ where $\nabla f(\mathbf{x}^*)=\mathbf{0}$, $H(\mathbf{x}^*)$ positive definite.

### General Line-Search Framework

$$\mathbf{x}_{n+1} = \mathbf{x}_n + \alpha_n\,\mathbf{d}_n$$

Requirements: (1) $\mathbf{d}_n$ is a **descent direction**: $\mathbf{g}_n^T\mathbf{d}_n < 0$. (2) Step $\alpha_n > 0$. (3) **Infinite travel**: $\sum_n \alpha_n = \infty$.

**Armijo (sufficient decrease) condition** for accepting $\alpha$:
$$f(\mathbf{x}_n + \alpha\,\mathbf{d}_n) < f(\mathbf{x}_n) + \sigma\alpha\,\mathbf{g}_n^T\mathbf{d}_n, \qquad \sigma \in [10^{-4}, 10^{-1}]$$

**Backtracking:** start with $\alpha' = \alpha_{n-1}\frac{\mathbf{g}_{n-1}^T\mathbf{d}_{n-1}}{\mathbf{g}_n^T\mathbf{d}_n}$; repeatedly shrink $\alpha \leftarrow \frac{1}{2}\alpha$ until Armijo holds.

### Variable Key

| Symbol | Meaning |
|--------|---------|
| $\mathbf{x}_n \in \mathbb{R}^d$ | Current iterate |
| $\mathbf{g}_n = \nabla f(\mathbf{x}_n) \in \mathbb{R}^d$ | Gradient at current iterate |
| $H(\mathbf{x}_n) \in \mathbb{R}^{d\times d}$ | Hessian (matrix of second partials) |
| $\hat{H}_n \in \mathbb{R}^{d\times d}$ | BFGS approximation to Hessian |
| $\mathbf{d}_n \in \mathbb{R}^d$ | Search direction |
| $\alpha_n > 0$ | Step length (scalar) |
| $\Delta\mathbf{x} = \mathbf{x}_n - \mathbf{x}_{n-1}$ | Last step vector |
| $\Delta\mathbf{g} = \mathbf{g}_n - \mathbf{g}_{n-1}$ | Change in gradient |
| $\mathbf{e}_i$ | $i$-th standard basis vector (for coordinate descent) |
| $B \subseteq \{1,\ldots,N\}$ | Minibatch index set |

### Method Comparison Table

| Property | **Coordinate Descent** | **Gradient (Steepest) Descent** | **Newton's ($n$D)** | **BFGS** | **SGD / Minibatch** |
|---|---|---|---|---|---|
| **Direction formula** | $\mathbf{d}_n = \pm\mathbf{e}_i$ (cycle through coordinate axes) | $\mathbf{d}_n = -\mathbf{g}_n$ | Solve $H(\mathbf{x}_n)\mathbf{d}_n = -\mathbf{g}_n$ | $\mathbf{d}_n = -\hat{H}_n^{-1}\mathbf{g}_n$ | $\mathbf{d}_n = -\frac{1}{ \vert B \vert }\sum_{i\in B}\nabla l_i(\mathbf{x}_n)$ |
| **Hessian/approximation update** | N/A | N/A | Recompute exact $H(\mathbf{x}_n)$ each step | $\hat{H}_n = \hat{H}_{n-1} + \dfrac{\Delta\mathbf{g}\,\Delta\mathbf{g}^T}{\Delta\mathbf{g}^T\Delta\mathbf{x}} - \dfrac{\hat{H}_{n-1}\Delta\mathbf{x}\,\Delta\mathbf{x}^T\hat{H}_{n-1}}{\Delta\mathbf{x}^T\hat{H}_{n-1}\Delta\mathbf{x}}$ | N/A |
| **Initial Hessian** | N/A | N/A | Exact $H(\mathbf{x}_0)$ | $\hat{H}_0 = I$ (first step = gradient descent) | N/A |
| **Derivatives needed** | None | $\nabla f$ | $\nabla f$ and $H$ (full Hessian) | $\nabla f$ only | $\nabla l_i$ for sampled $i\in B$ |
| **Cost per step** | $O(d)$ + 1D exact minimisation | $O(d)$ gradient + line search | $O(d^2)$ Hessian entries + $O(d^3)$ LU solve | $O(d^2)$ update via Sherman–Morrison–Woodbury | $O(d\cdot \vert B \vert )$ per step |
| **Convergence rate** | Slow (depends on problem structure) | Linear; rate $\approx \frac{\kappa-1}{\kappa+1}$ where $\kappa = \lambda_{\max}/\lambda_{\min}$ | Quadratic near min (when $H$ PD) | Superlinear; globally convergent (mild conditions + line search) | Sublinear in expectation; noisy |
| **Conditions** | $f$ differentiable (or subdifferentiable) | $f$ differentiable; $\nabla f$ Lipschitz | $H(\mathbf{x}^*)$ positive definite; $\mathbf{x}_0$ close to $\mathbf{x}^*$ | Curvature condition: $\Delta\mathbf{x}^T\Delta\mathbf{g}>0$; $\hat{H}_0$ PD | $f = (1/N)\sum_i l_i$; each $\nabla l_i$ computable |
| **Limitations** | Very slow for correlated variables; gets stuck on non-axis-aligned valleys | Zigzags in ill-conditioned (elongated) valleys; consecutive gradients orthogonal (with exact line search) | Fails if $H$ not PD (may race toward saddle/max); $O(d^3)$ expensive; fragile in non-convex regions | Curvature condition can fail for non-convex $f$ (fix: reset $\hat{H}=I$); less accurate than Newton | High variance; needs step-size schedule; not for high-precision problems |
| **When to use** | Separable or nearly-separable $f$; very simple fallback | Convex, well-conditioned $f$; quick implementation | Small $d$; convex $f$; Hessian cheap; need quadratic speed | **Default $n$D optimiser**; large $d$; Hessian unavailable | $d$ huge; $N$ very large; ML sum-structured objectives |

**BFGS curvature condition:** $\Delta\mathbf{x}^T\Delta\mathbf{g} > 0$ ensures the rank-2 update preserves positive-definiteness of $\hat{H}_n$.
- Always holds for convex $f$ (gradient is monotone).
- Can fail for non-convex $f$. **Fix:** reset $\hat{H}_n = I$ and continue.

**L-BFGS:** variant storing only the last $m$ update pairs $(\Delta\mathbf{x}, \Delta\mathbf{g})$ instead of the full $d\times d$ matrix — reduces storage from $O(d^2)$ to $O(md)$. Standard in large-scale ML.

---

### Worked Examples

**Gradient Descent** — minimise $f(x,y) = x^2 + 4y^2$. Gradient: $\mathbf{g} = (2x,\, 8y)^T$.

Start $\mathbf{x}_0 = (3, 1)^T$. $\mathbf{g}_0 = (6, 8)^T$. Direction: $\mathbf{d}_0 = -(6, 8)^T$.

**Exact line search** (minimise $f(\mathbf{x}_0 + \alpha\mathbf{d}_0)$ over $\alpha \geq 0$): For quadratic $f = \frac{1}{2}\mathbf{x}^T A\mathbf{x}$ with $A = \begin{pmatrix}2&0\\0&8\end{pmatrix}$:

$$\alpha^* = \frac{\mathbf{g}_0^T\mathbf{g}_0}{\mathbf{g}_0^T A\mathbf{g}_0} = \frac{36+64}{36\cdot2+64\cdot8} = \frac{100}{72+512} = \frac{100}{584} \approx 0.171$$

Wait — $A$ here acts as the coefficient not exactly the Hessian. For $f = x^2+4y^2$, Hessian $H=\begin{pmatrix}2&0\\0&8\end{pmatrix}$.

$\alpha^* = \frac{\|\mathbf{g}_0\|^2}{\mathbf{g}_0^T H \mathbf{g}_0} = \frac{100}{6^2\cdot2 + 8^2\cdot8} = \frac{100}{72+512} \approx 0.0171$

$\mathbf{x}_1 = (3-6\times0.0171,\; 1-8\times0.0171) = (2.897,\; 0.863)^T$

$\mathbf{g}_1 = (5.79, 6.90)^T$. Check: $\mathbf{g}_1^T\mathbf{d}_0 = 5.79(-6)+6.90(-8) = -34.7-55.2 \neq 0$... let me recalculate.

$\mathbf{x}_1 = (3-6\alpha, 1-8\alpha)$ with $\alpha\approx0.0171$: $x_1 \approx 2.897$, $y_1 \approx 0.863$.
$\mathbf{g}_1 = (2\times2.897, 8\times0.863) = (5.79, 6.90)^T$.
$\mathbf{g}_1^T\mathbf{d}_0 = (5.79)(-6)+(6.90)(-8) = -34.74-55.20 = -89.94 \neq 0$.

Key geometric fact (exact line search): consecutive gradients are **orthogonal** — $\mathbf{g}_{n+1}^T\mathbf{d}_n = 0$. Let's verify: $\mathbf{g}_1^T(-\mathbf{g}_0)=0$ iff $\mathbf{g}_1^T\mathbf{g}_0=0$. Indeed $5.79\times6+6.90\times8 = 34.7+55.2=89.9 \approx 0$? No...

Actually: $\frac{d}{d\alpha}f(\mathbf{x}_0+\alpha\mathbf{d}_0)|_{\alpha=\alpha^*} = \nabla f(\mathbf{x}_1)^T\mathbf{d}_0 = 0$. This is the **gradient orthogonality** with respect to $\mathbf{d}_0$, not $\mathbf{g}_0$: $\mathbf{g}_1^T\mathbf{d}_0 = 0$. So $\mathbf{g}_1\perp\mathbf{d}_0=(-\mathbf{g}_0)$.

$\mathbf{g}_1^T\mathbf{d}_0 = (5.79)(-6)+(6.90)(-8) = -34.74-55.20 = -89.94$. Not zero.

Let me recompute $\alpha^*$ carefully. $f(\mathbf{x}_0+\alpha\mathbf{d}_0) = (3-6\alpha)^2 + 4(1-8\alpha)^2$

$= 9-36\alpha+36\alpha^2 + 4(1-16\alpha+64\alpha^2)$

$= 9-36\alpha+36\alpha^2+4-64\alpha+256\alpha^2$

$= 13-100\alpha+292\alpha^2$

Minimise: $\frac{d}{d\alpha} = -100+584\alpha = 0 \Rightarrow \alpha^* = 100/584 \approx 0.1712$.

$\mathbf{x}_1 = (3-6(0.1712), 1-8(0.1712)) = (3-1.027, 1-1.370) = (1.973, -0.370)^T$.

$\mathbf{g}_1 = (2(1.973), 8(-0.370)) = (3.945, -2.959)^T$.

$\mathbf{g}_1^T\mathbf{d}_0 = (3.945)(-6)+(-2.959)(-8) = -23.67+23.67 = 0$. ✓ Orthogonality confirmed.

Now direction $\mathbf{d}_1 = -(3.945, -2.959)^T = (-3.945, 2.959)^T$. This direction is **not aligned** toward origin $(0,0)$ — it zigzags! This is the **zigzag phenomenon**: steepest descent overshoots along the valley when it's narrow (condition number of $H$ is $\kappa=8/2=4$ here).

**Contrast with Newton:** $\mathbf{d}_0 = -H^{-1}\mathbf{g}_0 = -\begin{pmatrix}1/2&0\\0&1/8\end{pmatrix}(6,8)^T = (-3,-1)^T$.

$\mathbf{x}_1 = (3,1)^T + (-3,-1)^T = (0,0)^T$. **One step to the exact minimum**, regardless of conditioning.

---

**Newton's Method ($n$D opt)** — minimise $f(x,y) = x^2+4y^2$ starting at $(3,1)^T$. $H=\begin{pmatrix}2&0\\0&8\end{pmatrix}$, $\mathbf{g}=(6,8)^T$.

Solve $H\mathbf{d}_0 = -\mathbf{g}_0$: $\begin{pmatrix}2&0\\0&8\end{pmatrix}\mathbf{d}_0 = \begin{pmatrix}-6\\-8\end{pmatrix} \Rightarrow \mathbf{d}_0 = (-3,-1)^T$.

$\mathbf{x}_1 = (0,0)^T$. $\nabla f(0,0) = \mathbf{0}$. Done in **one step**.

For non-quadratic $f$: Newton uses $\alpha_n=1$ by default. If $\alpha=1$ fails the Armijo condition (e.g. $H$ has small or negative eigenvalues), fall back to line search.

---

**BFGS** — same $f$, $\mathbf{x}_0=(3,1)^T$, $\hat{H}_0=I$.

**Step 0:** $\mathbf{d}_0 = -I^{-1}(6,8)^T = (-6,-8)^T$. Line search gives $\alpha_0\approx0.1712$ (same as gradient descent).

$\mathbf{x}_1 = (1.973, -0.370)^T$. $\mathbf{g}_1 = (3.945, -2.959)^T$.

$\Delta\mathbf{x} = (-1.027, -1.370)^T$. $\Delta\mathbf{g} = (-2.055, -10.959)^T$.

Curvature check: $\Delta\mathbf{x}^T\Delta\mathbf{g} = (-1.027)(-2.055)+(-1.370)(-10.959) = 2.110+15.014 = 17.12 > 0$. ✓

BFGS rank-2 update produces $\hat{H}_1$ incorporating curvature information. After sufficient steps, $\hat{H}_n\approx H$ and the method takes near-Newton steps — achieving **superlinear convergence** without ever computing the Hessian explicitly.

---

**SGD / Minibatch** — minimise regression loss $l(w) = \frac{1}{4}\sum_{i=1}^4(wx_i-y_i)^2$ with data $(x_i,y_i)$: $(1,2),(2,4),(3,5),(4,8)$. True slope $w^*\approx 1.9$.

Start $w_0=0$, step $\alpha=0.05$.

**Full gradient:** $\nabla l(w) = \frac{1}{2}\sum_i (wx_i-y_i)x_i = \frac{1}{2}[(-2)(1)+(-4)(2)+(-5)(3)+(-8)(4)] = \frac{1}{2}(-2-8-15-32) = -28.5$.

Update: $w_1 = 0 + 0.05\times28.5 = 1.425$.

**SGD (1 sample):** sample $i=2$: $\nabla l_2(w_0) = (0\cdot2-4)\cdot2 = -8$. $w_1^{\text{SGD}} = 0+0.05\times8=0.4$. Much noisier per step but each step costs $4\times$ less.

**Key:** noisy gradient = high variance but escapes shallow local minima; decay step size $\alpha_n \to 0$ ensures eventual convergence.

---

**Coordinate Descent** — minimise $f(x,y) = (x-3)^2 + 2(y-1)^2$. True min at $(3,1)$.

Start $(0,0)^T$.

- **Step 1** (fix $y=0$, minimise over $x$): $\partial f/\partial x = 2(x-3)=0 \Rightarrow x=3$. Move to $(3,0)$.
- **Step 2** (fix $x=3$, minimise over $y$): $\partial f/\partial y = 4(y-1)=0 \Rightarrow y=1$. Move to $(3,1)$. ✓

2-step convergence because $f$ is **separable** (no $xy$ cross-term). For $f(x,y)=(x-y)^2+x^2$, the Hessian has off-diagonal terms and coordinate descent would zigzag — many steps needed.

---

## Table 5: Numerical Integration — 1D Methods (Original + Composite)

**Goal:** Approximate $\displaystyle I = \int_a^b f(x)\,dx$.

### Variable Key

| Symbol | Meaning |
|--------|---------|
| $h = b-a$ | Strip width (single rule) |
| $n$ | Number of strips (composite rule) |
| $h = (b-a)/n$ | Strip width (composite rule) |
| $m = (a+b)/2$ | Midpoint of strip |
| $x_i = a + ih$, $i=0,\ldots,n$ | Equally-spaced nodes |
| $f_i = f(x_i)$ | Function value at node $i$ |
| $m_i = x_{i-1} + h/2$ | Midpoint of $i$-th strip |
| $\xi \in (a,b)$ | Unknown point in error term (guaranteed to exist by MVT/Taylor) |

### Single-Strip Rules

| Property | **Midpoint** | **Trapezium** | **Simpson's** | **Boole's** |
|---|---|---|---|---|
| **Formula** | $(b-a)\,f(m)$ | $\dfrac{b-a}{2}[f(a)+f(b)]$ | $\dfrac{b-a}{6}[f(a)+4f(m)+f(b)]$ | $\dfrac{b-a}{90}[7f_0+32f_1+12f_2+32f_3+7f_4]$ |
| **Nodes used** | $m = (a+b)/2$ only | $a$ and $b$ | $a,\;m=(a+b)/2,\;b$ | 5 equally-spaced: $a,\;a+h/4,\;a+h/2,\;a+3h/4,\;b$ |
| **$f$ evaluations** | 1 | 2 | 3 | 5 |
| **Error formula** | $-\dfrac{(b-a)^3}{24}f''(\xi)$ | $+\dfrac{(b-a)^3}{12}f''(\xi)$ | $-\dfrac{(b-a)^5}{2880}f^{(4)}(\xi)$ | $O\!\left((b-a)^7 f^{(6)}\right)$ |
| **Exact for poly. degree** | $\leq 1$ | $\leq 1$ | $\leq 3$ (one better than expected!) | $\leq 5$ |
| **Evaluates endpoints?** | No | Yes | Yes | Yes |

> **Why Simpson's is exact for degree $\leq 3$:** You might expect degree $\leq 2$ (fitting a quadratic with 3 nodes). But the error formula has $f^{(4)}(\xi)$ not $f^{(3)}(\xi)$: the 3rd-order error term vanishes by symmetry of the nodes around $m$. This is a free bonus from symmetric spacing.

> **Relationship:** Midpoint error $= -\tfrac{h^3}{24}f''$; trapezium error $= +\tfrac{h^3}{12}f'' = -2\times\text{midpoint error}$. Simpson's $=\tfrac{2\,\text{midpoint}+\text{trapezium}}{3}$, so both leading errors cancel — achieving $O(h^5)$ at no extra cost!

### Composite Rules ($n$ strips, strip width $h=(b-a)/n$)

| Property | **Composite Midpoint** | **Composite Trapezium** | **Composite Simpson's** ($n$ even) | **Monte Carlo** (nD) |
|---|---|---|---|---|
| **Formula** | $h\displaystyle\sum_{i=1}^{n}f(m_i)$ | $h\!\left[\tfrac{1}{2}f(a)+\displaystyle\sum_{i=1}^{n-1}f(x_i)+\tfrac{1}{2}f(b)\right]$ | $\dfrac{h}{3}\!\left[f_0 + 4(f_1+f_3+\cdots+f_{n-1}) + 2(f_2+\cdots+f_{n-2}) + f_n\right]$ | $A(R)\cdot\dfrac{1}{N}\displaystyle\sum_{i=1}^N f(X_i)$, $X_i\stackrel{\text{iid}}{\sim}\text{Uniform}(R)$ |
| **Error bound** | $\dfrac{(b-a)^3}{24n^2}\max \vert f'' \vert $ | $\dfrac{(b-a)^3}{12n^2}\max \vert f'' \vert $ | $\dfrac{(b-a)^5}{180n^4}\max \vert f^{(4)} \vert $ | $O(N^{-1/2})$ **independent of dimension $d$** |
| **Convergence rate** | $O(n^{-2})$ | $O(n^{-2})$ | $O(n^{-4})$ | $O(N^{-1/2})$ |
| **$f$ evaluations** | $n$ (midpoints only) | $n+1$ (all boundary + interior) | $2n+1$ (all $x_0,\ldots,x_n$) | $N$ |
| **Smoothness required** | $f\in C^2[a,b]$ | $f\in C^2[a,b]$ | $f\in C^4[a,b]$; **$n$ must be even** | $f$ integrable; no smoothness needed |
| **Evaluates endpoints?** | No (safe near singularities) | Yes | Yes | N/A |
| **Limitations** | Error bound requires bounded $f''$ | Same; trapezium $\approx 2\times$ larger error than midpoint | $n$ must be even; need bounded $f^{(4)}$; negative weights in high-order Newton–Cotes | Slow $O(N^{-1/2})$; probabilistic; only competitive for $d\geq 5$ |
| **When to use** | Endpoint singularity; $f$ only $C^2$; no endpoint evaluation | Same smoothness as midpoint; standard baseline | **Default 1D workhorse** for $C^4$ smooth $f$ | $d\geq 5$; irregular region; $f$ is a black-box simulation |

**Effect of doubling $n$:** Midpoint/trapezium error $\to$ error$/4$. Simpson's error $\to$ error$/16$. Simpson's is vastly more efficient for smooth $f$.

**When to match rule to smoothness:**

| Smoothness of $f$ | Best rule | Reason |
|---|---|---|
| $C^2$ only (kink, limited derivatives) | Composite midpoint or trapezium | Higher rules need $f^{(4)}$, which doesn't exist — convergence degrades |
| $C^4$ | Composite Simpson's | Best effort-to-accuracy tradeoff |
| $C^6$ | Boole's | Gets $O(n^{-6})$ |
| Very smooth, $d=1$ | Gaussian quadrature (not in syllabus) | Optimal convergence |
| $d\geq 5$, any smoothness | Monte Carlo | Curse of dimensionality: tensor grids cost $O(n^d)$ |
| Endpoint singularity (e.g. $f(a)=\infty$) | Composite midpoint | Avoids evaluating at dangerous endpoints |

---

### Worked Examples

**Midpoint (single strip)** — $\int_0^1 x^2\,dx = 1/3$.

$m=0.5$. Approx $= (1)(0.25) = 0.25$. Error $= 0.333-0.25 = 0.083$.

Check formula: $-\frac{1^3}{24}f''(\xi) = -\frac{1}{24}\cdot2 = -\frac{1}{12}\approx -0.083$ ✓ (formula gives the signed correction to add).

---

**Composite Midpoint ($n=4$)** — $\int_0^1 x^2\,dx$.

$h=0.25$. Midpoints: $0.125, 0.375, 0.625, 0.875$. $f$ values: $0.015625, 0.140625, 0.390625, 0.765625$.

$\approx 0.25\times(0.015625+0.140625+0.390625+0.765625) = 0.25\times1.3125 = 0.328125$

Error $= 0.333-0.328 = 0.005$. Error bound: $\frac{1}{24\times16}\times2 = \frac{1}{192}\approx0.0052$ ✓

---

**Trapezium (single strip)** — $\int_0^1 x^2\,dx$.

$\approx \frac{1}{2}[f(0)+f(1)] = \frac{1}{2}[0+1] = 0.5$. Error $= 0.333-0.5 = -0.167$.

Formula: $+\frac{1}{12}\times2 = +0.167$ ✓ (trapezium **overestimates** for convex $f$; correction is negative).

---

**Composite Trapezium ($n=4$)** — $\int_0^1 x^2\,dx$.

$h=0.25$. Nodes: $0,0.25,0.5,0.75,1$. Values: $0, 0.0625, 0.25, 0.5625, 1$.

$\approx 0.25\times[\tfrac{1}{2}(0)+0.0625+0.25+0.5625+\tfrac{1}{2}(1)] = 0.25\times[0+0.0625+0.25+0.5625+0.5] = 0.25\times1.375=0.34375$

Error $= 0.333-0.344=-0.010$. Error bound: $\frac{1}{12\times16}\times2=\frac{1}{96}\approx0.0104$ ✓

---

**Simpson's (single strip)** — $\int_0^1 x^2\,dx$.

$\approx \frac{1}{6}[f(0)+4f(0.5)+f(1)] = \frac{1}{6}[0+1+1] = \frac{1}{3}$. **Exact!** ($x^2$ has degree $\leq 3$.)

---

**Composite Simpson's ($n=4$)** — $\int_0^\pi \sin(x)\,dx = 2$.

$h=\pi/4$. Values: $f_0=\sin(0)=0$, $f_1=\sin(\pi/4)=\frac{\sqrt{2}}{2}\approx0.7071$, $f_2=\sin(\pi/2)=1$, $f_3=\sin(3\pi/4)\approx0.7071$, $f_4=\sin(\pi)=0$.

Pattern $1$–$4$–$2$–$4$–$1$:

$$\approx\frac{\pi/4}{3}[0+4(0.7071)+2(1)+4(0.7071)+0] = \frac{\pi}{12}[2.8284+2+2.8284] = \frac{\pi}{12}\times7.6569\approx2.0046$$

Error $\approx0.0046$. Compare: composite trapezium with $n=4$ gives error $\approx0.082$. Simpson's is $\sim18\times$ more accurate with the same nodes!

---

**Boole's Rule (single strip)** — $\int_0^1 e^x\,dx = e-1\approx1.71828$.

5 equally-spaced nodes: $f_0=1,\;f_1=e^{0.25}\approx1.2840,\;f_2=e^{0.5}\approx1.6487,\;f_3=e^{0.75}\approx2.1170,\;f_4=e\approx2.7183$.

$$\approx\frac{1}{90}[7(1)+32(1.2840)+12(1.6487)+32(2.1170)+7(2.7183)] = \frac{1}{90}[7+41.09+19.78+67.74+19.03] = \frac{154.64}{90}\approx1.71822$$

Error $\approx0.00006$. Highly accurate with just 5 evaluations.

---

## Table 6: Error, Uncertainty, and Accuracy

### Part A — Error Terminology

| Concept | Symbol/Formula | What it measures | Notes |
|---------|---------------|-----------------|-------|
| **Error** | $\tilde{u}-u$ | Signed difference: computed $-$ true | Can be positive or negative |
| **Absolute error** | $ \vert \tilde{u}-u \vert $ | Magnitude of error | For vectors: use $ \Vert \tilde{\mathbf{u}}-\mathbf{u} \Vert $ (Euclidean norm) |
| **Relative error** | $ \vert \tilde{u}-u \vert / \vert u \vert $ | Error scaled by true value | Dimensionless; meaningful only if $u\neq0$ |
| **Forward error** | $\tilde{f}(x)-f(x)$ | How wrong is the **output**? | Direct measure of answer quality |
| **Backward error** | $\min \vert \delta x \vert $ s.t. $\tilde{f}(x)=f(x+\delta x)$ | How wrong is the **input** that would give this output? | Small backward error $\Rightarrow$ computed answer is exact for a nearby input |
| **Truncation error** | e.g. $O(h^p)$ in integration | Error from approximating an $\infty$-procedure by finite one | Inherent to the mathematical method; exists even with exact arithmetic |
| **Roundoff error** | $\leq \varepsilon_{\text{mach}}\cdot \vert u \vert $ per operation | Error from finite floating-point representation | Inherent to hardware; not the same as truncation error |

**Machine epsilon $\varepsilon_{\text{mach}}$:** Smallest $\varepsilon>0$ such that $1+\varepsilon \neq 1$ in floating point. Equivalently: maximum relative rounding error when storing a real number.

| Precision | Bits | $\varepsilon_{\text{mach}}$ |
|-----------|------|--------------------------|
| Single (float32) | 23 mantissa bits | $\approx 6\times10^{-8}$ |
| Double (float64) | 52 mantissa bits | $\approx 1.1\times10^{-16}$ |

---

### Part B — Catastrophic Cancellation

Subtracting two nearly equal floating-point numbers cancels the leading significant digits, leaving only the (inaccurate) trailing digits.

**Example:** Compute $1-\cos(10^{-4})$ in double precision.

$\cos(10^{-4}) \approx 1 - 5\times10^{-9}$ (differs from 1 in the 9th decimal place). Subtraction: the 15 leading digits of $\cos(10^{-4})$ cancel with those of 1, leaving only the last few digits — relative error can be $\sim10^7$ times machine epsilon.

**Remedy:** Algebraic reformulation. Use $1-\cos x = 2\sin^2(x/2)$:

$$2\sin^2(5\times10^{-5})\approx2(5\times10^{-5})^2 = 5\times10^{-9}$$

This is computed directly with full double-precision accuracy. **Rule:** always rewrite to avoid subtracting nearly-equal quantities.

---

### Part C — Rates of Convergence

Let $\epsilon_n = x_n - x^*$ be the error at step $n$ of an iterative method. As $n\to\infty$:

| Type | Formal condition | Intuition | Digits gained per step | Example |
|------|-----------------|-----------|----------------------|---------|
| **Sublinear** | $ \vert \epsilon_{n+1} \vert / \vert \epsilon_n \vert \to 1$ | Convergence slows down | $<1$ per step | Naive fixed-point near flat region |
| **Logarithmic** | Sublinear + step ratios $\to1$ | Barely converging; extremely slow | $\approx0$ | Bisection in a pathological case |
| **Linear** (order 1) | $ \vert \epsilon_{n+1} \vert / \vert \epsilon_n \vert \to r$, $0<r<1$ | Error multiplied by constant $r<1$ each step | $-\log_{10}(r)$ per step (constant) | Bisection ($r=0.5$, $\approx 0.30$ digits/step); gradient descent |
| **Superlinear** | $ \vert \epsilon_{n+1} \vert / \vert \epsilon_n \vert \to0$ | Faster than any fixed linear rate | Grows with $n$ | Secant ($q\approx1.618$); BFGS |
| **Order $q>1$** | $ \vert \epsilon_{n+1} \vert / \vert \epsilon_n \vert ^q\to C>0$ | Error raised to power $q$ each step | $\sim q^n\times$(initial digits) | General |
| **Quadratic** ($q=2$) | $ \vert \epsilon_{n+1} \vert / \vert \epsilon_n \vert ^2\to C$ | Correct digits **double** each step | Doubles | Newton's (root-finding + optimisation) |
| **Cubic** ($q=3$) | $ \vert \epsilon_{n+1} \vert / \vert \epsilon_n \vert ^3\to C$ | Correct digits **triple** each step | Triples | Halley's method |

**Lemma 4.2 (required):** If $A_n$ converges with order $q$, then $A_{2n}$ (every other term) converges with order $q^2$.

*Practical use:* Compare method A (order $q$, costs 1 unit/step) vs method B (order $p$, costs 2 units/step). Method B is worth it only if $p > q^2$.

*Example:* Newton (order 2, costs 2 evals: $f+f'$) vs secant (order 1.618, costs 1 eval per step). For fair comparison, 2 secant steps ($q^2=1.618^2\approx2.618$) cost 2 evals and give order 2.618 $>$ 2. So **secant beats Newton per function evaluation** when $f'$ is expensive!

---

### Part D — Fundamental Tolerance Limits

Near a minimum, $f(\mathbf{x})-f(\mathbf{x}^*)\approx\tfrac{1}{2}\|\mathbf{x}-\mathbf{x}^*\|^2_H$ (quadratic behaviour). A perturbation $\delta$ in $\mathbf{x}$ causes $\sim\delta^2$ change in $f$. When $\delta^2\sim\varepsilon_{\text{mach}}$, the values $f(\mathbf{x})$ and $f(\mathbf{x}^*)$ are **indistinguishable in floating point**:

| Setting | Achievable tolerance $ \Vert \mathbf{x}-\mathbf{x}^* \Vert / \Vert \mathbf{x}^* \Vert $ | Reason |
|---------|----------------------------------------------------------------------|--------|
| **Root-finding** | $O(\varepsilon_{\text{mach}})\approx10^{-16}$ | $f$ is approx. linear near root: perturbation $\delta$ in $x$ gives $\sim f'(x^*)\delta$ in $f$ — both are representable |
| **Optimisation** | $O(\sqrt{\varepsilon_{\text{mach}}})\approx10^{-8}$ | $f$ is approx. quadratic near min: perturbation $\delta$ in $x$ gives $\sim\delta^2$ in $f$ — distinguishable only when $\delta\gtrsim\sqrt{\varepsilon}$ |

**Sensible defaults:** tolerance $10^{-8}$ for optimisation; $10^{-12}$ or smaller for root-finding.

---

### Part E — Termination Criteria

**Root-finding** — stop when any of:
1. **Step small:** $|x_n - x_{n-1}| < \text{tol}(1+|x_n|)$ — relative change in $x$ below tol
2. **Residual small:** $|f(x_n)| < \text{tol}$ — backward error small
3. **Budget exhausted:** $n = N_{\max}$
4. **Step ill-defined:** $f'(x_n)\approx 0$, singular Jacobian, etc.

**Optimisation** — stop when any of:
1. **Step small:** $\|\mathbf{x}_n-\mathbf{x}_{n-1}\|<\text{tol}(1+\|\mathbf{x}_n\|)$
2. **Value change small:** $|f(\mathbf{x}_n)-f(\mathbf{x}_{n-1})|<\text{tol}(1+|f(\mathbf{x}_n)|)$ *(caution: $f$ could creep to $-\infty$ with tiny steps)*
3. **Gradient small:** $\|\nabla f(\mathbf{x}_n)\|<\text{tol}(1+\|\nabla f(\mathbf{x}_0)\|)$ — **best default when gradient available**
4. **Budget exhausted** or step ill-defined

---

### Part F — A Posteriori Error Estimates

**For root-finding** — given candidate $\tilde{x}$:

$$|\tilde{x}-x^*|\approx\frac{|f(\tilde{x})|}{|f'(\tilde{x})|}$$

This is the size of the next Newton step. It is exact when $f$ is linear; reliable when $\tilde{x}$ is already close to $x^*$.

**Example:** $f(x)=x^2-2$, $\tilde{x}=1.5$. $|f(1.5)|=0.25$, $|f'(1.5)|=3$. Estimate: $0.25/3=0.083$. True error: $|1.5-\sqrt{2}|=0.086$. Good estimate!

---

### Part G — Quick Comparison: Optimisation vs Root-Finding

| | **Root-Finding** | **Optimisation** |
|---|---|---|
| Solving | $\mathbf{f}(\mathbf{x})=\mathbf{0}$ | $\nabla f(\mathbf{x})=\mathbf{0}$ |
| "Better" iterate | Smaller residual $ \Vert \mathbf{f}(\mathbf{x}_n) \Vert $ (not always monotone) | Smaller value $f(\mathbf{x}_n)$ — can always enforce with Armijo |
| Best achievable tol | $O(\varepsilon_{\text{mach}})\approx10^{-16}$ | $O(\sqrt{\varepsilon_{\text{mach}}})\approx10^{-8}$ |
| Newton variant | $x_{n+1}=x_n-f/f'$; quadratic near simple root; fails if $J$ singular | $\mathbf{d}=-H^{-1}\mathbf{g}$; quadratic near min with PD Hessian; fails if $H$ not PD |
| Second-order structure | **Jacobian** $J$ (generally non-symmetric) | **Hessian** $H$ (always symmetric) |
| Key extra requirement | Residual shrinks | Must ensure **descent** ($\mathbf{g}^T\mathbf{d}<0$), not just stationarity |
| Multiple solutions | Can have many roots — method finds one near $\mathbf{x}_0$ | Local minima vs global minima; global min generally NP-hard |

---

### Part H — Error Analysis for Composite Integration (Required Proofs)

**How to derive composite error bounds** from single-strip bounds:

Single-strip midpoint error on $[x_{i-1},x_i]$: $E_i = -\frac{h^3}{24}f''(\xi_i)$ for some $\xi_i\in(x_{i-1},x_i)$.

Sum over $n$ strips: $E = \sum_{i=1}^n E_i = -\frac{h^3}{24}\sum_{i=1}^n f''(\xi_i)$.

Since $f''$ is continuous on $[a,b]$ and $\frac{1}{n}\sum f''(\xi_i)$ is an average of values of $f''$, by the Intermediate Value Theorem there exists $\xi\in(a,b)$ such that $\frac{1}{n}\sum f''(\xi_i) = f''(\xi)$.

So $|E| = \frac{h^3}{24}\cdot n\cdot|f''(\xi)| \leq \frac{h^3}{24}\cdot n\cdot\max|f''|$.

Substitute $h = (b-a)/n$: $|E| \leq \frac{(b-a)^3}{24n^2}\max|f''|$. $\square$

The same argument applies to trapezium and Simpson's (substituting their single-strip error formulas).

---

---

## Table 7: Chain Rule and Product Rule — Jacobian and Gradient Forms

**Context:** When differentiating composed or multiplied functions, the form of the rule depends on the *types* of $f$ and $g$. This table covers all combinations that appear in this course.

### Notation Key

| Symbol | Meaning |
|--------|---------|
| Scalar function $\mathbb{R}\to\mathbb{R}$ | Plain 1D function $f(t)$; derivative is $f'(t)\in\mathbb{R}$ |
| Scalar field $\mathbb{R}^n\to\mathbb{R}$ | Multi-input, scalar output; gradient $\nabla f\in\mathbb{R}^n$, Jacobian $J_f = \nabla f^T \in\mathbb{R}^{1\times n}$ |
| Vector field $\mathbb{R}^n\to\mathbb{R}^m$ | Multi-input, vector output; Jacobian $J_f\in\mathbb{R}^{m\times n}$ with $(J_f)_{ij}=\partial f_i/\partial x_j$ |
| $\nabla f(\mathbf{x})\in\mathbb{R}^n$ | Gradient: column vector of partial derivatives (only defined for scalar fields) |
| $J_f(\mathbf{x})\in\mathbb{R}^{m\times n}$ | Jacobian matrix (defined for any $f:\mathbb{R}^n\to\mathbb{R}^m$) |

> **Relationship:** For a scalar field $f:\mathbb{R}^n\to\mathbb{R}$, the Jacobian is a $1\times n$ row vector, i.e. $J_f = \nabla f^T$. The Hessian $H(f) = J(\nabla f)$ — the Jacobian of the gradient.

---

### Chain Rule

| $f$ type | $g$ type | $(f\circ g)$ type | Jacobian rule $J_{f\circ g}(\mathbf{x})$ | Gradient rule $\nabla(f\circ g)(\mathbf{x})$ |
|---|---|---|---|---|
| Vector field $\mathbb{R}^n\to\mathbb{R}^m$ | Vector field $\mathbb{R}^p\to\mathbb{R}^n$ | Vector field $\mathbb{R}^p\to\mathbb{R}^m$ | $J_f(g(\mathbf{x}))\,J_g(\mathbf{x})$ | N/A (not scalar output) |
| Scalar field $\mathbb{R}^n\to\mathbb{R}$ | Vector field $\mathbb{R}^p\to\mathbb{R}^n$ | Scalar field $\mathbb{R}^p\to\mathbb{R}$ | $J_f(g(\mathbf{x}))\,J_g(\mathbf{x})$ | $J_g(\mathbf{x})^T\,\nabla f(g(\mathbf{x}))$ |
| Scalar function $\mathbb{R}\to\mathbb{R}$ | Scalar field $\mathbb{R}^p\to\mathbb{R}$ | Scalar field $\mathbb{R}^p\to\mathbb{R}$ | $f'(g(\mathbf{x}))\,J_g(\mathbf{x})$ | $f'(g(\mathbf{x}))\,\nabla g(\mathbf{x})$ |
| Scalar function $\mathbb{R}\to\mathbb{R}$ | Scalar function $\mathbb{R}\to\mathbb{R}$ | Scalar function $\mathbb{R}\to\mathbb{R}$ | $f'(g(x))\,g'(x)$ | $f'(g(x))\,g'(x)$ (same) |

**Key observation:** The Jacobian rule is always matrix multiplication in the form $J_f \cdot J_g$ (outer $\times$ inner). The gradient version replaces $J_g^T$ and $\nabla f$ since gradients are column vectors.

---

### Product Rule

| Operation | $f$ type | $g$ type | Result type | Jacobian rule $J_{fg}(\mathbf{x})$ | Gradient rule $\nabla(fg)(\mathbf{x})$ |
|---|---|---|---|---|---|
| Scalar $\times$ Vector | Scalar field $\mathbb{R}^n\to\mathbb{R}$ | Vector field $\mathbb{R}^n\to\mathbb{R}^m$ | Vector field $\mathbb{R}^n\to\mathbb{R}^m$ | $g(\mathbf{x})\,J_f(\mathbf{x}) + f(\mathbf{x})\,J_g(\mathbf{x})$ | N/A (not scalar output) |
| Scalar $\times$ Scalar | Scalar field $\mathbb{R}^n\to\mathbb{R}$ | Scalar field $\mathbb{R}^n\to\mathbb{R}$ | Scalar field $\mathbb{R}^n\to\mathbb{R}$ | $g(\mathbf{x})\,J_f(\mathbf{x}) + f(\mathbf{x})\,J_g(\mathbf{x})$ | $g(\mathbf{x})\,\nabla f(\mathbf{x}) + f(\mathbf{x})\,\nabla g(\mathbf{x})$ |
| Dot product $\mathbf{f}\cdot\mathbf{g}$ | Vector field $\mathbb{R}^n\to\mathbb{R}^m$ | Vector field $\mathbb{R}^n\to\mathbb{R}^m$ | Scalar field $\mathbb{R}^n\to\mathbb{R}$ | $\mathbf{f}(\mathbf{x})^T J_g(\mathbf{x}) + \mathbf{g}(\mathbf{x})^T J_f(\mathbf{x})$ | $J_g(\mathbf{x})^T\mathbf{f}(\mathbf{x}) + J_f(\mathbf{x})^T\mathbf{g}(\mathbf{x})$ |
| Standard 1D | Scalar function $\mathbb{R}\to\mathbb{R}$ | Scalar function $\mathbb{R}\to\mathbb{R}$ | Scalar function $\mathbb{R}\to\mathbb{R}$ | $f'(x)\,g(x) + f(x)\,g'(x)$ | $f'(x)\,g(x) + f(x)\,g'(x)$ (same) |

**Why the dot-product gradient transposes:** $\mathbf{f}\cdot\mathbf{g} = \mathbf{f}^T\mathbf{g}$ is a scalar. Differentiating via the chain rule gives a column-vector gradient, requiring $J^T$ (which is $n\times m$ times $m\times 1$) to produce an $n\times1$ result.

---

### Worked Examples

**Chain rule (scalar field ∘ vector field):** Let $f(\mathbf{u}) = \|\mathbf{u}\|^2 = u_1^2+u_2^2$ and $\mathbf{g}(\mathbf{x}) = (x_1 x_2,\, x_1+x_2)^T$. Find $\nabla(f\circ\mathbf{g})$ at $\mathbf{x}=(1,2)^T$.

$\mathbf{g}(1,2) = (2,3)^T$. $\nabla f(\mathbf{u}) = 2\mathbf{u}$, so $\nabla f(g(\mathbf{x})) = (4,6)^T$.

$J_g(\mathbf{x}) = \begin{pmatrix}x_2 & x_1 \\ 1 & 1\end{pmatrix} = \begin{pmatrix}2&1\\1&1\end{pmatrix}$ at $(1,2)$.

$$\nabla(f\circ g)(\mathbf{x}) = J_g(\mathbf{x})^T\nabla f(g(\mathbf{x})) = \begin{pmatrix}2&1\\1&1\end{pmatrix}^T\begin{pmatrix}4\\6\end{pmatrix} = \begin{pmatrix}2&1\\1&1\end{pmatrix}\begin{pmatrix}4\\6\end{pmatrix} = \begin{pmatrix}14\\10\end{pmatrix}$$

**Check dimensions:** $J_g^T\in\mathbb{R}^{2\times2}$, $\nabla f\in\mathbb{R}^2$ → result $\in\mathbb{R}^2$ ✓ (gradient of a scalar field on $\mathbb{R}^2$).

---

**Chain rule (scalar function ∘ scalar field):** Let $f(t)=e^t$ and $g(\mathbf{x})=x_1^2+x_2^2$. Find $\nabla(f\circ g)$.

$f'(t)=e^t$, $\nabla g(\mathbf{x}) = (2x_1,2x_2)^T$.

$$\nabla(f\circ g)(\mathbf{x}) = f'(g(\mathbf{x}))\,\nabla g(\mathbf{x}) = e^{x_1^2+x_2^2}\begin{pmatrix}2x_1\\2x_2\end{pmatrix}$$

This is the standard "outer derivative times inner gradient" rule.

---

**Product rule (dot product):** Let $\mathbf{f}(\mathbf{x}) = \mathbf{x}$ and $\mathbf{g}(\mathbf{x}) = A\mathbf{x}$ (with $A$ symmetric). Then $\mathbf{f}\cdot\mathbf{g} = \mathbf{x}^T A\mathbf{x}$.

$J_f = I$ (identity), $J_g = A$.

$$\nabla(\mathbf{f}\cdot\mathbf{g}) = J_g^T\mathbf{f} + J_f^T\mathbf{g} = A^T\mathbf{x} + I^T A\mathbf{x} = A\mathbf{x} + A\mathbf{x} = 2A\mathbf{x}$$

(using $A^T=A$ for symmetric $A$). This recovers the standard result $\nabla(\mathbf{x}^T A\mathbf{x})=2A\mathbf{x}$.

---

**Product rule (scalar × scalar):** Let $f(\mathbf{x}) = x_1^2$ and $g(\mathbf{x}) = x_2^3$. Find $\nabla(fg)$.

$\nabla f = (2x_1,0)^T$, $\nabla g = (0,3x_2^2)^T$.

$$\nabla(fg) = g\,\nabla f + f\,\nabla g = x_2^3\begin{pmatrix}2x_1\\0\end{pmatrix} + x_1^2\begin{pmatrix}0\\3x_2^2\end{pmatrix} = \begin{pmatrix}2x_1 x_2^3\\3x_1^2 x_2^2\end{pmatrix}$$

Verify directly: $fg = x_1^2 x_2^3$, $\partial/\partial x_1 = 2x_1 x_2^3$ ✓, $\partial/\partial x_2 = 3x_1^2 x_2^2$ ✓.

---

*End of Summary Tables*
