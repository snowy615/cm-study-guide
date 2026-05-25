# DAA Study Guide — Oxford HT 2026
*Sandra Kiefer · Based on all 8 lecture PDFs*

---

## Table of Contents
1. [Asymptotics & Program Cost](#1-asymptotics--program-cost)
2. [Divide and Conquer](#2-divide-and-conquer)
3. [Data Structures: Heaps & Priority Queues](#3-data-structures-heaps--priority-queues)
4. [Dynamic Programming](#4-dynamic-programming)
5. [Graph Decomposition](#5-graph-decomposition)
6. [Paths in Graphs](#6-paths-in-graphs)
7. [Greedy Algorithms](#7-greedy-algorithms)
8. [Stable Matching](#8-stable-matching)

---

## 1. Asymptotics & Program Cost

### 1.1 Why Efficiency Matters

A 10⁶× speed increase only adds 20 to the solvable problem size when T(n) = 2ⁿ, versus a 10⁶× increase for T(n) = n.  
**Lesson:** algorithmic efficiency dominates hardware improvements.

### 1.2 Insertion Sort

#### Pseudocode
```
Insertion-Sort(A)
  for j = 1 to A.length - 1
    key = A[j+1]
    // Insert A[j+1] into sorted sequence A[1..j]
    i = j
    while i > 0 and A[i] > key
      A[i+1] = A[i]
      i = i - 1
    A[i+1] = key
```

#### Correctness Proof (Loop Invariant Method)

**Outer loop invariant (I):** At the start of the j-th iteration, `A[1..j]` consists of the elements originally in `A[1..j]` in sorted order.

| Phase | Argument |
|-------|----------|
| **Initialisation** | When j=1, `A[1..1]` is a singleton — trivially sorted. |
| **Maintenance** | The inner while loop inserts `key = A[j+1]` into its correct place in `A[1..j]`, yielding a sorted `A[1..j+1]`. |
| **Termination** | Loop exits when j = n. Then `A[1..n]` is sorted — algorithm is correct. |

**Inner while loop invariant:**
- **(I1)** `A[1..i] · A[i+2..j+1]` is `a₁,...,aⱼ` in sorted order  
- **(I2)** All elements in `A[i+2..j+1]` are strictly greater than `key`

| Phase | Argument |
|-------|----------|
| **Initialisation** | i=j: I1 holds by assumption; I2 holds because `A[j+2..j+1]` is empty. |
| **Maintenance** | If `A[i] > key`, shift `A[i]` right; I1 and I2 still hold with i decremented. |
| **Termination** | Exits when i=0 or `A[i] <= key`. Both: inserting `key` at `A[i+1]` gives sorted `A[1..j+1]`. |

#### Complexity
| Case | Cost | Why |
|------|------|-----|
| Worst | Θ(n²) | Reverse-sorted input; tⱼ = j+1 |
| Best | Θ(n) | Already sorted; tⱼ = 1 |
| Average | Θ(n²) | tⱼ ≈ j/2 (same asymptotic class) |

**Properties:** stable, in-place, online.

---

### 1.3 Asymptotic Notation

#### Big-O (upper bound)
> **Definition.** f ∈ O(g) if ∃ n₀ ∈ ℕ⁺, c ∈ ℝ⁺ such that ∀ n ≥ n₀, f(n) ≤ c·g(n).

#### Big-Ω (lower bound)
> f(n) = Ω(g(n)) iff g(n) = O(f(n)).

#### Big-Θ (tight bound)
> f = Θ(g) iff ∃ c₁, c₂, n₀ such that ∀ n ≥ n₀: c₁·g(n) ≤ f(n) ≤ c₂·g(n).

#### Key Lemma (Properties of Big-O)
For f, g, h : ℕ → ℝ⁺:
1. Constant factors can be absorbed into O.
2. If f₁ ∈ O(g₁) and f₂ ∈ O(g₂), then f₁+f₂ ∈ O(max(g₁,g₂)).
3. If f₁ ∈ O(g₁) and f₂ ∈ O(g₂), then f₁·f₂ ∈ O(g₁·g₂).
4. O is transitive.
5. Every polynomial of degree ℓ ≥ 0 is in O(nˡ).
6. lgᶜ(n) ∈ O(nᵈ) for any constants c, d > 0.
7. nᶜ ∈ O(dⁿ) for any constant c > 0, d > 1.
8. For all a, b > 1: log_a(n) = Θ(log_b(n)) — so we write Θ(log n) without base.

---

## 2. Divide and Conquer

### 2.1 Merge Sort

#### Pseudocode
```
Merge-Sort(A, p, r)
  if r > p + 1
    q = floor((p + r)/2)
    Merge-Sort(A, p, q)
    Merge-Sort(A, q, r)
    Merge(A, p, q, r)
// Initial call: Merge-Sort(A, 1, n+1)

Merge(A, p, q, r)
  n1 = q - p;  n2 = r - q
  Create arrays L[1..n1+1] and R[1..n2+1]
  for i = 1 to n1:  L[i] = A[p+i-1]
  for j = 1 to n2:  R[j] = A[q+j-1]
  L[n1+1] = INF;  R[n2+1] = INF
  i = 1;  j = 1
  for k = p to r-1
    if L[i] <= R[j]
      A[k] = L[i];  i = i+1
    else
      A[k] = R[j];  j = j+1
```

#### Complexity
- **Merge:** Θ(n) where n = r-p.
- **Merge-Sort recurrence:** T(n) = 2T(n/2) + Θ(n), giving **T(n) = Θ(n log n)**.

**Properties:** stable (left-biased merge), NOT in-place (Θ(n) extra), NOT online.

#### Recurrence Analysis — Three Methods

**1. Guess-and-test (substitution)**
Guess T(n) ≤ an lg n. Base: T(2) = 2c+2d ≤ 2a if a ≥ c+d. Inductive step shows T(n) ≤ an lg n if a ≥ d. ✓

**2. Recursion tree**
Level ℓ has cost Cₗ(n) ∈ [d′n, dn]. With lg n levels + cn leaf cost:
d′n lg n + cn ≤ T(n) ≤ dn lg n + cn → **T(n) = Θ(n log n)**.

**3. Master Theorem**

> **Theorem (DPV/CLRS).** Suppose T(n) ≤ aT(⌈n/b⌉) + O(nᵈ), a > 0, b > 1, d ≥ 0. Then:
> - T(n) = O(nᵈ)          if d > log_b a
> - T(n) = O(nᵈ log n)    if d = log_b a
> - T(n) = O(n^{log_b a}) if d < log_b a

**Proof sketch:** Recursion tree with log_b(n)+1 levels. Level k has aᵏ nodes each of size n/bᵏ. Cost at level ℓ ≤ c·nᵈ·(a/bᵈ)ˡ. Sum the geometric series — three cases depending on ratio a/bᵈ vs 1.

**Merge-Sort:** a=b=2, d=1 → d = log₂2 → case 2 → T(n)=O(n log n). ✓

---

### 2.2 Binary Search

```
BinSearch(A, p, r, z)
  // A sorted in increasing order
  if p >= r:  return "No"
  q = floor((p + r)/2)
  if z = A[q]:  return "Yes"
  else if z < A[q]:  return BinSearch(A, p, q, z)
  else:  return BinSearch(A, q+1, r, z)
```
**T(n) = T(ceil(n/2)) + O(1) → T(n) = O(log n)** (Master Theorem: a=1, b=2, d=0 → case 2).

---

### 2.3 Selection (Median of Medians) — O(n) worst case

```
Select(A, i)
  1. Divide n elements into ceil(n/5) groups of <= 5.
  2. Find median of each group (by Insertion-Sort on each group) -> O(n).
  3. x = Select(medians, ceil(n/10))   // median of medians
  4. Partition A around x -> lower L (< x), upper U (> x).
  5. k = |L|
  6a. if i = k+1: return x
  6b. if i < k+1: return Select(L, i)
  6c. if i > k+1: return Select(U, i-k-1)
```

**Key bound:** |L|, |U| <= floor(7n/10 + 6).

**Proof of size bound:** At least ceil(ceil(n/5)/2) - 2 groups have medians >= x; each such group contributes >= 3 elements > x. So #elements > x >= 3(n/10 - 2) = 3n/10 - 6, giving |L| <= n - (3n/10-6) = 7n/10 + 6. Symmetric for upper.

**Recurrence:** T(n) ≤ T(ceil(n/5)) + T(floor(7n/10+6)) + bn.

**Proof T(n) = O(n) — guess T(n) ≤ cn:**
Substitute: T(n) ≤ c(n/5+1) + c(7n/10+6) + bn = 9cn/10 + 7c + bn = cn + (-cn/10 + 7c + bn).
This is ≤ cn if c ≥ 10bn/(n-70). For n ≥ 140: n/(n-70) ≤ 2, so c ≥ 20b suffices.
Base: c = max{T(n)/n for n ≤ 140, 20b}. Induction closes. ✓

---

### 2.4 Integer Multiplication (Karatsuba)

Naive 4-multiplication split: T(n) = 4T(n/2) + O(n) → O(n²).

**Karatsuba (Gauss trick):** Use only xLyL, xRyR, (xL+xR)(yL+yR) to recover xLyR+xRyL for free.
```
T(n) = 3T(n/2) + O(n)  →  T(n) = O(n^{log₂ 3}) ≈ O(n^{1.585})
```

---

### 2.5 Strassen's Matrix Multiplication

Standard O(n³). Strassen uses 7 sub-multiplications P1–P7:
```
P1 = A(F-H),       P5 = (A+D)(E+H)
P2 = (A+B)H,       P6 = (B-D)(G+H)
P3 = (C+D)E,       P7 = (A-C)(E+F)
P4 = D(G-E)

XY top-left     = P5 + P4 - P2 + P6
XY top-right    = P1 + P2
XY bottom-left  = P3 + P4
XY bottom-right = P1 + P5 - P3 - P7
```
**T(n) = 7T(n/2) + O(n²) → T(n) = O(n^{log₂7}) ≈ O(n^{2.81})**.

---

### 2.6 Sorting Lower Bound & Counting Sort

#### Lower bound for comparison-based sorting
> **Theorem.** Every comparison-based sorting algorithm takes Ω(n log n) in the worst case.

**Proof.** Every algorithm corresponds to a decision tree with n! leaves (one per permutation). A binary tree of depth d has <= 2ᵈ leaves. So depth >= log(n!) = Ω(n log n) (Stirling). ✓

#### Counting Sort (non-comparison — Θ(n+k))
```
CountingSort(A, k)
  // keys in {0,...,k}
  for i = 0 to k:  C[i] = 0
  for j = 1 to n:  C[A[j]]++            // frequency
  for i = 1 to k:  C[i] += C[i-1]       // cumulative
  for j = n downto 1:                    // stable: right-to-left
    B[C[A[j]]] = A[j]
    C[A[j]]--
```
**Time:** Θ(n + k). Stable. Not comparison-based.

---

### 2.7 Fast Fourier Transform (FFT) ★

Multiplies polynomials in O(n log n) arithmetic operations.

**Key idea:** A(x) = B(x²) + xC(x²).  Evaluate at pairs ±v simultaneously.  
Use n-th roots of unity ωₙ = e^{2πi/n} recursively.

```
Recursive-FFT(a[0..n-1])
  if n = 1: return a
  u = Recursive-FFT(a[0, 2, ..., n-2])   // even-index coefficients
  v = Recursive-FFT(a[1, 3, ..., n-1])   // odd-index coefficients
  w = e^{2*pi*i/n};  x = 1
  for k = 0 to n/2-1:
    y[k]       = u[k] + x*v[k]
    y[k + n/2] = u[k] - x*v[k]
    x = x * w
  return y
```
**T(n) = 2T(n/2) + O(n) → T(n) = O(n log n)**.

**Interpolation:** FFT computes y = V(ω)·a. Since V(ω)⁻¹ = V(ω⁻¹)/n, IFFT = same procedure with ω⁻¹, divide by n.

> **Theorem.** Product of two degree-n polynomials can be computed in O(n log n) arithmetic operations.

---

## 3. Data Structures: Heaps & Priority Queues

### 3.1 Binary Heap

A **max-heap** is an essentially complete binary tree satisfying:
> **Max-Heap Property:** A[i] <= A[floor(i/2)] for all 1 < i <= heap-size.

Array layout: root = A[1], left(i) = A[2i], right(i) = A[2i+1], parent(i) = A[floor(i/2)].  
**Height:** floor(log n).

---

### 3.2 Max-Heapify

**Pre:** left and right subtrees of i are max-heaps.  
**Post:** subtree rooted at i is a max-heap.

```
Max-Heapify(A, i)
  l = 2i;  r = 2i+1;  n = A.heap-size
  if l <= n and A[l] > A[i]:  largest = l
  else:  largest = i
  if r <= n and A[r] > A[largest]:  largest = r
  if largest != i
    swap A[i] and A[largest]
    Max-Heapify(A, largest)
```

**Time:** T(n) <= T(2n/3) + Θ(1) → **O(log n)** (Master, case 2 or height argument).

---

### 3.3 Make-Max-Heap

```
Make-Max-Heap(A)
  A.heap-size = n
  for i = ceil((n+1)/2) - 1 downto 1
    Max-Heapify(A, i)
```

**Time: O(n)** — not O(n log n).

**Proof:** T(n) <= Σ_{h=0}^{floor(log n)} (n/2^h)*(c*h) = cn * Σ_{h=0}^∞ h/2^h = cn * 2 = O(n).
(Uses Σ kxᵏ = x/(1-x)² at x=1/2.)

---

### 3.4 Heapsort

```
Heapsort(A)
  Make-Max-Heap(A)
  for i = A.heap-size downto 2
    swap A[1] and A[i]
    A.heap-size -= 1
    Max-Heapify(A, 1)
```

**Loop invariant:** A[i+1..n] is sorted; elements in A[1..i] are <= all in A[i+1..n].

**Time:** O(n) + O(n)·O(log n) = **O(n log n)**.  
**Properties:** in-place, NOT stable, worst-case O(n log n).

---

### 3.5 Priority Queue Operations

| Operation | Time |
|-----------|------|
| `Heap-Maximum(A)` → return A[1] | Θ(1) |
| `Heap-Extract-Max(A)` | O(log n) |
| `Heap-Increase-Key(A, i, key)` | O(log n) |
| `Heap-Insert(A, key)` | O(log n) |

```
Heap-Extract-Max(A)
  if A.heap-size < 1: error
  max = A[1]
  A[1] = A[A.heap-size];  A.heap-size -= 1
  Max-Heapify(A, 1)
  return max

Heap-Increase-Key(A, i, key)
  if key < A[i]: error
  A[i] = key
  while i > 1 and A[Parent(i)] < A[i]
    swap A[i] and A[Parent(i)];  i = Parent(i)

Heap-Insert(A, key)
  A.heap-size += 1
  A[A.heap-size] = -INF
  Heap-Increase-Key(A, A.heap-size, key)
```

---

## 4. Dynamic Programming

### 4.1 The DP Paradigm

1. **Define subproblems** (ordered smallest → largest).
2. **Optimal substructure:** optimal solution built from optimal sub-solutions.
3. **Solve bottom-up**, storing results.

---

### 4.2 Change-Making Problem

**Subproblem:** C[u] = min coins to make value u.  
**Recurrence:** C[u] = 1 + min{ C[u - xᵢ] : 1 <= i <= n, u >= xᵢ }

```
Change-Making(x₁,...,xₙ; v)
  C[0] = 0
  for u = 1 to v
    C[u] = 1 + min{ C[u - xᵢ] : 1 <= i <= n, u >= xᵢ }
  return C[v]
```
**Time:** O(nv).

---

### 4.3 Knapsack

#### With repetition (unbounded)
**Recurrence:** K[w] = max{ K[w - wᵢ] + vᵢ : i ∈ {1..n}, wᵢ <= w }
```
K[0] = 0
for w = 1 to W
  K[w] = max{ K[w - wᵢ] + vᵢ : wᵢ <= w }
return K[W]
```
**Time:** O(nW).

#### Without repetition (0-1 knapsack)
**Subproblem:** K[w, j] = max value, capacity w, items 1..j.  
**Recurrence:** K[w, j] = max{ K[w-wⱼ, j-1] + vⱼ,  K[w, j-1] }
```
for j = 0 to n:  K[0,j] = 0
for w = 0 to W:  K[w,0] = 0
for j = 1 to n
  for w = 1 to W
    if wⱼ > w:  K[w,j] = K[w, j-1]
    else:        K[w,j] = max{K[w-wⱼ, j-1]+vⱼ,  K[w, j-1]}
return K[W, n]
```
**Time:** O(nW).

---

### 4.4 Longest Increasing Subsequence (LIS)

**Subproblem:** L[j] = length of longest IS ending at j.  
**Recurrence:** L[j] = 1 + max{ L[i] : i < j, A[i] < A[j] }

```
Longest-Increasing-Subsequence(A)
  L[1]=1;  P[1]=NIL;  k=1
  for j = 2 to n
    L[j]=1;  P[j]=NIL
    for i = 1 to j-1
      if A[i] < A[j] and L[i] >= L[j]
        L[j] = 1+L[i];  P[j] = i
    if L[j] > L[k]:  k = j
  // Reconstruct via backpointers P
  Create B[1..L[k]]
  for j = L[k] downto 1
    B[j] = A[k];  k = P[k]
```
**Time:** O(n²).

#### Faster LIS — O(n log r) with Binary Search

Maintain K[1..r]: K[i] = index of smallest tail of IS of length i (A[K[1]] < ... < A[K[r]]).

```
Longest-Increasing-Subsequence-BS(A)
  r=1;  K[1]=1;  P[1]=NIL;  K[0]=NIL
  for j = 2 to n
    // Binary search for smallest i s.t. A[j] <= A[K[i]]
    i0=1;  i1=r+1
    while i0 < i1
      im = floor((i0+i1)/2)
      if A[j] <= A[K[im]]:  i1 = im
      else:  i0 = im+1
    if i0 > r:  r = r+1
    K[i0] = j;  P[j] = K[i0-1]
  // Reconstruct
  k = K[r]
  for j = r downto 1
    B[j] = A[k];  k = P[k]
```
**Time:** O(n log r), where r = LIS length.

---

### 4.5 Edit Distance (Levenshtein)

**Subproblem:** E[i,j] = edit distance between x[1..i] and y[1..j].  
**Recurrence:** E[i,j] = min{ E[i-1,j]+1,  E[i,j-1]+1,  E[i-1,j-1]+δ(i,j) }  
where δ(i,j) = 1 if x[i] ≠ y[j], else 0.  
**Base cases:** E[0,j]=j, E[i,0]=i.  
**Time:** O(mn).

---

### 4.6 TSP via DP

**Subproblem:** C[S,j] = shortest path starting at 1, ending at j, through all nodes of S.  
**Recurrence:** C[S,j] = min{ C[S\{j},i] + d_{ij} : i ∈ S\{1,j} }  
**Base:** C[{1},1] = 0.  
**Answer:** min{ C[{1..n},j] + dⱼ₁ : j ≠ 1 }  
**Time:** O(n² 2ⁿ).

> **Theorem.** TSP is NP-hard. DP does NOT produce a polynomial-time solution.

**Failure of Optimal Substructure example (Longest Simple Path):**  
l[a,c] ≠ l[a,b] + l[b,c] in general because paths may share vertices.  
**Lesson:** DP requires optimal substructure; not all problems have it.

---

## 5. Graph Decomposition

### 5.1 Graph Representations

| Representation | Space | Edge lookup |
|----------------|-------|-------------|
| Adjacency matrix | O(|V|²) | O(1) |
| Adjacency lists | O(|V|+|E|) | O(degree) |

---

### 5.2 Depth-First Search (DFS)

```
DFS(V, E)
  for u in V:  colour[u] = white;  pi[u] = NIL
  time = 0
  for u in V:
    if colour[u] = white:  DFS-Visit(u)

DFS-Visit(u)
  time++; d[u] = time; colour[u] = grey
  for v in Adj[u]
    if colour[v] = white
      pi[v] = u;  DFS-Visit(v)
  time++; f[u] = time; colour[u] = black
```
**Time:** Θ(|V| + |E|).

#### Parenthesis Theorem
> **Theorem.** For any u,v exactly one holds:
> 1. [d[u],f[u]] and [d[v],f[v]] are disjoint (neither is ancestor of the other).
> 2. [d[v],f[v]] ⊂ [d[u],f[u]] → v is a descendant of u.
> 3. [d[u],f[u]] ⊂ [d[v],f[v]] → u is a descendant of v.

**Corollary:** v is a descendant of u iff d[u] < d[v] < f[v] < f[u].

#### White Path Theorem
> v is a descendant of u iff at time d[u], ∃ all-white path u → v.

**Proof:** "Only if": immediate (DFS visits all white vertices reachable from grey u). "If": induction on path length. Length 1: v ∈ Adj[u], will be discovered. Length ℓ+1: by IH, uₗ is a descendant of u → d[uₗ] < f[u]. When uₗ discovered, v is either white (→ descendant of uₗ, hence u) or already stamped (Parenthesis Theorem → v still in interval of u). ✓

#### Edge Classification
| Type | Colour of v when (u,v) explored | Time stamps |
|------|----------------------------------|-------------|
| Tree | white | d[u] < d[v] < f[v] < f[u] |
| Back | grey | d[v] < d[u] < f[u] < f[v] |
| Forward | black (desc) | d[u] < d[v] < f[v] < f[u] |
| Cross | black (non-desc) | d[v] < f[v] < d[u] < f[u] |

#### Cycle Detection
> **Lemma.** A directed graph has a cycle iff DFS has a back edge.

**Proof:**
- ⇐: Back edge (u,v) → v is an ancestor of u → cycle via DFS path v→u.
- ⇒: Let ⟨v₀,...,vₖ,v₀⟩ be cycle. Let vᵢ = first discovered node. All others reachable from vᵢ → descendants → edge (vᵢ₋₁,vᵢ) is a back edge. ✓

---

### 5.3 Topological Sort

```
Topological-Sort(V, E)
  // G must be a DAG
  Run DFS(V, E), computing f[v] for all v
  Output V in decreasing f-order
```
**Time:** Θ(|V| + |E|).

#### Correctness
> **Proposition.** In a DAG: (u,v) ∈ E ⟹ f[u] > f[v].

**Proof.** When (u,v) explored, u is grey.
- v **white**: v becomes descendant of u → Parenthesis Thm: f[v] < f[u].
- v **black**: f[v] already stamped < current time < f[u].
- v **grey**: impossible in a DAG (would be back edge → cycle). ✓

---

### 5.4 Strongly Connected Components (Kosaraju)

```
SCC(G)
  1. DFS(G)  →  compute f[u] for all u
  2. Compute G^T (transpose)
  3. DFS(G^T), processing vertices in decreasing f-order (from step 1)
  4. Each DFS tree in step 3 = one SCC
```
**Time:** Θ(|V| + |E|).

#### Key Lemmas
1. **SCC structure:** If edge from C to C', no edge back. ⟹ SCC-DAG is a DAG.  
   **Proof:** If path C→C' and C'→C existed, they'd be one SCC. Contradiction.

2. **Finishing time lemma:** Edge from SCC C to SCC C' ⟹ f[C] > f[C'].  
   **Proof:** Same argument as Topological Sort correctness.

3. **Correctness:** Vertex u with max f[u] is in the SCC C with max f[C]. In G^T, no outgoing edges from C → DFS from u visits exactly C. Repeat for remaining SCCs in decreasing f-order. ✓

---

## 6. Paths in Graphs

### 6.1 Breadth-First Search (BFS)

**Problem:** single-source shortest paths, unit weights.

```
BFS(V, E, s)
  d[s]=0;  pi[s]=nil
  for u in V - {s}:  d[u]=INF;  pi[u]=nil
  Q = {};  Enqueue(Q, s)
  while Q != {}
    u = Dequeue(Q)
    for v in Adj[u]
      if d[v] = INF
        d[v] = d[u]+1;  pi[v] = u;  Enqueue(Q, v)
```
**Time:** O(|V| + |E|).

#### Correctness (three-lemma proof)

**Lemma 1:** d[v] < ∞ ⟹ ∃ path of length d[v] from s to v.  
*Proof by induction:* Base d[v]=0: v=s trivially. Step: pi[v]=u, d[v]=d[u]+1; IH gives path to u; add edge (u,v).

**Corollary 1:** d[v] >= δ(s,v).

**Lemma 2 (queue):** Values in queue Q differ by at most 1; earlier-enqueued vertices have ≤ d-values.  
*Proof by induction on operations.*

**Lemma 3 (upper bound):** ∃ path of length ℓ from s to v ⟹ d[v] ≤ ℓ.  
*Proof by induction on ℓ.* Three cases for when v is enqueued relative to vₗ₀.

> **Theorem 1.** d[v] = δ(s,v) for all v.  
> (Follows from Corollary 1 + Lemma 3)

> **Theorem 2.** pi[v] is the predecessor of v on a shortest path (when d[v] > 0).

---

### 6.2 Dijkstra's Algorithm

**Problem:** single-source shortest paths, non-negative weights.

```
Dijkstra(V, E, w, s)
  for v in V:  d[v]=INF;  pi[v]=nil
  d[s] = 0
  Q = Make-Min-Queue(V)   // keyed by d[v]
  while Q != {}
    u = Extract-Min(Q)
    for v in Adj[u]
      if d[u] + w(u,v) < d[v]          // relaxation step
        d[v] = d[u] + w(u,v)
        pi[v] = u
        Decrease-Key(Q, v, d[v])
```
**Time (min-heap):** O((|V|+|E|) log |V|).  
*Improved to O(|E| + |V| log |V|) with Fibonacci heap.*

#### Loop Invariants
- **I1:** d[v] >= δ(s,v) ∀v (always overestimate).
- **I2:** ∀v ∈ S = V\Q: d[v] = δ(s,v) (finished vertices: exact).

**Init:** I1 clear; I2 vacuous (S=∅).  
**Termination:** S=V, I2 → d[v]=δ(s,v) ∀v. ✓

**Maintenance of I1:** If d[v] changes: d[v] = d[u]+w(u,v) >= δ(s,u)+w(u,v) >= δ(s,v). ✓ (triangle inequality)

**Maintenance of I2 — Convergence property:**
> **Lemma.** If s→...→x→y is a shortest path and d[x]=δ(s,x) before relaxing (x,y), then d[y]=δ(s,y) after.  
> **Proof.** After: d[y] ≤ d[x]+w(x,y) = δ(s,x)+w(x,y) = δ(s,y). Also d[y] >= δ(s,y) by I1. Done. ✓

Full I2 maintenance: let y = first vertex outside S on shortest path to v, x = its predecessor in S. Convergence property → d[y]=δ(s,y). Since d[v] ≤ d[y] = δ(s,y) ≤ δ(s,v) and d[v] >= δ(s,v): d[v]=δ(s,v). ✓

---

### 6.3 Bellman-Ford Algorithm

**Problem:** single-source shortest paths, arbitrary weights; detects negative cycles.

```
Bellman-Ford(V, E, w, s)
  for v in V:  d[v]=INF;  pi[v]=nil
  d[s] = 0
  for i = 1 to |V|-1
    for each edge (u,v) in E:
      if d[u] + w(u,v) < d[v]
        d[v] = d[u]+w(u,v);  pi[v]=u
  // Negative-cycle check:
  for each edge (u,v) in E:
    if d[u]+w(u,v) < d[v]:  return false
  return true
```
**Time:** Θ(|V|·|E|).

#### Correctness (no negative cycles)

**Invariant:**
- **I1:** d[v] >= δ(s,v) ∀v (maintained by each relaxation step, same proof as Dijkstra).
- **I2:** Before iteration i, for every shortest path ⟨v₀,...,vₖ⟩: d[vⱼ]=δ(s,vⱼ) for j < i.

**Init:** I2 holds for j=0 (d[s]=0). **Termination:** i=|V|, I2 gives d[v]=δ(s,v) for paths with ≤ |V|-1 edges. ✓

**Maintenance of I2:** d-values non-increasing across iterations. Processing edge (vᵢ₋₁,vᵢ): d[vᵢ] ≤ d[vᵢ₋₁]+w = δ(s,vᵢ₋₁)+w = δ(s,vᵢ). ✓

#### Negative cycle detection
> **Theorem.** If negative cycle c = ⟨v₁,...,vₙ,v₁⟩ is reachable, algorithm returns false.

**Proof.** w(c) < 0 ⟹ Σᵢ[d[vᵢ]+w(vᵢ,vᵢ₊₁)] = Σd[vᵢ] + w(c) < Σd[vᵢ] = Σd[vᵢ₊₁]. So ∃vᵢ: d[vᵢ]+w(vᵢ,vᵢ₊₁) < d[vᵢ₊₁] → check returns false. ✓

---

### 6.4 Shortest Paths in DAGs

```
DAG-Shortest-Path(V, E, w, s)
  Topological-Sort(V, E)
  for v in V:  d[v]=INF;  pi[v]=nil
  d[s] = 0
  for u in V (topological order)
    for v in Adj[u]
      if d[u]+w(u,v) < d[v]
        d[v] = d[u]+w(u,v);  pi[v]=u
```
**Time:** Θ(|V| + |E|). Works with negative weights (no negative cycles possible in DAG).

**Correctness:** Invariant: before processing vₖ, d[vᵢ]=δ(s,vᵢ) for i ≤ k. Each vertex's predecessor in the shortest path comes earlier in topological order → already has exact distance → convergence property applies. ✓

---

### 6.5 Floyd-Warshall (All-Pairs Shortest Paths)

**Subproblem:** d[i,j;k] = shortest i→j path using only intermediate nodes in {1,...,k}.  
**Base:** d[i,j;0] = w(i,j) if edge exists, else ∞.  
**Recurrence:** d[i,j;k+1] = min( d[i,j;k],  d[i,k+1;k] + d[k+1,j;k] )

```
Floyd-Warshall(V, E, w)
  for i,j:  initialise d[i,j;0]
  for k = 0 to |V|-1
    for i = 1 to |V|
      for j = 1 to |V|
        d[i,j;k+1] = min(d[i,j;k], d[i,k+1;k] + d[k+1,j;k])
```
**Time:** O(|V|³). Requires no negative-weight cycles.

---

### 6.6 Shortest Path Algorithm Summary

| Algorithm | Weights | Complexity | Notes |
|-----------|---------|------------|-------|
| BFS | Unit (=1) | O(\|V\|+\|E\|) | — |
| DAG-SP | Arbitrary | Θ(\|V\|+\|E\|) | DAG only |
| Dijkstra | Non-negative | O((\|V\|+\|E\|)log\|V\|) | Greedy |
| Bellman-Ford | Arbitrary | Θ(\|V\|\|E\|) | Detects neg. cycles |
| Floyd-Warshall | No neg. cycles | O(\|V\|³) | All-pairs |

---

## 7. Greedy Algorithms

### 7.1 The Greedy Approach
At each step, make the locally optimal choice — never reconsider.  
Correct for some problems; fails for others.

---

### 7.2 Minimum Spanning Tree — Key Definitions

- **Cut** (S, V\S): partition of vertices.
- **Light edge** crossing a cut: minimum weight among all crossing edges.
- **Safe edge:** (u,v) is safe for A ⊆ MST if A∪{(u,v)} ⊆ some MST.

#### Cut Lemma (core of MST correctness)
> **Lemma.** If A ⊆ some MST, (S,V\S) is a cut **respecting** A (no edge of A crosses it), and (u,v) is a **light edge** crossing it, then (u,v) is **safe** for A.

**Proof.** Let T be an MST with A ⊆ T. T has a unique u-v path p; p must cross the cut — let (x,y) be one such crossing edge.  
Construct T' = T - {(x,y)} ∪ {(u,v)}:
- **T' is a tree:** adding (u,v) creates exactly one cycle (with p); removing (x,y) from that cycle leaves connected + |V|-1 edges.
- **w(T') ≤ w(T):** w(T') = w(T) - w(x,y) + w(u,v) ≤ w(T) since (u,v) is light.
- **A ⊆ T':** (x,y) was removed but (x,y) ∉ A (A respects the cut).
So T' is MST containing A∪{(u,v)}: (u,v) is safe. ✓

#### Generic MST
```
Generic-MST(V, E, w)
  A = {}
  while A is not a spanning tree
    find safe edge (u,v) for A
    A = A ∪ {(u,v)}
  return A
```
**Loop invariant:** A ⊆ some MST. Init: A=∅ trivially ⊆ any MST. Maintenance: only safe edges added. Termination: A is spanning tree + subset of MST → A is an MST. ✓

---

### 7.3 Kruskal's Algorithm

```
Kruskal(V, E, w)
  A = {}
  for v in V:  Make-Set(v)
  Sort E by weight (increasing)
  for each (u,v) in sorted order
    if Find-Set(u) != Find-Set(v)
      A = A ∪ {(u,v)};  Union(u,v)
  return A
```
**Correctness:** Added edge always connects two distinct components → it's a light edge crossing the cut (component(u), V\component(u)), which respects A. By Cut Lemma, safe. ✓  
**Time:** O(|E| log |E|).

**Disjoint-set implementations:**
| Implementation | Time for Kruskal |
|----------------|-----------------|
| Linked list | O(|E| + |V|²) |
| Weighted linked list | O(|E| log |E|) |
| Disjoint-set forest | O(|E| log |E|) effectively |

---

### 7.4 Prim's Algorithm

```
Prim(V, E, w, r)
  for u in V:  key[u]=INF;  pi[u]=nil;  Insert(Q, u)
  Decrease-Key(Q, r, 0)
  while Q != {}
    u = Extract-Min(Q)
    for v in Adj[u]
      if v in Q and w(u,v) < key[v]
        pi[v] = u;  Decrease-Key(Q, v, w(u,v))
```
**Correctness:** Extracted vertex u has minimum key — so edge (pi[u], u) is the light edge crossing cut (V\Q, Q), respecting A. By Cut Lemma, safe. ✓  
**Time:** O(|E| log |V|) (min-heap).

*Note: Prim and Dijkstra have essentially the same structure; key difference is key[v] = min edge weight (Prim) vs d[u]+w(u,v) (Dijkstra).*

---

### 7.5 Activity Selection

**Problem:** max-size set of non-overlapping activities [sᵢ, fᵢ).

#### Greedy Choice Lemma
> **Lemma.** ∃ optimal solution containing the activity with minimum finish time.

**Proof (exchange argument):** Let a = activity with min finish time. Any optimal solution OPT contains some earliest-finishing activity a'. Replace a' with a: f(a) ≤ f(a') so compatibility is preserved. New solution has same size. ✓

```
Activity-Selection(s, f)
  Sort by finish time
  A = {1};  k = 1
  for j = 2 to n
    if s[j] >= f[k]     // compatible with last chosen
      A = A ∪ {j};  k = j
  return A
```
**Time:** O(n log n).  
**Correctness:** Apply greedy choice lemma iteratively: at each step the greedy choice is safe → by induction the full greedy solution is optimal.

---

### 7.6 Tree Facts

> **Lemma.** Undirected G = tree iff each pair of vertices connected by unique path.

> **Lemma.** G = tree ⟹ |E| = |V| - 1.

**Proof:**
- **G connected ⟹ |E| ≥ |V|-1:** n vertices, 0 edges = n components. Each edge reduces #components by ≤ 1. Need ≥ n-1 edges to connect.
- **G acyclic ⟹ |E| ≤ |V|-1:** Induction on n. If >1 component: IH for each. If 1 component (connected tree): remove one edge → 2 components; apply IH; |E| = |E₁|+|E₂|+1 ≤ (|V₁|-1)+(|V₂|-1)+1 = |V|-1. ✓

> **Lemma.** Every connected graph has a spanning tree.  
> **Proof:** Greedily add edges without cycles until no more can be added. Resulting T is maximal acyclic; disconnectedness would contradict maximality (G is connected). ✓

---

## 8. Stable Matching

### 8.1 Problem Definition

**Input:** n hospitals H, n students S. Each hospital ranks all students; each student ranks all hospitals.

**Perfect matching M:** every hospital and student in exactly one pair.

**Unstable pair (h,s) ∉ M:** h prefers s over its match **AND** s prefers h over its match.

**Stable matching:** perfect matching with no unstable pairs.

---

### 8.2 Gale-Shapley Algorithm

```
Gale-Shapley(preference lists)
  M = {}
  while (some hospital h is unmatched and hasn't proposed to everyone)
    s = first student on h's list not yet proposed to
    if s is unmatched
      M = M ∪ {h-s}
    else if s prefers h to current partner h'
      M = (M - {h'-s}) ∪ {h-s}
    // else: s rejects h (nothing changes)
  return M
```

**Observations:**
- Students only ever "trade up" (move to a preferred hospital).
- Hospitals propose in decreasing preference order.

---

### 8.3 Proofs of Correctness

#### Termination
> **Lemma.** At most n² iterations.

**Proof.** Each hospital proposes to each student at most once → ≤ n² proposals. ✓

#### Perfect Matching
> **Lemma.** Output M is a perfect matching.

**Proof (contradiction).** Suppose h is unmatched at termination → h proposed to every student. Once matched, students remain matched. So all n students were proposed to and are matched. But |H| = |S| = n → n matches exist. If h unmatched, some student must be too, but we said all students matched → contradiction. ✓

#### Stability
> **Lemma.** M is stable.

**Proof.** Let (h,s) ∉ M. Two cases:
1. **h never proposed to s:** h proposes in preference order, so h prefers its current match over s. Not an unstable pair.
2. **h proposed to s:** s rejected h (or later replaced h with someone better). Since students only trade up, s's current match is preferred over h. Not an unstable pair. ✓

> **Theorem (Gale-Shapley 1962).** The algorithm always finds a stable matching for any instance.

---

### 8.4 Optimality

- **Hospital-optimal assignment:** each hospital gets its best valid partner (best across all stable matchings). GS returns this.
- **Student-pessimal assignment:** each student gets their worst valid partner. GS also returns this (as a corollary).

**Truthfulness:** Hospitals cannot gain by misreporting. Students **can** sometimes gain by lying.

---

### 8.5 Stable Roommate Counterexample

Stable roommate problem (non-bipartite) need not have a solution:
```
Preferences:  A: B > C > D
              B: C > A > D
              C: A > B > D
              D: any
```
Every possible matching has an unstable pair → **no stable roommate matching exists**.

---

## Quick Reference: Complexity Table

| Algorithm | Time | Space |
|-----------|------|-------|
| Insertion Sort | O(n²) worst; O(n) best | O(1) |
| Merge Sort | Θ(n log n) | Θ(n) |
| Heapsort | O(n log n) | O(1) |
| Counting Sort | Θ(n+k) | Θ(k) |
| Binary Search | O(log n) | O(1) |
| Select (MoM) | O(n) | O(1) |
| Karatsuba | O(n^1.585) | — |
| Strassen | O(n^2.81) | — |
| FFT | O(n log n) | O(n) |
| Make-Max-Heap | O(n) | O(1) |
| DFS | Θ(|V|+|E|) | O(|V|) |
| BFS | O(|V|+|E|) | O(|V|) |
| Topological Sort | Θ(|V|+|E|) | O(|V|) |
| SCC (Kosaraju) | Θ(|V|+|E|) | O(|V|) |
| Dijkstra | O((|V|+|E|)log|V|) | O(|V|) |
| Bellman-Ford | Θ(|V||E|) | O(|V|) |
| DAG Shortest Path | Θ(|V|+|E|) | O(|V|) |
| Floyd-Warshall | O(|V|³) | O(|V|²) |
| Kruskal | O(|E|log|E|) | O(|V|) |
| Prim | O(|E|log|V|) | O(|V|) |
| Gale-Shapley | O(n²) | O(n²) |
| Change-Making (DP) | O(nv) | O(v) |
| Knapsack (DP) | O(nW) | O(nW) |
| LIS (O(n²) DP) | O(n²) | O(n) |
| LIS (binary search) | O(n log r) | O(n) |
| Edit Distance | O(mn) | O(mn) |
| TSP (DP) | O(n² 2ⁿ) | O(n 2ⁿ) |

---

## Key Proof Techniques Summary

| Technique | Used in |
|-----------|---------|
| **Loop invariant** (init, maintenance, termination) | Insertion Sort, Generic MST, Heapsort, BFS, Dijkstra, Bellman-Ford, DAG-SP |
| **Exchange argument** | Cut Lemma (MST), Activity Selection, Gale-Shapley stability |
| **Mathematical induction** | White Path Thm, BFS lemmas, Tree size, Insertion Sort inner loop |
| **Proof by contradiction** | GS perfect matching, Cycle detection lemma |
| **Guess-and-test / substitution** | Merge Sort, Select O(n) proof |
| **Recursion tree** | Merge Sort, Master Theorem proof |
| **Geometric series** | Make-Max-Heap O(n), Master Theorem |
| **Decision tree** | Ω(n log n) sorting lower bound |
| **Triangle inequality** | Dijkstra I1 maintenance |
| **Stirling's approximation** | n! = Ω(n log n), TSP brute force |
