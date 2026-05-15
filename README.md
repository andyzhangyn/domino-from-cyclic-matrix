# Domino from Cyclic Matrix

A Python implementation of an algorithm that computes invariants of the formal Brauer group of a supersingular abelian $g$-fold with cyclic Frobenius action.

---

## Mathematical Background

### Setting

Let $X$ be a supersingular abelian $g$-fold over an algebraically closed field $k$ of characteristic $p$. The covariant Dieudonné module $M = H^1_{\mathrm{crys}}(X/W)$ is the first crystalline cohomology of $X$: a free $W(k)$-module of rank $2g$ carrying the structure of an $F$-crystal.

This program works under the assumption that $F$ acts **cyclically** on $M$: there exists a basis $e_1, \ldots, e_{2g}$ such that

$$F(e_i) = p^{a_i} u_i \, e_{i+1 \bmod 2g}$$

for units $u_i \in W(k)$. The **exponent sequence** $a = (a_1, \ldots, a_{2g})$ records the $p$-adic valuations of the entries of this cyclic matrix. Because $X$ is supersingular of dimension $g$:

- $a_i \in \{0, 1\}$ for all $i$
- $\sum a_i = g$

Two exponent sequences differing by a cyclic rotation define isomorphic $F$-crystals, so we work with **rotation equivalence classes**.

### From $H^1$ to $H^2$: the exterior square

The second crystalline cohomology is computed as the exterior square of the $F$-crystal:

$$H^2_{\mathrm{crys}}(X/W) = \Lambda^2 H^1_{\mathrm{crys}}(X/W) = \Lambda^2 M.$$

### Nygaard reconstruction and the formal Brauer group

For a Mazur–Ogus variety $X/k$ (i.e. torsion-free crystalline cohomology and Hodge–de Rham degeneration), the **Nygaard reconstruction theorem** takes as input the $F$-crystal $H^2_{\mathrm{crys}}(X/W)$ and recovers the **formal Brauer group** $\widehat{\mathrm{Br}}_X$.

The formal group $\widehat{\mathrm{Br}}_X$ is associated to the Dieudonné module $H^2(X, W\mathcal{O}_X)$ (the slope-0 part of $H^2_{\mathrm{crys}}$), which is precisely the Cartier–Dieudonné module of $p$-typical curves on $\widehat{\mathrm{Br}}_X$.

For **supersingular** abelian varieties this reconstruction is especially clean: $H^2(X, W\mathcal{O}_X)$ has no finite-free Dieudonné part and consists entirely of **domino** pieces. Consequently $\widehat{\mathrm{Br}}_X$ is **unipotent**: a finite iterated extension of $\widehat{\mathbb{G}}_a$.

### The domino $U$

Denote by $U$ the domino part of $H^2(X, W\mathcal{O}_X)$, i.e.

$$U := \operatorname{Dom}\!\bigl(H^2(X, W\mathcal{O}_X) \xrightarrow{d} H^2(X, W\Omega^1_X)\bigr).$$

This algorithm computes the following invariants of $U$ (and hence of $\widehat{\mathrm{Br}}_X$):

| Invariant | Description |
|-----------|-------------|
| **a-number** | Number of cyclic $0 \to 1$ transitions in the exponent sequence; equals $\dim_k \operatorname{Hom}(\alpha_p, X[p])$ |
| **dim** | $\dim U$: depends only on $g$ and the Newton polygon of $H^2$ (all Newton polygons coincide in the supersingular case) |
| **p-exp** | Smallest power of $p$ annihilating $U$ |
| **type-seq** | By Ekedahl, every domino carries a unique decreasing filtration whose $i$-th associated graded piece is a direct sum of elementary dominoes $U_i$. The type sequence is the nondecreasing list of integers $i$ appearing in this graded, counted with multiplicity |
| **Isog** | Isogeny type of $\widehat{\mathrm{Br}}_X$, computed via the Greene–Kleitman algorithm (leaf-peeling + partition conjugation) |
| **σ** | **Generalized Artin invariant**: sum of all type-seq entries across all $\Lambda^2$-components. Generalizes the classical Artin invariant from the $g = 2$ case, where $\dim U = 1$ and the type sequence consists of a single integer $\in \{1, 2\}$ |

---

## Algorithm

### Step 1 — Exponent sequences

Given $g \geq 1$, enumerate all binary sequences of length $2g$ with exactly $g$ ones, **up to cyclic rotation**. The canonical representative of each class is chosen as the lexicographically smallest rotation.

### Step 2 — $\Lambda^2$ decomposition

The exterior square $\Lambda^2 M$ decomposes into $g$ components indexed by a shift $i \in \{1, \ldots, g\}$. For $a = (a_1, \ldots, a_{2g})$:

- For $1 \leq i \leq g-1$: component of length $2g$ with $j$-th entry $a_j + a_{j+i \bmod 2g}$
- For $i = g$: component of length $g$ with $j$-th entry $a_j + a_{j+g}$

Each component has entries in $\{0, 1, 2\}$ and sum equal to its length, reflecting that $\Lambda^2 M$ has all slopes $1/2$.

### Step 3 — Lattice path (b-sequence)

Each component $a$ is converted to a **closed lattice path** $b = (b_1, \ldots, b_n)$ by:

$$b_1 = 0, \qquad b_{k+1} = b_k + a_k - 1.$$

Since $\sum a_i = n$, the path closes ($b_n = b_1$ before normalization). We shift so $\min(b_i) = 0$. Each step satisfies $b_{k+1} - b_k = a_k - 1 \in \{-1, 0, +1\}$.

### Step 4 — Indecomposable decomposition

Because the path is closed, $b$ is treated as a **cyclic sequence**. It is decomposed into maximal contiguous arcs of nonzero values (zeros are the "ground level"); arcs may wrap around the endpoint. Each arc is called **indecomposable**.

### Step 5 — dim, p-exp, and type-seq

For each indecomposable arc $(b_1, \ldots, b_m)$:

- **dim** $= \#\{i : b_{i+1} = b_i - 1\}$ (with $b_{m+1} := 0$)
- **p-exp** $= \max_i b_i$
- **type-seq**: a rooted tree built recursively,

  $$\operatorname{type\text{-}seq}(b_1, \ldots, b_m) = \operatorname{Tree}\!\bigl(m,\; \text{children from type-seq of indecomp. parts of } (b_1{-}1,\ldots,b_m{-}1)\bigr)$$

  terminating when all values reach zero (leaf node). Node labels are arc lengths; the multiset of all labels, sorted nondecreasing, is the type sequence of $U$.

For a full $\Lambda^2$-component, dim and p-exp are summed/maximized over all arcs.

### Step 6 — Isog partition

The type-seq forest is processed by the **Greene–Kleitman algorithm** (iterative leaf-peeling):

1. Count all current leaf nodes; append count to a list.
2. Remove all leaves; repeat until the forest is empty.

The resulting counts, sorted in decreasing order, form a partition of the total node count. Its **conjugate partition** (transpose of the Young diagram) is $\operatorname{Isog}(\widehat{\mathrm{Br}}_X)$.

*Example*: $[4, 2, 1] \mapsto [3, 2, 1, 1]$.

The **generalized Artin invariant** $\sigma$ is the sum of all type-seq entries across all $\Lambda^2$-components.

---

## Implementation

```
domino.py        — main algorithm (Steps 1–6)
visualize.py     — plots the distribution of σ for g = 2..12
output_*.log     — precomputed results for various g
artin_invariant.png / sum_frequency.png — generated figures
```

### `domino.py`

| Function | Role |
|----------|------|
| `step1(g)` | Canonical exponent sequences up to rotation |
| `step2(a)` | $\Lambda^2$ decomposition |
| `step3(a)` | Exponent sequence → closed lattice path |
| `step4(b)` | Cyclic indecomposable decomposition |
| `step5(b)` | dim, p-exp, type-seq trees |
| `step6(trees)` | Isog via Greene–Kleitman + conjugation |
| `a_number(a)` | Number of cyclic $0 \to 1$ transitions |
| `analyze(g)` | Verbose output for a single $g$ |
| `analyze_compact(g)` | One line per exponent sequence |

### Usage

```bash
# Compact output (default)
python domino.py 3

# Verbose output with full per-component detail
python domino.py -v 3

# Multiple values
python domino.py 1 2 3 4 5

# Redirect to file
python domino.py 1 2 3 4 > output.log
```

### Output format

**Compact** (one line per exponent sequence):
```
g = 3
   1.  a-seq=[0, 0, 0, 1, 1, 1]  a-number=1  type-seq=[2, 4]  sum=6  Isog=[2]
   2.  a-seq=[0, 0, 1, 0, 1, 1]  a-number=2  type-seq=[1, 2, 4]  sum=7  Isog=[1, 1]
   ...
```

**Verbose** (`-v`): shows the b-sequence, indecomposable arcs, and per-component invariants for each $\Lambda^2$-component.

---

## Visualization

`visualize.py` reads the precomputed log files and plots the distribution of $\sigma$ for $g = 2, \ldots, 12$, with bars color-coded by a-number. Vertical lines mark the superspecial value $\sigma = g(g-1)/2$ and the supergeneral value $\sigma = g^2(g-1)/2$.

```bash
python visualize.py   # produces sum_frequency.png
```

![Artin invariant distribution](artin_invariant.png)
