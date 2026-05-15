# Domino from Cyclic Matrix

A Python implementation of an algorithm that computes invariants of the Dieudonné module of a supersingular abelian $g$-fold with cyclic Frobenius action.

---

## Mathematical Background

### Setting

Let $X$ be a supersingular abelian $g$-fold over an algebraically closed field $k$ of characteristic $p$. Its covariant Dieudonné module $M$ is a free $W(k)$-module of rank $2g$, equipped with a Frobenius-semilinear operator $F$.

This program works under the assumption that $F$ acts **cyclically**: there exists a basis $e_1, \ldots, e_{2g}$ of $M$ such that

$$F(e_i) = p^{a_i} u_i \, e_{i+1 \bmod 2g}$$

for units $u_i \in W(k)$. The **exponent sequence** $a = (a_1, \ldots, a_{2g})$ records the $p$-adic valuations of the $2g$ entries of this cyclic matrix. Because $X$ is supersingular of dimension $g$, the exponent sequence satisfies:

- $a_i \in \{0, 1\}$ for all $i$
- $\sum a_i = g$

Two exponent sequences differing by a cyclic rotation define isomorphic Dieudonné modules, so we work with **rotation equivalence classes**.

### Invariants computed

The algorithm computes invariants of $M' := \operatorname{Dom}(E_1^{02})$, where $E_1^{02} = H^2(X, W\mathcal{O}_X)$ is a term in the slope spectral sequence, and $\operatorname{Dom}$ denotes the domino part of the differential $d : H^2(X, W\mathcal{O}_X) \to H^2(X, W\Omega^1_X)$.

For each exponent sequence (and each of its $\Lambda^2$-components), the program outputs:

| Invariant | Description |
|-----------|-------------|
| **a-number** | Number of cyclic $0 \to 1$ transitions in the exponent sequence; equals $\dim_k \operatorname{Hom}(\alpha_p, X[p])$ |
| **dim** | Dimension of the relevant piece of $M'$ |
| **p-exp** | The $p$-exponent: the maximum height of the corresponding lattice path |
| **type-seq** | A nondecreasing sequence of positive integers encoding the recursive level structure of the lattice path arcs; equivalently, the node sizes of a rooted tree |
| **Isog** | A nonincreasing partition derived from the type-seq tree via leaf-peeling and conjugation; encodes the isogeny structure |
| **σ** | Generalized Artin invariant: the sum of all type-seq node sizes across all $\Lambda^2$-components |

---

## Algorithm

### Step 1 — Exponent sequences

Given $g \geq 1$, enumerate all binary sequences of length $2g$ with exactly $g$ ones, **up to cyclic rotation**. The canonical representative of each class is chosen as the lexicographically smallest rotation.

### Step 2 — $\Lambda^2$ decomposition

The exterior square of the cyclic Dieudonné module decomposes into $g$ components. For a sequence $a = (a_1, \ldots, a_{2g})$ and shift $i \in \{1, \ldots, g\}$:

- For $1 \leq i \leq g-1$: component of length $2g$ with $j$-th entry $a_j + a_{j+i \bmod 2g}$
- For $i = g$: component of length $g$ with $j$-th entry $a_j + a_{j+g}$

Each component has entries in $\{0, 1, 2\}$ and sum equal to its length, reflecting that $\Lambda^2$ of a supersingular module has all slopes $1/2$ after this decomposition.

### Step 3 — Lattice path (b-sequence)

Each component is converted to a **closed lattice path** $b = (b_1, \ldots, b_n)$ by:

$$b_1 = 0, \qquad b_{k+1} = b_k + a_k - 1$$

Since $\sum a_i = n$, the path closes: $b_n = b_1$ before normalization. We then shift so $\min(b_i) = 0$. Each step satisfies $b_{k+1} - b_k = a_k - 1 \in \{-1, 0, +1\}$.

### Step 4 — Indecomposable decomposition

Because the path is closed, $b$ is treated as a **cyclic sequence**. The path is decomposed into maximal contiguous arcs of nonzero values (zeros are the "ground level"). Each such arc is called **indecomposable**; arcs may wrap around the end of the sequence.

### Step 5 — dim, p-exp, and type-seq

For each indecomposable arc $(b_1, \ldots, b_m)$:

- **dim** $= \#\{i : b_{i+1} = b_i - 1\}$ (where $b_{m+1} := 0$): the number of descent steps, counting the final return to zero
- **p-exp** $= \max_i b_i$: the height of the arc
- **type-seq**: a rooted tree built recursively:

  $$\operatorname{type\text{-}seq}(b_1, \ldots, b_m) = \operatorname{Tree}\!\bigl(m,\; \text{children from type-seq of indecomp. parts of } (b_1{-}1,\ldots,b_m{-}1)\bigr)$$

  The recursion terminates when all values reach zero (producing a leaf). The integer labels stored at each node are the arc lengths.

For a full component, dim and p-exp are summed/maximized over all indecomposable arcs, and the type-seq trees form a forest.

### Step 6 — Isog partition

The type-seq forest is processed by **iterative leaf-peeling**:

1. Count all current leaf nodes; append count to a list
2. Remove all leaves; repeat until the forest is empty

This yields a sequence of counts. Sorted in decreasing order, this is a partition of the total number of nodes. Its **conjugate partition** (transpose of the Young diagram) is $\operatorname{Isog}(M')$.

*Example*: partition $[4,2,1]$ conjugates to $[3,2,1,1]$.

The **generalized Artin invariant** $\sigma$ is the sum of all type-seq node sizes across all $\Lambda^2$-components of a given exponent sequence.

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
| `step6(trees)` | Isog via leaf-peeling + conjugation |
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
