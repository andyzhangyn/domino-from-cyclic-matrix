#!/usr/bin/env python3
"""
Computes invariants of Dieudonné modules of supersingular abelian g-folds.
Steps 1-6 as described in the algorithm notes (Domino from Cyclic Matrix).
"""

import sys
from itertools import combinations
from typing import List, Optional, Tuple


# ─────────────────────────────────────────────────────────────
# Tree for type-seq
# ─────────────────────────────────────────────────────────────

class Tree:
    """Node in the type-seq rooted tree. size = length of the indecomposable interval."""

    def __init__(self, size: int, children: Optional[List["Tree"]] = None):
        self.size = size
        self.children: List[Tree] = children if children is not None else []

    def is_leaf(self) -> bool:
        return not self.children

    def num_nodes(self) -> int:
        return 1 + sum(c.num_nodes() for c in self.children)

    def __repr__(self) -> str:
        if self.is_leaf():
            return f"Leaf({self.size})"
        return f"Node({self.size}, [{', '.join(map(repr, self.children))}])"


# ─────────────────────────────────────────────────────────────
# Step 1: canonical exponent sequences up to cyclic rotation
# ─────────────────────────────────────────────────────────────

def step1(g: int) -> List[Tuple[int, ...]]:
    """
    All binary sequences of length 2g with exactly g ones,
    one representative per cyclic-rotation equivalence class.
    """
    n = 2 * g
    seen: set = set()
    result: List[Tuple[int, ...]] = []

    for ones_pos in combinations(range(n), g):
        ones_set = set(ones_pos)
        seq = tuple(1 if i in ones_set else 0 for i in range(n))
        canon = min(seq[r:] + seq[:r] for r in range(n))
        if canon not in seen:
            seen.add(canon)
            result.append(canon)

    return result


# ─────────────────────────────────────────────────────────────
# Step 2: Λ² decomposition
# ─────────────────────────────────────────────────────────────

def step2(a: Tuple[int, ...]) -> List[Tuple[int, ...]]:
    """
    Λ² decomposition of exponent sequence a (length 2g).
    Returns g components: g-1 of length 2g (i=1..g-1) and 1 of length g (i=g).
    Component i: entry j = a[j] + a[j+i]  (indices mod 2g).
    """
    g = len(a) // 2
    n = len(a)
    components = []

    for i in range(1, g):
        components.append(tuple(a[j % n] + a[(j + i) % n] for j in range(n)))

    # i = g: length-g component
    components.append(tuple(a[j] + a[j + g] for j in range(g)))

    return components


# ─────────────────────────────────────────────────────────────
# Step 3: b-sequence
# ─────────────────────────────────────────────────────────────

def step3(a: Tuple[int, ...]) -> Tuple[int, ...]:
    """
    b[0] = 0,  b[k+1] = b[k] + a[k] - 1.
    Shift so min(b) = 0.
    Since Σa[i] = len(a), the path closes: b[n] = b[0] before shifting.
    Property: b[i+1] - b[i] = a[i] - 1 ∈ {-1, 0, +1}.
    """
    b = [0] * len(a)
    for k in range(1, len(a)):
        b[k] = b[k - 1] + a[k - 1] - 1
    shift = min(b)
    return tuple(x - shift for x in b)


# ─────────────────────────────────────────────────────────────
# Step 4: cyclic indecomposable decomposition
# ─────────────────────────────────────────────────────────────

def step4(b: Tuple[int, ...]) -> List[Tuple[int, ...]]:
    """
    b is treated as a cyclic sequence (it closes: b[n] = b[0] before shift).
    Return maximal non-zero contiguous runs in the cycle; runs may wrap end→start.
    """
    n = len(b)
    if n == 0:
        return []

    # All non-zero: the entire cycle is one interval
    if all(v != 0 for v in b):
        return [b]

    # All zero: no intervals
    if all(v == 0 for v in b):
        return []

    # Find the first zero as the "seam" and walk the cycle from there
    seam = next(i for i, v in enumerate(b) if v == 0)
    intervals: List[Tuple[int, ...]] = []
    i = (seam + 1) % n

    while i != seam:
        if b[i] != 0:
            run: List[int] = []
            while i != seam and b[i] != 0:
                run.append(b[i])
                i = (i + 1) % n
            intervals.append(tuple(run))
        else:
            i = (i + 1) % n

    return intervals


# ─────────────────────────────────────────────────────────────
# Step 5: dim, p-exp, type-seq
# ─────────────────────────────────────────────────────────────

def _dim(interval: Tuple[int, ...]) -> int:
    """
    dim = #{i : iv[i+1] = iv[i] - 1}  with  iv[m] := 0 by convention.
    Counts 'down steps' in the interval, including the final drop to 0.
    """
    ext = list(interval) + [0]
    return sum(1 for i in range(len(interval)) if ext[i + 1] == ext[i] - 1)


def _type_seq(interval: Tuple[int, ...]) -> Tree:
    """
    Recursively build type-seq Tree for an indecomposable interval.
    type_seq([b1,...,bm]) = Tree(m, children)
    where children are type_seq trees of the cyclic-indecomposable parts of (b1-1,...,bm-1).
    """
    b_minus1 = tuple(x - 1 for x in interval)
    sub_intervals = step4(b_minus1)
    return Tree(len(interval), [_type_seq(s) for s in sub_intervals])


def step5(b: Tuple[int, ...]) -> Tuple[int, int, List[Tree]]:
    """
    Returns (dim, p_exp, type_seq_trees) for the b-sequence.
      dim      = total descent count across all indecomposable intervals
      p_exp    = maximum value in b  (0 if b is all zeros)
      trees    = one Tree per indecomposable interval
    """
    intervals = step4(b)
    dim = sum(_dim(iv) for iv in intervals)
    p_exp = max(b) if any(v != 0 for v in b) else 0
    trees = [_type_seq(iv) for iv in intervals]
    return dim, p_exp, trees


# ─────────────────────────────────────────────────────────────
# Step 6: Isog via leaf-peeling
# ─────────────────────────────────────────────────────────────

def _count_leaves(forest: List[Tree]) -> int:
    return sum(1 if t.is_leaf() else _count_leaves(t.children) for t in forest)


def _peel(tree: Tree) -> Optional[Tree]:
    """Remove all leaf nodes. Returns None if tree is itself a leaf."""
    if tree.is_leaf():
        return None
    new_ch = [c for c in (_peel(ch) for ch in tree.children) if c is not None]
    return Tree(tree.size, new_ch)


def _conjugate(partition: List[int]) -> List[int]:
    """Transpose the Young diagram of a partition."""
    if not partition:
        return []
    return [sum(1 for p in partition if p >= j) for j in range(1, max(partition) + 1)]


def step6(trees: List[Tree]) -> Tuple[List[int], List[int]]:
    """
    Leaf-peeling loop → pre-conjugate partition → Isog(M) = conjugate.
    Returns (pre_conj, isog).
    """
    counts: List[int] = []
    forest = list(trees)

    while forest:
        counts.append(_count_leaves(forest))
        forest = [t for t in (_peel(r) for r in forest) if t is not None]

    pre_conj = sorted(counts, reverse=True)
    return pre_conj, _conjugate(pre_conj)


# ─────────────────────────────────────────────────────────────
# Driver
# ─────────────────────────────────────────────────────────────

def _all_sizes(trees: List[Tree]) -> List[int]:
    """Collect all node sizes from a forest, recursively."""
    sizes = []
    for t in trees:
        sizes.append(t.size)
        sizes.extend(_all_sizes(t.children))
    return sizes


def a_number(a: Tuple[int, ...]) -> int:
    """Number of cyclic 0→1 transitions in the exponent sequence."""
    n = len(a)
    return sum(1 for i in range(n) if a[i] == 0 and a[(i + 1) % n] == 1)


def _combined(g: int, a: Tuple[int, ...]) -> Tuple[List[int], List[int]]:
    """Return (combined type-seq sizes sorted, combined Isog sorted desc) for one exponent sequence."""
    combined_sizes: List[int] = []
    combined_isog: List[int] = []
    for comp in step2(a):
        b = step3(comp)
        _, _, trees = step5(b)
        combined_sizes.extend(sorted(_all_sizes(trees)))
        _, isog = step6(trees)
        combined_isog.extend(isog)
    return sorted(combined_sizes), sorted([x for x in combined_isog if x > 0], reverse=True)


def analyze(g: int) -> None:
    print(f"\n{'=' * 62}")
    print(f"  g = {g}   (sequences of length {2 * g}, sum = {g})")
    print(f"{'=' * 62}")

    seqs = step1(g)
    print(f"\nStep 1: {len(seqs)} canonical exponent sequence(s) up to cyclic rotation\n")

    for si, a in enumerate(seqs, 1):
        print(f"{'─' * 55}")
        print(f"Sequence {si}: {list(a)}")

        components = step2(a)
        print(f"Step 2 (Λ²): {len(components)} component(s)")

        combined_sizes: List[int] = []
        combined_isog: List[int] = []

        for ci, comp in enumerate(components, 1):
            label = "length g" if len(comp) == g else "length 2g"
            print(f"\n  ── Component {ci} ({label}) ──")
            print(f"     a-seq    : {list(comp)}")

            b = step3(comp)
            print(f"     b-seq    : {list(b)}")

            indecomp = step4(b)
            print(f"     Indecomp : {[list(iv) for iv in indecomp]}")

            dim, p_exp, trees = step5(b)
            print(f"     dim(M)   = {dim}")
            print(f"     p-exp(M) = {p_exp}")

            sizes = sorted(_all_sizes(trees))
            print(f"     type-seq : {sizes}")

            _, isog = step6(trees)
            print(f"     Isog(M)  = {isog}")

            combined_sizes.extend(sizes)
            combined_isog.extend(isog)

        print(f"\n  ── Combined for sequence {si} ──")
        print(f"     type-seq : {sorted(combined_sizes)}")
        print(f"     Isog     : {sorted([x for x in combined_isog if x > 0], reverse=True)}")

    print()


def analyze_compact(g: int) -> None:
    print(f"\ng = {g}")
    for si, a in enumerate(step1(g), 1):
        ts, isog = _combined(g, a)
        print(f"  {si:2d}.  a-seq={list(a)}  a-number={a_number(a)}  type-seq={ts}  sum={sum(ts)}  Isog={isog}")
    print()


def main() -> None:
    compact = True
    if len(sys.argv) > 1:
        args = sys.argv[1:]
    else:
        try:
            raw = input("Enter g value(s) (integers ≥ 1, space-separated) [-v for verbose]: ")
            args = raw.split()
        except EOFError:
            return

    if "-v" in args:
        compact = False
        args = [a for a in args if a != "-v"]

    for arg in args:
        try:
            g = int(arg)
        except ValueError:
            print(f"Skipping '{arg}': not an integer")
            continue
        if g < 1:
            print(f"Skipping g={g}: must be ≥ 1")
            continue
        if compact:
            analyze_compact(g)
        else:
            analyze(g)


if __name__ == "__main__":
    main()
