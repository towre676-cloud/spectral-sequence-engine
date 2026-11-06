# Spectral Sequence Engine

**Exact, inspectable spectral sequences for finite filtered chain complexes.**

A Python research toolkit for moving from a filtered boundary operator to genuine
spectral-sequence pages and differentials, with rigorous finite-field arithmetic,
persistence cross-checks, and portable computational replay receipts.

**What it calculates now:**

| Calculation | Status |
|---|---|
| Filtered simplicial complexes and ordinary persistence | Implemented |
| Arbitrary finite filtered chain complexes over **Fₚ** (prime p) | Implemented |
| Associated graded `E₀` and relative-homology `E₁` | Implemented |
| Actual `d₀, d₁, d₂, …` matrices and successive `Eᵣ` pages | Implemented |
| Exact quotient representatives, including hidden higher differentials | Implemented |
| `dᵣ² = 0` and `Eᵣ₊₁ = H(Eᵣ,dᵣ)` | Checked |
| Finite-filtration convergence and barcode-length identities | Checked |
| Filtered chain maps and `dᵣ` naturality | Implemented and checked |
| Deterministic JSON receipts with tamper and replay checks | Implemented |
| Total integral homology with free and torsion invariant factors | Implemented (separate chain-complex module) |
| Prime-field multiplicative pages of a supplied filtered DGA | Implemented and checked |
| Ungraded exact-couple validation and successive derived couples | Implemented |
| Integral spectral pages, graded geometric exact couples, cohomology spectral sequences, extension problems | **Not implemented** |

## See a nonzero higher differential

```sh
python spectral_sequence.py --example d3
```

For a degree-one cell born at filtration 3 whose boundary is a degree-zero
cell born at 0, the generator survives `E₁`, `E₂`, and `E₃`, then the exact
`d₃` matrix has rank 1, so `E₄` is zero.

```text
E_1: (0,0):1, (3,-2):1; rank(d_1)=0
E_2: (0,0):1, (3,-2):1; rank(d_2)=0
E_3: (0,0):1, (3,-2):1; rank(d_3)=1
E_4: (zero); rank(d_4)=0
Exact verification: PASS
```

The **hidden-d₂** case is more interesting: `examples/hidden_d2.json` has a
same-filtration differential which hides the higher obstruction until a
representative is corrected. It is not a manually entered differential—the
engine derives it from the full chain complex.

## Start

Requires **Python 3.10+**; the runtime and tests use only the standard library.

```sh
python -m unittest discover -s tests -v
python spectral_sequence.py --example triangle
python spectral_sequence.py examples/hidden_d2.json --certificate receipt.json
python spectral_sequence.py --replay receipt.json
```

Optional installable console command:

```sh
python -m pip install -e .
spectral-sequence --example d2
```

## Python API

```python
from spectral_sequence import FilteredChainComplex

complex_ = FilteredChainComplex([
    ('a', 0, 0, {}),
    ('b', 1, 2, {'a': 1}),
], prime=3)

assert complex_.page(2) == {(0, 0): 1, (2, -1): 1}
assert complex_.differential(2)[(2, -1)]['rank'] == 1
assert complex_.page(3) == {}
assert complex_.verify()
receipt = complex_.certificate()
```

Each differential entry contains `target`, `columns` (column-major matrix in
computed quotient bases), and exact `rank`. Keys are `(p,q)` for filtration
`p` and complementary grading `q`.

The original `filtered_complex.FilteredComplex` API for simplices and binary
persistence remains available. Call `to_chain_complex(prime=3)` to upgrade a
simplicial filtration to signed oriented boundaries over **F₃**, or use
`page(r)` / `differential(r)` directly from the legacy object.

Further modules `integral.IntegralChainComplex`, `multiplicative.FilteredDGA`,
and `exact_couple.ExactCouple` provide integral **total** homology, filtered
DGA page products, and ungraded derived exact couples respectively. See
[extension examples and precise limitations](docs/EXTENSIONS.md).

Filtered chain maps are in `morphisms.FilteredChainMap`, which rejects maps
that fail to preserve degree, filtration, or the chain-map square, and checks
naturality on every computed page.

## Definitions and verification

We use an **increasing homological** filtration:

\[
E^1_{p,q}=H_{p+q}(F_p C/F_{p-1} C), \qquad
 d_r:E^r_{p,q}\longrightarrow E^r_{p-r,q+r-1}.
\]

For each bidegree the engine constructs subspaces of the **original** chain
groups, not just dimensions. It computes exact quotients and induced maps.
Its independent barcode reduction must recover the same total differential
ranks as the number of finite intervals of length `r`. It also checks
`E_∞`'s graded Betti ranks against the homology of the total complex.

See [the mathematical conventions](docs/MATH.md),
[replay and trust boundaries](docs/CERTIFICATES.md), and
[the expansion roadmap](docs/ROADMAP.md).

## Scope and research integrity

The primary spectral-page engine is **finite and field-valued**. Its separate
integral routine computes total homology over Z, not integral spectral pages. Its exact equalities
are derived by machine-executed modular linear algebra and deterministic
replay—not a Lean/Coq proof, and not a solution to arbitrary spectral-sequence
convergence or extension problems. Complexity grows with the number of
cells and filtration span; large sparse complexes will need optimized
algorithms and memory controls.

Historical exploratory claims about sheaves, motivic phenomena, physics,
cryptography, or "Lorentz differentials" are **not** accepted inputs or
established results. Building trustworthy mathematical primitives comes first.

Contributions should include a minimal example, a test demonstrating the
mathematical claim, and a clear statement of what has and has not been certified.
