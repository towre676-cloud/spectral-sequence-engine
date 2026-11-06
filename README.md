# Spectral Sequence Engine

**Exact, inspectable spectral sequences for finite filtered chain complexes.**

A Python research toolkit for moving from a filtered boundary operator to genuine
spectral-sequence pages and differentials, with rigorous finite-field arithmetic,
persistence cross-checks, and portable computational replay receipts.

**What it calculates now (v0.5.0):**

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
| Higher Bockstein spectral sequence from integral chain lifts | Implemented, with independent prime-primary torsion-rank cross-checks |
| First-quadrant double complexes and both filtration spectral sequences | Implemented, with independently computed E1/E2 |
| Singly graded exact couples and their derived differentials | Implemented |
| **Arbitrary finite simplicial-map adapters** and automatic base-skeletal filtrations | **Implemented**, with degenerate simplex handling and oriented chain-map naturality |
| Inverse images of base skeleta, vertex fibers and subcomplexes | Implemented and checked |
| Independent `E1` relative-homology computation from geometric same-level boundaries | Implemented |
| Staircase triangulations of finite simplicial products | Implemented |
| Simplicial mapping tori and circle-base local coefficient systems | Implemented and independently checked using fiber monodromy |
| Total integral simplicial homology and Bockstein of constructed spaces | Implemented (distinct from integral filtered pages) |
| Serre-style Hopf fibration from a published 36-tetrahedron simplicial map | Reproduced over F2/F3/F5/F7 |
| General integral filtered pages, arbitrary geometric Serre/Adams sequences, extension problems | **Not implemented** |

## v0.5.0: geometry-to-spectral-sequence pipeline

Give the engine an **arbitrary finite simplicial map** (domain facets, base facets
and a vertex map) and it constructs the base-skeletal filtration automatically.
It validates simplex images, checks oriented chain maps and graded naturality,
computes actual `d_r` matrices, and independently reconstructs `E1` relative
homology from the filtration blocks. These calculations are **filtration
spectral sequences**; without fibration hypotheses we make **no generic Serre
spectral sequence claim**.

```python
from simplicial_maps import SimplicialComplex, SimplicialMap
K = SimplicialComplex([('x','y','z')])
B = SimplicialComplex([('a','b')])
f = SimplicialMap(K, B, {'x':'a', 'y':'b', 'z':'b'})
assert f.chain_complex(3).verify()
assert f.chain_map(3).verify_naturality()
assert f.relative_e1_page(3) == f.chain_complex(3).page(1)
```

```sh
python -m simplicial_maps examples/hopf_simplicial_map.json --prime 3 \
    --certificate geometric-receipt.json
python -m simplicial_maps --replay geometric-receipt.json
python -m simplicial_maps --example klein --prime 3
python -m research.geometric_benchmarks
python -m research.random_geometry_audit --trials 1000 --seed 20261106
```

The **Hopf fibration** and the independently constructed **trivial projection**
`S1 × S2 -> S2` have the same 168-simplex f-vector and the same `E2`
dimensions, but their calculated `d2` ranks are **1 and 0**, respectively.
The geometric mapping-torus constructor additionally distinguishes the torus
and Klein bottle through nontrivial `pi1`-action on fiber homology. The Klein
bottle has `H1(Z) = Z ⊕ Z/2`, confirmed independently using integral simplex
boundaries and a rank-one mod-2 first Bockstein; an orientation-reversing
`S2` mapping torus yields `H2(Z) = Z/2`. The mod-2 and odd-prime page
comparisons agree with explicit monodromy invariant/coinvariant calculations.

See [geometric definitions, exact benchmark data and scope boundaries](docs/GEOMETRIC_FILTRATIONS.md).

## Research landmark: the actual triangulated Hopf fibration

Our `research/hopf_fibration.py` passes the published simplicial Hopf map
`S³ → S²` (with fiber `S¹`) to the general-purpose geometric adapter and
retains an independent filtration reference. The adapter closes its 36
tetrahedral facets under faces (168 simplices) and filters the total space
by inverse images of the **base skeleta**.
Without manually entering page differentials, the engine obtains exactly the
four expected one-dimensional E2 positions `(0,0), (0,1), (2,0), (2,1)`, a
**nonzero d2** from `(2,0)` to `(0,1)`, and E3 containing just `(0,0)`
and `(2,1)`. This reproduces the E0/E1/E2/E3 **dimensions** reported in the
[Macaulay2 worked Hopf-fibration example](https://www.macaulay2.com/doc/Macaulay2-1.25.11/share/doc/Macaulay2/SpectralSequences/html/___Computing_spthe_sp__Serre_sp__Spectral_sp__Sequence_spassociated_spto_spa_sp__Hopf_sp__Fibration.html)
from actual simplicial data, now checked over F2/F3/F5/F7.

```sh
python -m research.hopf_fibration
python -m unittest tests.test_hopf_fibration -v
```

The [research report](docs/RESEARCH_2026_11_06.md) documents provenance,
mathematical interpretation, exact algorithms, limitations and next targets.
This is a **specific geometric example**, not an implementation of general
Serre spectral sequences or discovery of a new Hopf-fibration theorem.

### Higher Bockstein differentials from integer lifts

```python
from integral import IntegralChainComplex
from bockstein import LiftedBockstein
c = IntegralChainComplex([1, 1], {1: [[16]]})
b = LiftedBockstein(c, prime=2)
assert b.differential(1)[1] == ((0,),)
assert b.differential(4)[1] == ((1,),)
assert b.page(5) == {}
assert b.verify()
```

The Bockstein engine constructs actual integer lifts and beta_r maps and
cross-checks their ranks against independently calculated integral torsion.
See `bockstein.py` for the separate rank-only interface and deterministic
certificate/replay. `bicomplex.py` and `graded_couple.py` supply new structural
machinery for double-complex and graded exact-couple calculations.

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
DGA page products, and ungraded derived exact couples respectively. Version
0.4.0 added `bockstein.LiftedBockstein`, `bicomplex.DoubleComplex`,
`bicomplex.tensor_double_complex`, and `graded_couple.GradedExactCouple`.
Version 0.5.0 adds `simplicial_maps.SimplicialMap`, `simplicial_product`,
`mapping_torus_projection`, and `circle_local_system_e2`. See
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

The primary filtered spectral-page engine is **finite and field-valued**.
The general geometric-map adapter produces a base-skeletal filtration spectral
sequence, **not necessarily** a Serre spectral sequence; only the designated
Hopf, product and mapping-torus examples carry the specified geometric bundle
interpretations. Its separate
integral routines compute total homology over Z and a prime-field Bockstein spectral
sequence from that integral data, **not general integral filtered spectral pages**. Its exact equalities
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
