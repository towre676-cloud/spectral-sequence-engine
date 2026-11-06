## 0.5.0 — Reusable simplicial geometry and twisted circle bundles (November 6, 2026)

- General finite simplicial-map engine accepting domain/base generating facets and a checked vertex map; canonical ordered simplices, degeneracies and induced base-skeletal filtrations.
- Concrete inverse-image subcomplexes, fibers, independent relative E1 computation, oriented filtered chain-map comparison and page naturality.
- Staircase simplex-product triangulations, nontrivial mapping tori, and independently computed circle local-system E2 from homology monodromy.
- Geometric Hopf refactor: compare the published reference filtration and automatically generated simplex births; verify Hopf S3 and trivial S1×S2 have identical 168-simplex f-vectors and E2 but rank-one vs rank-zero d2.
- Torus/Klein and twisted S2 bundle calculations: integral Z/2 torsion and mod-2 Bockstein vs odd-characteristic behavior; connected two-fold covering by exchanging discrete fiber sheets.
- Human-readable and tagged JSON inputs; deterministic replay certificates; CLI, 20+ new tests, 100 seeded random-map unit checks and a separate 1,000-map reproducible audit, experimental reporting and explicit research boundaries.
- Does **not** claim general Serre fibrations, integral filtered pages, Adams computations, or formal proof checking.

## 0.4.0 — Hopf geometry, lifted Bocksteins, graded couples, double complexes (November 6, 2026)

- Independently reproduced the published simplicial Hopf-fibration spectral-sequence example: 36 tetrahedra, 168 faces, exact E0/E1/E2/E3, nonzero rank-one d2, and total S3 homology over four prime fields.
- Exact lifted Bockstein pages and actual higher beta_r matrices, independently cross-checked against integral p-primary torsion invariants, with SHA-256 receipts and deterministic replay.
- Independent double-complex E1/E2 calculation and total-chain convergence, two filtrations, tensor/Koszul sign adapter, and hidden d2–d5 staircase regressions.
- Singly graded exact couples with checked degree shifts, derivative maps, and nonzero hidden graded d2.
- Added bounded multi-degree randomized integral tests, published geometric fixture, research limitations and provenance documentation. No claim to general geometric Serre/Adams or integral filtered spectral sequences.

## 0.3.0 — Additional exact algebraic machinery (November 6, 2026)

- Finite free integral chain-complex total homology, including torsion invariant factors, without third-party runtime dependencies.
- Exact-couple validation and iterative derived couples with a nonzero derived d₂ witness.
- Filtered associative DGAs over prime fields, chain and page products, graded Leibniz, associative/commutative checks, and representative-independence verification.
- Expanded tests, packaging, CI compilation, documentation, and honest scope boundaries.

# Change log

## 0.2.0 — exact filtered chain-complex expansion (proposed)

- Extended the original F₂ simplicial/persistence prototype into a general
  finite filtered chain-complex spectral-sequence engine over prime fields.
- Added actual page quotients, higher differentials, naturality maps,
  finite convergence checks and independent barcode comparison.
- Added reproducible JSON receipts, CLI examples, documentation and
  deterministic regression tests.
- Preserved the original `filtered_complex.FilteredComplex` API.

Not an independent formal verification or published package release.
