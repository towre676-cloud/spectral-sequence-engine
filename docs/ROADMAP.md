# Roadmap and proof-authority ledger

## Landed in this expansion

1. Finite filtered chain complexes over prime fields with checked `∂²=0`.
2. Full chain-level `E₀`, `E₁`, and `Eᵣ` quotient spaces and actual `dᵣ`
   matrices, including representative corrections and higher-page cases.
3. Chain maps checked for filtration and chain commutation, with page maps
   checked for naturality.
4. Ordinary persistence reduction as an independent structural cross-check.
5. Computation-replay JSON, CLI, examples, repeatable random regression tests.
6. Compatibility entry points for the original simplicial F₂ API.

## Landed in v0.3.0 and v0.4.0

- Finite free integral total homology with torsion invariant factors.
- Un-graded and singly graded exact couples with genuine derived differentials.
- Prime-field filtered differential graded associative algebras and page products.
- First-quadrant bicomplexes, two filtrations, independent E1/E2 checks, and tensor-complex adapters.
- Bockstein spectral sequences from finite free integral complexes, including actual beta_r matrices from integer lifts and independent torsion fingerprints.
- A published 36-tetrahedron simplicial Hopf-map calculation reproducing the nonzero Serre d2 and S3 total homology over F2/F3/F5/F7.

## Next: structural mathematics

- Filtered mapping cones, distinguished long exact sequences, exact couples,
  and independently checked connecting morphisms.
- Mayer–Vietoris filtration from simplicial covers, and examples where
  the two independent constructions agree.
- Efficient sparse page reductions and memory caps with benchmark receipts.
- Multiplicative extension and geometric cup-products **only when a filtered
  differential graded algebra and diagonal approximation are actually supplied**.
- Integral **filtered** spectral pages and torsion-aware filtration extensions,
  clearly separating associated graded groups from reconstructed groups; total
  integer homology and mod-p Bocksteins are now implemented separately.
- Functorial invariants for filtered chain homotopies and cone comparison.

## Next: applications and external evidence

- Explicit chain-level adapters for cell complexes, cubical complexes,
  selected filtered complexes in network science, and crystallographic
  cell complexes with trustworthy incidence data.
- Compare computed field-persistence output against independent packages
  such as GUDHI or Ripser where their supported domains overlap.
- Independent verifier and formal kernel for finite-field matrix operations,
  followed by public cross-platform reproducibility receipts.

## Not established

No general cohomological descent, sheaf/motivic obstruction theory,
Adams/Adams–Novikov computation, all-space convergence theorem, extension
classification, physical prediction, or cryptographic security property
follows from this release. Such claims require separate objects, theorems,
proofs, fixtures, and validations.
