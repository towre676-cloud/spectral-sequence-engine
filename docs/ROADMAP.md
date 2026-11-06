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

## Next: structural mathematics

- Filtered mapping cones, distinguished long exact sequences, exact couples,
  and independently checked connecting morphisms.
- Mayer–Vietoris filtration from simplicial covers, and examples where
  the two independent constructions agree.
- Efficient sparse page reductions and memory caps with benchmark receipts.
- Multiplicative structure and the `dᵣ` Leibniz rule **only when a filtered
  differential graded algebra is actually supplied**.
- Integral Smith normal forms and torsion-aware filtration extensions,
  clearly separating associated graded groups from reconstructed groups.
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
