# Spectral Sequence Engine

Initial exact computational topology foundations over the field **F₂**.
This is a bounded research prototype, not a general spectral-sequence solver.

The core represents finite filtered simplicial complexes, verifies the
filtration and face conditions, computes one-parameter persistent-homology
intervals by boundary-column reduction, and computes the dimensions of
`E¹_{p,q} = H_{p+q}(F_p,F_{p-1}; F₂)` through associated-graded boundary matrices.

## Reproduce

```sh
python -m unittest discover -s tests -v
```

The filled-triangle fixture has Betti ranks:
- level 0: β₀=3;
- level 1: β₀=1, β₁=1;
- level 2: β₀=1.

The model checks these against persistence intervals and the E¹ page. It
does **not** compute the first differential d¹, higher differentials, page
convergence, extension problems, sheaf cohomology, or multiparameter persistence.

## Next mathematics

1. Implement d¹ and prove d¹d¹=0 with indexed relative representatives.
2. Implement page-to-page exact quotient and lifting for dʳ.
3. Verify convergence against the total homology on complete finite examples.
4. Introduce Mayer–Vietoris filtrations and cross-check independent software.
5. Add certificate exports containing rational/F₂ matrices and replay receipts.
6. Only then connect to structural network or crystal applications.

The April 2024 historical excavation motivated this project; speculative
"Lorentz differentials", unspecified sheaves and fairness claims from those
conversations are expressly **not** mathematical inputs.
