# Mathematical specification

All chain groups are finite-dimensional over F_ℓ for **prime** ℓ. The
increasing filtration is by subcomplexes `... ⊂ F_(p-1) C ⊂ F_p C ⊂ ...`.
The boundary preserves filtration, decreases homological degree by one, and
satisfies `∂² = 0` exactly modulo ℓ. We set `F_p = 0` below the minimum grade.

## Pages as actual quotient spaces

For `r ≥ 0`, define

```
Z^r_{p,n} = { x ∈ F_p C_n : ∂x ∈ F_(p-r) C_(n-1) }
B^r_{p,n} = F_p C_n ∩ ∂(F_(p+r) C_(n+1)).
```

Then the implementation represents (for `r ≥ 1`)

```
E^r_{p,q} = Z^r_(p,p+q) /
            ( Z^(r-1)_(p-1,p+q) + B^(r-1)_(p,p+q) ).
```

A basis is found using modular kernel elimination, intersection with a
filtration subspace, and elimination modulo an explicitly constructed
*denominator* span. No floating-point rank thresholds are used.

At `r = 0`, `E⁰` is the associated graded chain group `F_p/F_(p-1)`
and `d₀` is the same-grade part of `∂`. At `r ≥ 1`, apply the *original*
`∂` to a quotient representative, then express its image in the target's
boundary subspace and quotient representatives. The quotient coordinates
are the matrix column of `dᵣ`.

The bidegree is `(-r, r-1)`:

```
               d_r
 E^r_(p,q) -------------> E^r_(p-r,q+r-1).
```

This is a **homological** indexing convention. Cohomological conventions
(reverse arrows/index signs) cannot be imported without conversion.

## Why the finite filtration stabilizes

A finite chain complex has minimum and maximum generator births, `m` and
`M`. A differential with `r > M-m` has no nonzero target of the required
filtration; hence there can be no additional page changes. The engine checks
stabilization at `E_(M-m+1)` (and at the next page) and compares

```
Σ_{p+q=n} dim E∞_(p,q) = dim H_n(C).
```

This is **graded-dimension convergence** for this finite filtered model;
it does not by itself identify an integral group extension, a product,
filtered chain equivalence, or a general unbounded convergent sequence.

## Independent persistence identity

Ordinary persistence uses a separate filtered boundary-column reduction
(not the page quotient algorithm). If `I=[birth,death)` is a finite bar
(including length-zero intervals), exact finite-field algebra gives

```
Σ_(p,q) rank d_r^(p,q) = #{finite persistence bars of length r}.
```

Each infinite interval gives an essential homology generator. These are
valuable independent consistency checks, but the same underlying boundary
input is shared, and the code is not formally verified.

## A hidden differential rather than a staged declaration

Take generators over F₃:

| Generator | Degree | Birth | Boundary |
|---|---:|---:|---|
| a | 0 | 0 | 0 |
| b | 0 | 1 | 0 |
| u | 1 | 1 | a+b |
| v | 1 | 2 | b |

Here `d₀(u) = b` on the associated graded, so the grade-one classes
cancel. The grade-two generator `v` is not a `d₁` cycle on the nose,
but its representative can be corrected by `v-u`, giving

```
∂(v-u) = -a ≠ 0.
```

Thus the corrected representative supports a nonzero `d₂` from the
filtration-two degree-one class to the filtration-zero degree-zero class.
The engine finds this automatically from `Z²` and the quotient model.

## What comes next

The natural mathematical expansions are sparse matrix kernels and quotient
representatives, integral Smith normal form with torsion/extension
tracking, multiplicative spectral sequences, filtered cochain complexes,
exact couples, Mayer–Vietoris, and external/independently checked proof
certificates. Each requires additional mathematical structure beyond the
current finite prime-field linear algebra.
