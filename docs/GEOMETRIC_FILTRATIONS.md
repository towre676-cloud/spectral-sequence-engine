# Finite simplicial-map spectral sequences (v0.5.0)

## Mathematical construction and limits

Let `f: K -> B` be a finite simplicial map specified by a vertex map. We
allow degeneracies: two distinct vertices of a simplex may have the same
image. For each simplex `sigma` let `f(sigma)` be the set of its **distinct**
image vertices. The map is simplicial exactly when `f(sigma)` is a simplex
of `B` for every `sigma` of `K`.

Set

`F_p K = {sigma in K : dim f(sigma) <= p}`, with `F_-1 K = empty`.

**Proposition.** These `F_p K` are subcomplexes, increasing in `p`, and
`F_dim(B) K = K`. Indeed, for a face `tau` of `sigma`, the image vertices
of `tau` are a subset of those of `sigma`, so
`dim f(tau) <= dim f(sigma)`. The oriented simplicial boundary therefore
preserves the filtration over every coefficient field.

The exact filtered chain-complex engine constructs

`E^1_{p,q} = H_{p+q}(F_p K, F_{p-1} K; F_l)`

and differentials `d_r: E^r_{p,q} -> E^r_{p-r,q+r-1}` (homological increasing
filtration). Finite filtration guarantees the associated graded total
homology is recovered on the stable page, subject to the exact computation.
The map itself also induces a **filtered chain map**
`C_*(K) -> C_*(B)` where the target `C_*(B)` has a skeletal filtration.
For a simplex whose image loses vertices the chain image is zero; for an
injective vertex restriction it is the target simplex with the signed
orientation permutation. Chain-map commutation and all-page naturality
are checked, including over odd characteristic.

**Critical caveat.** An arbitrary simplicial map need not be a Serre
fibration, nor is its `E2` page generally equal to
`H_p(B; H_q(fiber))`. The latter requires appropriate geometric
hypotheses, the homology local system and comparison between the skeletal
filtration and a Serre construction. For an arbitrary map we claim only
its **base-skeletal filtration spectral sequence** and convergence to
`H_*(K)`. Even a point mapping into one vertex of a circle supplies a
simple counterexample to inferring the tensor-product `E2` page from the
base and a naive fiber.

## Input, APIs and reproducibility

The module `simplicial_maps.py` accepts generating **facets** of `K` and
`B`, automatically generates all nonempty faces, rejects duplicate facets,
invalid vertices, partial vertex maps and images which are not base
simplices, then constructs exact oriented chain complexes. Vertices may
be integers, strings, or nonempty nested tuples thereof. The plain JSON
schema uses integer/string vertices; certificate JSON represents nested
tuples with explicit type tags.

```python
from simplicial_maps import SimplicialComplex, SimplicialMap
K = SimplicialComplex([('x','y','z')])
B = SimplicialComplex([('a','b')])
f = SimplicialMap(K, B, {'x':'a', 'y':'b', 'z':'b'})
assert f.births[('x','y','z')] == 1
assert f.chain_complex(prime=3).verify()
assert f.relative_e1_page(prime=3) == f.chain_complex(3).page(1)
assert f.chain_map(3).verify_naturality()
```

The independently calculated `relative_e1_page` recomputes homology of
the same-birth boundary blocks; it is **not** inferred from the higher-page
quotient algorithm. `preimage_skeleton(p)`, `fiber_over_vertex(v)` and
`preimage_subcomplex(L)` return concrete simplicial subcomplexes. The
`SimplicialComplex.integral_chain_complex()` method builds the fully
signed **unfiltered** integral chain complex for the existing torsion and
Bockstein engines.

```sh
python -m simplicial_maps examples/hopf_simplicial_map.json --prime 3 \
  --certificate /tmp/hopf-map-receipt.json
python -m simplicial_maps --replay /tmp/hopf-map-receipt.json
python -m simplicial_maps --example product-sphere --prime 3
python -m simplicial_maps --example klein --prime 3
python -m research.geometric_benchmarks
python -m research.random_geometry_audit --trials 1000 --seed 20261106
```

The receipt contains the exact facet/vertex-map input, f-vectors,
filtration simplex counts, successive page dimensions, differential ranks,
Betti numbers, a checksum of the underlying filtered-chain certificate,
and a SHA-256 digest. Replay reconstructs the complexes and recomputes all
reported numbers. These are deterministic **computational receipts**, not
cryptographic signatures or machine-checked formal proofs.

## Independent geometric experiments

### Hopf versus the product S1 × S2

The published 36-tetrahedron Hopf map is a triangulated `S3 -> S2`
with a circle fiber. Our geometric adapter reconstructs its 168 simplices,
including f-vector `(12,48,72,36)`. The fixture independently enumerates
faces and verifies that **every simplex** has the same induced filtration
birth as the old hand-written Hopf reference. The original mathematical
source is K. V. Madahar and K. S. Sarkaria, "A minimal triangulation of the
Hopf map and its application," *Geometriae Dedicata* (2000), reproduced in
[Macaulay2 documentation](https://www.macaulay2.com/doc/Macaulay2-1.25.11/share/doc/Macaulay2/SpectralSequences/html/___Computing_spthe_sp__Serre_sp__Spectral_sp__Sequence_spassociated_spto_spa_sp__Hopf_sp__Fibration.html).

The trivial `S1 × S2 -> S2` map is **generated automatically**, not copied
from the Hopf triangulation. The standard staircase triangulation of each
product of simplices yields exactly 168 simplices and the same
`(12,48,72,36)` f-vector. For **both** maps and each tested field
`F2,F3,F5` the `E2` page is one-dimensional in precisely
`(0,0),(0,1),(2,0),(2,1)`. But the computed
`d2: E2_(2,0) -> E2_(0,1)` is rank one for Hopf and zero for the product.
The results have the expected total homology:

| Map | rank d2 | nonzero Betti numbers |
| --- | --- | --- |
| Hopf S3 -> S2 | 1 | b0=1, b3=1 |
| Trivial S1 × S2 -> S2 | 0 | b0=b1=b2=b3=1 |

This is a standard topological distinction, *not* a new theorem. It is a
research-grade benchmark that requires the engine to construct and detect
an actual higher transgression from geometry.

### Mapping tori with genuine local systems

For any finite simplicial automorphism `g: F -> F`, we construct a mapping
torus over an `n`-gon (`n>=3`) by triangulating each simplex prism using
staircase paths and gluing the final top face by `g`. This construction
is a combinatorial model of a **bundle** over the circle.

For `V_q = H_q(F;F_l)` with monodromy `T_q=g_*`, the circle local-system
calculation is

`E2_(0,q) = coker(T_q-I)`, `E2_(1,q) = ker(T_q-I)`.

The finite base has dimension 1, so all `d_r` with `r>=2` vanish.
The module independently derives the monodromy matrix on actual homology
classes via an oriented chain map, calculates both invariant and
coinvariant dimensions, and compares that result to the exact geometric
filtration `E2` page. **These local-system calculations are certified for
the constructed mapping-torus family, not for arbitrary simplicial maps.**

- Identity automorphism on `S1`: torus, `H1(Z)=Z², H2(Z)=Z`.
- Reflection on `S1`: Klein bottle, `H1(Z)=Z ⊕ Z/2`, `H2(Z)=0`.
- Swap of two points: connected double covering of the base circle.
- Reflection on `S2`: twisted `S2`-bundle over `S1`, `H2(Z)=Z/2`,
  `H3(Z)=0`, versus the untwisted `S2 × S1` product.

The circle-reflection monodromy is `-1` on `H1(S1)`. Thus in odd
characteristic the Klein bottle has only `(0,0)` and `(1,0)` on its `E2`
page; in characteristic two the sign becomes `+1` and four entries remain.
This illustrates why suppressing `pi1`-action on the homology of the
fiber is invalid even for a genuine bundle.

The **integral** and **Bockstein** layers provide a further independent
check: the Klein bottle has a nonzero rank-one first mod-2 Bockstein
`H2(K; F2)->H1(K; F2)` caused by the `Z/2` torsion, whereas the torus
has no such differential. Thus the two spaces, despite identical ordinary
mod-2 Betti numbers, are distinguished by an integral-chain computation.

## Algorithmic limits, verification and next science

This implementation is exact, dense and intended for *finite and moderate*
triangulations. Facet closure may be exponential in the dimension of a
facet; page calculations require finite-field elimination; integral
reductions may grow large intermediate coefficients. No memory/runtime
promise is made for unrestricted input.

The regression suite checks: original published Hopf ranks in four
fields; product comparison; square triangulation; oriented signs and
collapsed simplices; subcomplex inverse images; monodromy at odd and even
primes; independent integral total homology and Bockstein; twisted S2
mapping torus; SHA-256/replay and deliberate tampering; and deterministic
random valid and invalid simplicial maps. `research/random_geometry_audit.py`
replays a separately seeded 1,000-map sweep, records its source SHA-256
fingerprint and checks all-page identities and naturality. These checks are not a formal
proof of software correctness or any new general theorem.

**Next steps:** chain-level certificates checked by a separate minimal
checker; computational detection/verification of fibrations (distinct from
merely testing map compatibility); general local coefficient systems over
higher-dimensional bases with explicit transport; rigorous filtered
cohomology products and cup-square structure; sparse/bounded reductions;
and a route from finite cell data to exact-couple or Serre constructions
under explicit hypotheses. General integral filtered spectral sequences,
extensions and Adams computations remain open.
