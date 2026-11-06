# Extensions: integral homology, exact couples, and multiplicative pages

## Mathematical scope

The repository now contains three **distinct**, deliberately bounded computation layers.

1. `integral.IntegralChainComplex` computes the **total homology** of a finite free chain complex over **Z**, including free ranks and torsion invariant factors. Integer row and column operations are unimodular. The outgoing boundary is diagonally reduced while tracking a full integer kernel basis, after which the incoming boundary is expressed in that kernel and reduced to cyclic factors. A gcd/lcm normalization reports invariant factors.
2. `exact_couple.ExactCouple` represents a **finite-dimensional ungraded** exact couple over a prime field. It validates all three exactness equations and can iterate derived exact couples. The differential on a given E term is `j ∘ k`; the first derived E is `ker(jk) / im(jk)`. Subsequent differentials arise from the derived maps, rather than supplied page data.
3. `multiplicative.FilteredDGA` validates an associative filtered differential graded algebra over a prime field, including graded Leibniz, optional graded commutativity, degree, and filtration. It computes induced products from the exact chain representatives of `FilteredChainComplex` on each `E_r` and checks descent across page denominators and the page-level Leibniz identity.

**Do not confuse the three layers.** The integral routine does not compute integral spectral sequences. The exact couple is ungraded; filtration-index shifts and geometric connecting maps have not been built. The multiplicative engine computes products from a user-supplied DGA and does not invent the cup-product data of a topological space. No comparison with an independent computer algebra implementation is being claimed.

## A genuinely hidden derived-couple differential

Use basis `a,b` for `D` and `x,y` for `E` over F₃:

- `i(b)=a`, `i(a)=0`;
- `j(b)=x`, `j(a)=0`;
- `k(y)=a`, `k(x)=0`.

The original `d₁=jk=0`, while the next exact couple has `D'=span(a)` and `E'=span(x,y)` with `j'(a)=x`, `k'(y)=a`. Therefore `d₂'(y)=x`, and the third page vanishes. The implementation verifies exactness of each derived couple, rather than only checking dimensions.

```python
from exact_couple import ExactCouple
c = ExactCouple(3, 2, 2,
    [(0,0), (1,0)], [(0,0), (1,0)], [(0,0), (1,0)])
assert [page['d_rank'] for page in c.derived_pages()] == [0, 1, 0]
```

## Integral torsion (total homology, not spectral pages)

```python
from integral import IntegralChainComplex
# C₂ = Z --(6)--> C₁ = Z --(0)--> C₀ = Z
c = IntegralChainComplex([1,1,1], {2: [[6]]})
assert c.homology()[1] == {'free_rank': 0, 'torsion': (6,)}
```

For a product of two cyclic groups of coprime orders, the routine normalizes them to a single invariant factor. Empty `torsion` denotes a torsion-free group, and `free_rank` is the number of copies of Z. Because coefficients can grow during integer row and column operations, large matrices require a production-grade Smith algorithm and additional benchmarks.

## Multiplicative exterior algebra

```python
from spectral_sequence import FilteredChainComplex
from multiplicative import FilteredDGA
c = FilteredChainComplex([
    ('one',0,0,{}), ('x',1,1,{}), ('y',1,1,{}), ('xy',2,2,{})
], prime=3)
products = {('one',v):{v:1} for v in c.names}
products.update({(v,'one'):{v:1} for v in c.names})
products[('x','y')] = {'xy':1}
products[('y','x')] = {'xy':-1}
a = FilteredDGA(c, products, graded_commutative=True)
assert a.page_product(1, 1, 1, 1, 1) == ((0,), (1,), (2,), (0,))
assert a.verify_page_well_defined()
assert a.verify_page_leibniz()
```

## Next research and nonclaims

- General **graded** exact couples for Serre, Adams, and Bockstein constructions require correct multi-index shifts and geometrically constructed `i,j,k` maps.
- Integral spectral sequences require saturated subquotient and module presentation algorithms, not simply replacing Fₚ arithmetic with Z arithmetic.
- Multiplicative extension problems need chain-level geometrically meaningful cup products, filtration-compatible diagonals, and additional data beyond `E_∞`.
- Cross-check integral torsion and filtered multiplicative products with independent software; add property-based randomized tests and large sparse benchmark suites.
- The current checks are executable exact arithmetic, **not** formally verified proofs and **not** a general Serre/Adams spectral-sequence solver.
