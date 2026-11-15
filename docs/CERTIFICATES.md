# Computational receipts and trust boundaries

Calling `complex_.certificate()` runs all internal checks and produces a
JSON object (`spectral-sequence-replay-v1`) with:

- Normalized prime-field generators, filtration levels and exact differentials.
- SHA-256 of canonical normalized input JSON.
- Dimension of every nonzero `(p,q)` group on each page from `E₀` through
  the stable page.
- Explicit column-major matrices, targets and ranks for every `dᵣ`.
- Ordinary persistence intervals including zero-length intervals.
- Total homology Betti dimensions, stable-page index and check indicator.
- SHA-256 of the canonical JSON report (excluding its own digest field).

`replay_certificate(receipt)` validates the declared format and digest,
reconstructs the original filtered chain complex, recomputes the complete
report and requires byte-equivalent normalized mathematical content.

**Threat model:** SHA-256 catches accidental tampering; it is *not* a
digital signature, and no outside trust anchor attests to receipt authorship.
Replay uses the same Python implementation and is not a logically
independent checker or proof assistant. Exact finite-field arithmetic
prevents floating-point conditioning issues but does not remove the
possibility of implementation errors. The tests compare independent
barcode reduction with the quotient-page algorithm and include exhaustive
invariants and deterministic generated examples; they are not a proof of
correctness for arbitrary inputs.

An ideal future version will include a minimized independent checker and
formal definitions/certificates in Lean or Coq, with explicit size bounds
and validated canonical encoding of matrices.

## Integral homology class receipts

The CRO witness method from [τ-Crystal](https://github.com/towre676-cloud/tau_crystal/tree/82eef62af64e020b369d3bbfd8ad2c4626ac67a5/tools/cro)
also applies to any `IntegralChainComplex`. The receipt certifies a specified
chain in the given bases. Use the same complex and chain to replay a receipt.

```python
from integral import IntegralChainComplex, integral_mapping_cone

circle = IntegralChainComplex([0, 1], {})
cone = integral_mapping_cone(circle, circle, {1: [[4]]})
receipt = cone.class_certificate(1, [2])
assert receipt["status"] == "torsion" and receipt["order"] == 2
assert cone.verify_class_certificate(1, [2], receipt)
```

For a cycle z in degree n, put B = d_(n+1). A boundary receipt supplies
Bx = z. An infinite-order receipt supplies an integer functional ℓ with
ℓB = 0 and ℓz ≠ 0. For exact finite order h, a lift supplies Bx = hz.
Modular functionals ℓ_i, each killing B modulo m_i, give the lower bound
via lcm_i(m_i / gcd(m_i, ℓ_i z)) = h. Every annihilating multiplier is
divisible by each of those factors. A noncycle receipt gives a nonzero
coordinate of d_n z. The hash binds the complex, degree and coordinates.

Generation uses the existing unimodular diagonalization, with row
operations retained for witnesses. If U B V = D and u = Uz, a nonzero u
in a zero row gives an infinite-order witness. Otherwise the least
annihilating multiplier is the least common multiple of the divisibility
requirements from the nonzero diagonal rows. Replay checks the witness
equations directly without running the reduction or computing homology.

In the degree-four circle cone, filter the target generator a at level
zero and the shifted source generator b at level one so db = 4a.
The class of 2a has exact order two. It reduces to zero in characteristic
two and bounds over every odd prime field. Modulo four, a coordinate
functional detects it. The filtered field pages and the lifted mod-two
Bockstein are tested separately. The latter has zero β₁ and rank-one β₂,
which detects the ambient order-four summand. The integral receipt
distinguishes the individual classes a, 2a and 4a. It certifies total
integral homology of the cone and does not compute integral spectral pages.
