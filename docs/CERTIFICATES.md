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
