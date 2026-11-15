# Independent verification of integral class receipts

`integral_receipt_checker.py` can run alone with Python 3.10 or later. It
imports only the standard library. It accepts existing
`integral-homology-class/1` receipts and implements its own input
validation, normalization, and witness checks. PR #5's generator and
in-process verifier remain available through `IntegralChainComplex`.

## Grammar and input binding

The verifier takes the based complex, degree and chain as separate,
caller-supplied inputs. A CLI input JSON object has exactly these keys:

```json
{"ranks":[0,1,1],"boundaries":{"2":[[4]]},"degree":1,"chain":[2]}
```

Let `ranks[n]` be the nonnegative rank of the free chain group in degree
`n`, and let ranks outside the supplied array be zero. Boundary keys are
canonical decimal strings for degrees `1` through `len(ranks)-1` in JSON
and integer degrees in the Python API. Each matrix is row-major, with
`ranks[n-1]` rows of length `ranks[n]`. Omitted maps expand to zero
matrices of these shapes. Zero-row matrices have no rows even if the
domain is nonzero. The verifier checks every consecutive composition
and requires the input degree to be in the ranks array and the chain to
have exactly `ranks[degree]` coordinates.

The receipt is a JSON object with the three common keys below and exactly
the additional keys specified for its status. All numeric entries are
integers. The CLI rejects duplicate keys, floating-point literals and
non-finite numbers. Booleans are rejected wherever integers are required.
Vectors are JSON arrays or Python lists or tuples.

| Common key | Value |
| --- | --- |
| `kind` | Literal string `integral-homology-class/1` |
| `input_sha256` | Lowercase 64-hex-digit SHA-256 of the canonical input |
| `status` | One of the four strings below |

| Status | Additional fields | Equations |
| --- | --- | --- |
| `noncycle` | `row`, `value` | Zero-based row of \(d_nz\) with the stated nonzero integer value |
| `boundary` | `lift` | \(d_nz=0\) and \(Bx=z\) |
| `infinite` | `functional` | \(d_nz=0\), \(\ell B=0\), and \(\ell z\ne0\) |
| `torsion` | `order`, `lift`, `witnesses` | \(d_nz=0\), \(Bx=hz\), and modular image orders have lcm \(h\) |

Here \(n\) is the degree, \(z\) the input chain, \(B=d_{n+1}\),
\(x\) the lift, and \(h\) the integer order, with \(h\ge2\). Lift
length is `ranks[n+1]` and functional length is `ranks[n]`.
`witnesses` is a nonempty array of objects with exactly `modulus` and
`functional`. Each integer modulus \(m\ge2\) gives the congruence
\(\ell B=0\pmod m\) and the image order \(m/\gcd(m,\ell z)\).
Composite moduli, redundant witnesses, and signed coordinates are valid.

For the digest, serialize the array
`[ranks, normalized_boundaries, degree, chain]` as ASCII JSON without
whitespace, with normalized boundary keys ordered by **numeric** degree
before converting to JSON strings. Include every interior boundary map.
Use ordinary decimal integer encoding. In Python, this is
`json.dumps(payload, sort_keys=True, separators=(",", ":"))` with
integer keys in `normalized_boundaries`. Receipt fields and all metadata
outside this payload are outside the digest. Equivalent explicit and
implicit zero maps produce the same digest.

## Why the equations suffice

An integer row annihilating \(B\) descends to a map from homology to
\(\mathbb Z\), so a nonzero evaluation detects infinite order. A row
annihilating \(B\) modulo \(m\) descends to a map to \(\mathbb Z/m\).
If an integer multiplier kills \([z]\), it also kills each modular image
of \([z]\), so each image order divides that multiplier. For torsion, the
lift makes \(h\) an annihilating multiplier and the modular rows force
every annihilating multiplier to be divisible by \(h\). The checker
tracks the lcm in `covered` and keeps \(h=\texttt{covered}\cdot
\texttt{remaining}\). A new image order \(s\) increases the lcm by
\(s/\gcd(s,\texttt{covered})\). Ending with remaining equal to one
proves exactness without factoring any modulus or order.

The witness forms are complete for finite free integer complexes. Apply
a unimodular diagonalization \(UBV=D\) for this existence argument.
A nonzero coordinate of \(Uz\) on a zero row gives the infinite-order
witness. Otherwise each nonzero diagonal row forces divisibility by
\(|D_{ii}|/\gcd(D_{ii},(Uz)_i)\). The lcm clears all diagonal equations,
which yields an integer lift through \(V\). The corresponding rows of
\(U\) give the modular witnesses. The checker itself uses only products,
integer comparison, division, gcd and hashing.

Acceptance establishes the status of a specific chain in the supplied
based complex, assuming correct execution of the checker and Python's
integer operations. Caller provenance determines whether that complex and
chain describe the intended object. SHA-256 binds the normalized input and
provides no authentication or signature. This is an independently
implemented Python checker with a mathematical justification. Machine
formalization, general integral page certificates, and resource ceilings
for hostile input are outside this implementation.

## Reproduction

```sh
python -m unittest discover -s tests -p 'test_integral_receipt_checker.py' -v
python -m unittest discover -s tests -v
python -m integral_receipt_checker input.json receipt.json
```

The CLI exits 0 for a valid receipt, 1 for a rejected receipt, and 2 for
unreadable files or invalid JSON or base input. The Python entry point is
`verify_class_receipt(ranks, boundaries, degree, chain, receipt)`. It
raises `ValueError` for invalid base data and returns a Boolean for
receipts decoded to ordinary Python values.

The tests check generated receipts with both verifiers, disable the
generator's reduction during replay, and run the CLI from a directory
containing only the checker and data. An independent oracle uses gcds of
integer minors in explicit coordinates on the cycle lattice. Nonzero
determinantal divisors recover the invariant factors and free rank.
For a cycle, adjoining its column increases the rank exactly when it
has infinite order. At equal ranks, the ratio of the final nonzero
determinantal divisors before and after adjoining the column is the
index of the boundary lattice in the enlarged lattice, hence the exact
class order. The test uses seed 20261115 with 96 sampled cycles in four
non-diagonal presentations, alongside deterministic checks.

## Filtered integral complexes

For a finite increasing filtration by based subcomplexes, a witness in
\(F_p C\) certifies a class there. A witness in \(C\) for the same cycle
certifies its image under inclusion. A boundary lift in \(F_p C\) proves
the class already dies by level \(p\). To certify the first death level
in a discrete filtration, add a nonzero class witness at the preceding
level, verifying the same cycle and the inclusion between the two
levels. Relative and associated-graded homology can use the same
checker when the quotient has a verified free basis, as with a filtration
by subsets of the chain basis. A filtered receipt must bind the
filtration, level and inclusions or quotient maps used to derive the
based complex.

For higher pages, use total degree \(n=p+q\) and define
\(Z_r^{p,n}=\{x\in F_p C_n: dx\in F_{p-r}C_{n-1}\}\). For \(r\ge1\),
the page term has numerator \(Z_r^{p,n}\) and denominator
\(Z_{r-1}^{p-1,n}+d Z_{r-1}^{p+r-1,n+1}\). This suggests reusing the
witness equations in a verified presentation of that quotient. A
complete presentation certificate must show that the proposed integral
basis spans the full numerator lattice and that the relation columns
generate the entire denominator. One route is to certify integral
kernel bases with unimodular changes of coordinates and verify the
projected denominator generators. Merely checking that columns are
cycles or valid relations can miss an integral index. A separate
verification of the induced differential must check that representatives
map by \(d\), that relations map to relations, and that the target
coordinates are valid. No such presentation or differential checker is
implemented here.

The extension issue persists even if every page is verified. In degree
one take generators \(a,b\) at filtration levels zero and one, with
degree-two generators \(x,y\) at the respective levels. Let
\(dx=2a\), \(dy=2b+\epsilon a\), with \(\epsilon\in\{0,1\}\).
Both associated-graded complexes give \(E_1^{0,1}=E_1^{1,0}=\mathbb Z/2\)
with all other terms zero. Every higher differential lowers total
degree by one and hence vanishes. When \(\epsilon=0\), total \(H_1\)
is \((\mathbb Z/2)^2\); when \(\epsilon=1\), it is \(\mathbb Z/4\)
with \(2[b]=-[a]\). The test certifies the two respective orders of
\([b]\). Total receipts can check a proposed lifted class or a proposed
relation by proving that its difference bounds. Reconstructing the
whole filtered homology group also requires certified generating and
relation data across the filtration.
