"""Integral homology of finite free chain complexes using unimodular reductions.

Integral homology and class certificates. This does NOT compute integral spectral pages,
Bockstein differentials, hidden extensions, or coefficients in a general ring.
No external dependencies. Reduction is cubic-ish for small examples but integer
coefficient growth can be prohibitive on large complexes.
"""
from fractions import Fraction
from hashlib import sha256
import json
from math import gcd, lcm


def _diagonalize(matrix, cols, track_rows=False):
    """D -> U D V diagonal, tracking V and optionally U.

    The diagonal need not be Smith-divisibility ordered. This suffices to
    calculate an integral kernel basis and a presentation for the cokernel.
    """
    a=[list(row) for row in matrix]
    rows=len(a)
    if any(len(row)!=cols for row in a):
        raise ValueError('matrix width mismatch')
    v=[[int(i==j) for j in range(cols)] for i in range(cols)]
    u=[[int(i==j) for j in range(rows)] for i in range(rows)] if track_rows else None
    def swap_rows(i,j):
        a[i],a[j]=a[j],a[i]
        if u is not None: u[i],u[j]=u[j],u[i]
    def swap_cols(i,j):
        for m in (a,v):
            for row in m: row[i],row[j]=row[j],row[i]
    def row_sub(i,j,factor):
        for t in range(cols): a[i][t]-=factor*a[j][t]
        if u is not None:
            for t in range(rows): u[i][t]-=factor*u[j][t]
    def col_sub(i,j,factor):
        for m in (a,v):
            for row in m: row[i]-=factor*row[j]
    k=0
    while k<min(rows,cols):
        entries=[(abs(a[i][j]),i,j) for i in range(k,rows) for j in range(k,cols) if a[i][j]]
        if not entries: break
        _,i,j=min(entries)
        swap_rows(k,i); swap_cols(k,j)
        while True:
            altered=False
            for i in range(k+1,rows):
                if a[i][k]:
                    quotient=a[i][k]//a[k][k]
                    row_sub(i,k,quotient)
                    if a[i][k]:
                        swap_rows(i,k)
                        altered=True
                        break
            if altered: continue
            for j in range(k+1,cols):
                if a[k][j]:
                    quotient=a[k][j]//a[k][k]
                    col_sub(j,k,quotient)
                    if a[k][j]:
                        swap_cols(k,j)
                        altered=True
                        break
            if altered: continue
            if any(a[i][k] for i in range(k+1,rows)) or any(a[k][j] for j in range(k+1,cols)):
                continue
            break
        k+=1
    if any(a[i][j] for i in range(rows) for j in range(cols) if i!=j):
        raise ArithmeticError('integer reduction did not diagonalize')
    return (a,v,k,u) if track_rows else (a,v,k)


def _integer_inverse(unimodular):
    n=len(unimodular)
    a=[[Fraction(z) for z in row]+[Fraction(int(i==j)) for j in range(n)]
       for i,row in enumerate(unimodular)]
    for k in range(n):
        pivot=next((i for i in range(k,n) if a[i][k]),None)
        if pivot is None: raise ArithmeticError('non-unimodular basis transformation')
        a[k],a[pivot]=a[pivot],a[k]
        t=a[k][k]
        a[k]=[x/t for x in a[k]]
        for i in range(n):
            if i!=k and a[i][k]:
                s=a[i][k]; a[i]=[x-s*y for x,y in zip(a[i],a[k])]
    result=[row[n:] for row in a]
    if any(x.denominator!=1 for row in result for x in row):
        raise ArithmeticError('basis transformation must be unimodular')
    return [[int(x) for x in row] for row in result]


def _invariant_factors(diagonal):
    """Convert an arbitrary direct sum of cyclics to invariant factors."""
    parts=sorted(abs(x) for x in diagonal if abs(x)>1)
    while True:
        altered=False
        for i in range(len(parts)):
            for j in range(i+1,len(parts)):
                if parts[j]%parts[i]:
                    a,b=parts[i],parts[j]
                    g=gcd(a,b)
                    parts[i],parts[j]=g,(a//g)*b
                    parts=sorted(x for x in parts if x>1)
                    altered=True
                    break
            if altered: break
        if not altered: return tuple(parts)


class IntegralChainComplex:
    """Finite free integral chain complex specified by generators and d matrices.

    ``ranks[n]`` is the number of C_n generators. ``boundaries[n]`` is a
    ROW-MAJOR integer matrix d_n:C_n -> C_(n-1). Missing maps are zero.
    """
    def __init__(self,ranks,boundaries):
        if not isinstance(ranks,(list,tuple)) or any(type(n) is not int or n<0 for n in ranks):
            raise ValueError('ranks must be nonnegative integers')
        self.ranks=tuple(ranks)
        if not isinstance(boundaries,dict):
            raise ValueError('boundaries must be a dictionary indexed by degree')
        if any(type(n) is not int or n<1 or n>=len(ranks) for n in boundaries):
            raise ValueError('boundary degree outside chain complex')
        self.boundaries={}
        for n in range(1,len(ranks)):
            matrix=boundaries.get(n,[[0]*ranks[n] for _ in range(ranks[n-1])])
            if len(matrix)!=ranks[n-1] or any(len(row)!=ranks[n] for row in matrix):
                raise ValueError('boundary matrix dimensions do not match ranks')
            if any(type(x) is not int for row in matrix for x in row):
                raise ValueError('boundary coefficients must be integers')
            self.boundaries[n]=[list(row) for row in matrix]
        for n in range(2,len(ranks)):
            a,b=self.boundaries[n-1],self.boundaries[n]
            if any(sum(a[i][k]*b[k][j] for k in range(ranks[n-1]))
                   for i in range(ranks[n-2]) for j in range(ranks[n])):
                raise ValueError(f'd_{n-1} d_{n} != 0')

    def homology(self):
        """Degree -> {'free_rank': int, 'torsion': tuple(invariant_factors)}."""
        answer={}
        for n,width in enumerate(self.ranks):
            a=self.boundaries.get(n,[])
            _,v,rnk=_diagonalize(a,width)
            inverse=_integer_inverse(v)
            right=self.boundaries.get(n+1,[[0]*self.ranks[n+1] for _ in range(width)] if n+1<len(self.ranks) else [[] for _ in range(width)])
            # Change the C_n basis by V, then restrict d_(n+1) to ker d_n.
            b=[[sum(inverse[i][k]*right[k][j] for k in range(width))
                for j in range(len(right[0]) if right else 0)] for i in range(width)]
            if any(x for row in b[:rnk] for x in row):
                raise ArithmeticError('boundary image does not lie in integer kernel')
            dim=width-rnk
            diagonal,_,incoming_rank=_diagonalize(b[rnk:],len(right[0]) if right else 0)
            torsion=_invariant_factors([diagonal[i][i] for i in range(incoming_rank)])
            answer[n]={'free_rank':dim-incoming_rank,'torsion':torsion}
        return answer

    def _class_input(self, degree, chain):
        if type(degree) is not int or not 0 <= degree < len(self.ranks):
            raise ValueError("degree outside chain complex")
        if not _integer_vector(chain, self.ranks[degree]):
            raise ValueError("chain must be an integer coordinate vector")
        preimage = json.dumps(
            [self.ranks, self.boundaries, degree, list(chain)],
            sort_keys=True, separators=(",", ":"),
        )
        width = self.ranks[degree]
        incoming_rank = self.ranks[degree + 1] if degree + 1 < len(self.ranks) else 0
        incoming = self.boundaries.get(degree + 1, [[] for _ in range(width)])
        outgoing_value = [_dot(row, chain) for row in self.boundaries.get(degree, [])]
        return sha256(preimage.encode("ascii")).hexdigest(), incoming, incoming_rank, outgoing_value

    def class_certificate(self, degree, chain):
        """Certify that a chain is a noncycle, boundary, or class of exact order."""
        digest, incoming, width, outgoing = self._class_input(degree, chain)
        base = {"kind": "integral-homology-class/1", "input_sha256": digest}
        for i, value in enumerate(outgoing):
            if value:
                return {**base, "status": "noncycle", "row": i, "value": value}
        if not any(chain):
            return {**base, "status": "boundary", "lift": [0] * width}
        diagonal, right, rank, left = _diagonalize(incoming, width, track_rows=True)
        values = [_dot(row, chain) for row in left]
        for i in range(rank, len(chain)):
            if values[i]:
                return {**base, "status": "infinite", "functional": left[i]}
        order = 1
        witnesses = []
        for i in range(rank):
            modulus = abs(diagonal[i][i])
            denominator = modulus // gcd(modulus, values[i])
            if denominator > 1:
                order = lcm(order, denominator)
                witnesses.append({"modulus": modulus, "functional": left[i]})
        solution = [order * values[i] // diagonal[i][i] if i < rank else 0
                    for i in range(width)]
        lift = [_dot(row, solution) for row in right]
        if order == 1:
            return {**base, "status": "boundary", "lift": lift}
        return {**base, "status": "torsion", "order": order, "lift": lift, "witnesses": witnesses}

    def verify_class_certificate(self, degree, chain, receipt):
        """Replay a receipt using matrix products and integer arithmetic only."""
        digest, incoming, width, outgoing = self._class_input(degree, chain)
        if (not isinstance(receipt, dict) or receipt.get("kind") != "integral-homology-class/1"
                or receipt.get("input_sha256") != digest):
            return False

        def fields(*extra):
            return set(receipt) == {"kind", "input_sha256", "status", *extra}

        status = receipt.get("status")
        if status == "noncycle":
            row, value = receipt.get("row"), receipt.get("value")
            return (fields("row", "value") and type(row) is int and 0 <= row < len(outgoing)
                    and type(value) is int and value != 0 and outgoing[row] == value)
        if any(outgoing):
            return False
        columns = list(zip(*incoming))

        def annihilates(functional, modulus=None):
            if not _integer_vector(functional, len(chain)):
                return False
            evaluations = [_dot(functional, col) for col in columns]
            return not any(x if modulus is None else x % modulus for x in evaluations)

        lift = receipt.get("lift")
        if status == "infinite":
            functional = receipt.get("functional")
            return fields("functional") and annihilates(functional) and _dot(functional, chain) != 0
        if not _integer_vector(lift, width):
            return False
        projected = [_dot(row, lift) for row in incoming]
        if status == "boundary":
            return fields("lift") and projected == list(chain)
        order, witnesses = receipt.get("order"), receipt.get("witnesses")
        if (status != "torsion" or not fields("order", "lift", "witnesses")
                or type(order) is not int or order < 2 or not isinstance(witnesses, list)
                or not witnesses or projected != [order * z for z in chain]):
            return False
        lower = 1
        for witness in witnesses:
            if not isinstance(witness, dict) or set(witness) != {"modulus", "functional"}:
                return False
            modulus, functional = witness["modulus"], witness["functional"]
            if type(modulus) is not int or modulus < 2 or not annihilates(functional, modulus):
                return False
            lower = lcm(lower, modulus // gcd(modulus, _dot(functional, chain)))
        return lower == order


def _dot(left, right):
    return sum(a * b for a, b in zip(left, right, strict=True))


def _integer_vector(vector, width):
    return (isinstance(vector, (list, tuple)) and len(vector) == width
            and all(type(x) is int for x in vector))


def integral_mapping_cone(source, target, maps):
    """Cone(f)_n = target_n ⊕ source_(n−1), ordered with target first.

    maps[n] is the row-major matrix of the integral chain map f_n. Omitted
    entries are zero. The cone validates the chain-map equations via d² = 0.
    """
    def rank(chain, degree):
        return chain.ranks[degree] if 0 <= degree < len(chain.ranks) else 0

    def rows(matrix, height, width):
        return [[0] * width for _ in range(height)] if matrix is None else matrix

    top = max(len(target.ranks), len(source.ranks) + bool(source.ranks))
    map_degrees = max(len(target.ranks), len(source.ranks))
    if not isinstance(maps, dict) or any(
        type(n) is not int or not 0 <= n < map_degrees for n in maps
    ):
        raise ValueError("invalid map degrees")
    for n, block in maps.items():
        if (not isinstance(block, (list, tuple)) or len(block) != rank(target, n)
                or any(not _integer_vector(row, rank(source, n)) for row in block)):
            raise ValueError("invalid integral chain-map block")
    ranks = [rank(target, n) + rank(source, n - 1) for n in range(top)]
    boundaries = {}
    for n in range(1, top):
        a, b = rank(target, n), rank(source, n - 1)
        c, d = rank(target, n - 1), rank(source, n - 2)
        target_boundary = rows(target.boundaries.get(n), c, a)
        source_boundary = rows(source.boundaries.get(n - 1), d, b)
        image = rows(maps.get(n - 1), c, b)
        boundaries[n] = ([list(left) + list(right) for left, right in zip(target_boundary, image, strict=True)]
                         + [[0] * a + [-x for x in row] for row in source_boundary])
    return IntegralChainComplex(ranks, boundaries)
