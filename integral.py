"""Integral homology of finite free chain complexes using unimodular reductions.

Total integral homology only: this does NOT compute integral spectral pages,
Bockstein differentials, hidden extensions, or coefficients in a general ring.
No external dependencies. Reduction is cubic-ish for small examples but integer
coefficient growth can be prohibitive on large complexes.
"""
from fractions import Fraction
from math import gcd


def _diagonalize(matrix, cols):
    """D -> U D V diagonal, tracking the integer-unimodular V.

    The diagonal need not be Smith-divisibility ordered. This suffices to
    calculate an integral kernel basis and a presentation for the cokernel.
    """
    a=[list(row) for row in matrix]
    rows=len(a)
    if any(len(row)!=cols for row in a):
        raise ValueError('matrix width mismatch')
    v=[[int(i==j) for j in range(cols)] for i in range(cols)]
    def swap_rows(i,j): a[i],a[j]=a[j],a[i]
    def swap_cols(i,j):
        for m in (a,v):
            for row in m: row[i],row[j]=row[j],row[i]
    def row_sub(i,j,factor):
        for t in range(cols): a[i][t]-=factor*a[j][t]
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
    return a,v,k


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
