"""Exact Bockstein ranks and explicit lifted differentials for free Z-chain complexes.

The mod-p Bockstein spectral sequence has E_1[n] = H_n(C; F_p),
 d_r: E_r[n] -> E_r[n-1], and E_infty[n] = (H_n(C; Z)/tors) ⊗ F_p.
For finite free C, a Z/p^r summand of H_n contributes one rank to
 d_r: E_r[n+1] -> E_r[n].

BocksteinSpectralSequence gives independent rank fingerprints from integral
homology. LiftedBockstein computes exact higher connecting differentials
using explicit integer lifts. Matrix entries depend on chosen quotient bases.
"""
from __future__ import annotations

from collections import Counter
from hashlib import sha256
import json

from integral import IntegralChainComplex
from spectral_sequence import _prime, _rank, _Span, _kernel, _linear_combination


def _valuation(m: int, prime: int) -> int:
    """Exponent of a prime in positive integer m."""
    if m < 1:
        raise ValueError('valuation requires a positive integer')
    count = 0
    while m % prime == 0:
        m //= prime
        count += 1
    return count


class BocksteinSpectralSequence:
    """Prime-primary Bockstein fingerprints of an integral chain complex.

    Results are *graded ranks*, not matrices or canonical generators. The
    page n-degree dimensions follow the higher universal coefficient exact
    sequences; differentials are characterized by their ranks.
    """
    def __init__(self, chain: IntegralChainComplex, prime: int):
        if not isinstance(chain, IntegralChainComplex):
            raise ValueError('chain must be an IntegralChainComplex')
        if not _prime(prime):
            raise ValueError('coefficient must be a prime field')
        self.chain = chain
        self.prime = prime
        self.integral_homology = chain.homology()
        self.exponents = {}
        for n, group in self.integral_homology.items():
            self.exponents[n] = Counter(
                exp for m in group['torsion']
                if (exp := _valuation(m, prime)) > 0)
        self.max_exponent = max((max(v, default=0) for v in self.exponents.values()), default=0)

    def _alive(self, n: int, r: int) -> int:
        return sum(multiplicity for exp, multiplicity in self.exponents.get(n, {}).items()
                   if exp >= r)

    def page(self, r: int) -> dict[int, int]:
        """n -> dim over F_p of E_r[n], r >= 1."""
        if type(r) is not int or r < 1:
            raise ValueError('Bockstein pages begin at r=1')
        return {n: size for n in range(len(self.chain.ranks))
                if (size := self.integral_homology[n]['free_rank']
                    + self._alive(n, r) + self._alive(n-1, r))}

    def differential_ranks(self, r: int) -> dict[int, int]:
        """source degree n -> rank of beta_r: E_r[n] -> E_r[n-1]."""
        if type(r) is not int or r < 1:
            raise ValueError('Bockstein differentials begin at r=1')
        return {n+1: count for n in range(len(self.chain.ranks))
                if (count := self.exponents[n].get(r, 0))}

    def torsion_from_differentials(self) -> dict[int, tuple[int, ...]]:
        """Reconstruct p-primary elementary divisors from beta ranks."""
        factors = {}
        for n in range(len(self.chain.ranks)):
            factors[n] = tuple(self.prime**r
                               for r in range(1, self.max_exponent+1)
                               for _ in range(self.differential_ranks(r).get(n+1, 0)))
        return factors

    def verify(self) -> bool:
        """Mod-p chain homology, UCT, page homology, and torsion checks."""
        p = self.prime
        chain = self.chain
        for n, width in enumerate(chain.ranks):
            outgoing = chain.boundaries.get(n, [])
            incoming = chain.boundaries.get(n+1, [])
            out_cols = [tuple(row[j] % p for row in outgoing)
                        for j in range(width)] if n else []
            in_cols = [tuple(row[j] % p for row in incoming)
                       for j in range(chain.ranks[n+1])] if n+1 < len(chain.ranks) else []
            computed = width - _rank(out_cols, p, chain.ranks[n-1] if n else 0) - _rank(in_cols, p, width)
            if computed != self.page(1).get(n, 0):
                raise ArithmeticError('UCT disagrees with direct mod-p homology')
        for r in range(1, self.max_exponent+2):
            page = self.page(r)
            ranks = self.differential_ranks(r)
            following = self.page(r+1)
            for n in range(len(chain.ranks)):
                expected = page.get(n, 0) - ranks.get(n, 0) - ranks.get(n+1, 0)
                if expected != following.get(n, 0):
                    raise ArithmeticError(f'Bockstein page identity failed at r={r}, n={n}')
        limit = self.page(self.max_exponent+1)
        free = {n: h['free_rank'] for n, h in self.integral_homology.items()
                if h['free_rank']}
        if limit != free:
            raise ArithmeticError('Bockstein E_infty does not recover free homology')
        for n, factors in self.torsion_from_differentials().items():
            original = tuple(sorted(p**v for exponent, count in self.exponents[n].items()
                                    for v in [exponent] for _ in range(count)))
            if tuple(sorted(factors)) != original:
                raise ArithmeticError('Bockstein torsion reconstruction failed')
        return True

    def certificate(self) -> dict:
        """Deterministic self-checking summary; replay is not independent proof."""
        self.verify()
        result = {
            'schema': 'bockstein-ranks-v1', 'prime': self.prime,
            'chain': {'ranks': list(self.chain.ranks),
                      'boundaries': {str(n): self.chain.boundaries[n]
                                     for n in sorted(self.chain.boundaries)}},
            'homology': {str(n): {'free_rank': h['free_rank'], 'torsion': list(h['torsion'])}
                         for n, h in self.integral_homology.items()},
            'pages': [
                {'r': r, 'dimensions': {str(n): x for n,x in self.page(r).items()},
                 'differential_ranks': {str(n): x for n,x in self.differential_ranks(r).items()}}
                for r in range(1, self.max_exponent+2)],
            'stabilizes_at': self.max_exponent+1,
            'nonclaim': 'No chain-level higher Bockstein maps or canonical representatives',
        }
        packed = json.dumps(result, sort_keys=True, separators=(',', ':'))
        result['sha256'] = sha256(packed.encode()).hexdigest()
        return result


class LiftedBockstein:
    """Compute *actual* higher Bockstein maps using integral chain lifts.

    For a chosen E_r class x (integer lift with d(x) divisible by p**r),
    beta_r[x] is represented by d(x)/p**r mod p. A beta_r cycle is lifted
    to the next page by subtracting integer boundary preimages and previous
    Bockstein image witnesses. The algorithm is finite, exact, and retains
    explicit representatives, but does not return canonical matrices.
    """
    def __init__(self, chain: IntegralChainComplex, prime: int):
        self.fingerprint=BocksteinSpectralSequence(chain,prime)
        self.chain=chain
        self.prime=prime
        self._pages=None

    def _width(self,n):
        return self.chain.ranks[n] if 0 <= n < len(self.chain.ranks) else 0

    def _boundary(self,degree,coefficients):
        if len(coefficients)!=self._width(degree):
            raise ArithmeticError('integral chain width mismatch')
        if degree<=0:
            return ()
        matrix=self.chain.boundaries[degree]
        return tuple(sum(a*x for a,x in zip(row,coefficients)) for row in matrix)

    def _initial_page(self):
        result={}
        p=self.prime
        for degree in range(len(self.chain.ranks)):
            width=self._width(degree)
            incoming=_Span(width,p)
            for j in range(self._width(degree+1)):
                unit=tuple(int(i==j) for i in range(self._width(degree+1)))
                incoming.add(tuple(v%p for v in self._boundary(degree+1,unit)))
            columns=[]
            for j in range(width):
                unit=tuple(int(i==j) for i in range(width))
                columns.append(tuple(v%p for v in self._boundary(degree,unit)))
            ker=_kernel(columns,p,self._width(degree-1))
            representatives=[]
            for vec in ker:
                if incoming.add(vec):
                    representatives.append(vec)
            result[degree]=tuple(representatives)
        return result

    def _ambient(self,degree,previous,representatives=()):
        """Ambient mod-p span: ordinary boundaries, earlier betas, E_r reps."""
        p=self.prime
        span=_Span(self._width(degree),p)
        corrections=[]
        for j in range(self._width(degree+1)):
            unit=tuple(int(i==j) for i in range(self._width(degree+1)))
            image=tuple(v%p for v in self._boundary(degree+1,unit))
            if span.add(image):
                corrections.append(('boundary',unit,0))
        for old_r,source in previous.get(degree,[]):
            actual=self._boundary(degree+1,source)
            divisor=p**old_r
            if any(value%divisor for value in actual):
                raise ArithmeticError('earlier Bockstein witness is not divisible')
            image=tuple((value//divisor)%p for value in actual)
            if span.add(image):
                corrections.append(('old',source,old_r))
        count=len(corrections)
        for vec in representatives:
            if not span.add(tuple(z%p for z in vec)):
                raise ArithmeticError('page representative is not independent modulo earlier images')
        return span,corrections,count

    def _quotient_image(self,degree,source,r,previous,targets):
        p=self.prime
        differential=self._boundary(degree,source)
        if any(x%(p**r) for x in differential):
            raise ArithmeticError('Bockstein input lacks required lift')
        image=tuple((x//(p**r))%p for x in differential)
        ambient,_,start=self._ambient(degree-1,previous,targets)
        return ambient.coordinates(image)[start:]

    def _lift_cycle(self,degree,weights,representatives,r,previous):
        """Correct beta_r cycle to an integer lift modulo p**(r+1)."""
        p=self.prime
        width=self._width(degree)
        chain=tuple(sum(a*v[i] for a,v in zip(weights,representatives))
                    for i in range(width))
        actual=self._boundary(degree,chain)
        divisor=p**r
        if any(v%divisor for v in actual):
            raise ArithmeticError('not an E_r cycle')
        image=tuple((v//divisor)%p for v in actual)
        span,corrections,_=self._ambient(degree-1,previous)
        coeff=span.coordinates(image)
        change=[0]*width
        for multiplier,(kind,source,old_r) in zip(coeff,corrections):
            factor=multiplier*(divisor if kind=='boundary' else p**(r-old_r))
            for i,v in enumerate(source):
                change[i]+=factor*v
        new=tuple(a-b for a,b in zip(chain,change))
        if any(v%(p**(r+1)) for v in self._boundary(degree,new)):
            raise ArithmeticError('integral lift correction failed')
        if tuple(x%p for x in new)!=tuple(x%p for x in chain):
            raise ArithmeticError('integral lift changed its mod-p homology class')
        return new

    def _compute(self):
        if self._pages is not None:
            return
        p=self.prime
        bound=self.fingerprint.max_exponent+1
        reps=self._initial_page()
        previous={}
        pages=[]
        for r in range(1,bound+1):
            matrices={}
            for degree in range(len(self.chain.ranks)):
                columns=tuple(self._quotient_image(degree,chain,r,previous,
                                reps.get(degree-1,())) for chain in reps.get(degree,()))
                matrices[degree]=columns
            observed={n:len(v) for n,v in reps.items() if v}
            if observed!=self.fingerprint.page(r):
                raise ArithmeticError(f'lifted E_{r} dimension disagrees with UCT')
            for degree,cols in matrices.items():
                rank=_rank(cols,p,len(reps.get(degree-1,())))
                expected=self.fingerprint.differential_ranks(r).get(degree,0)
                if rank!=expected:
                    raise ArithmeticError(f'lifted beta_{r} rank disagrees with p-primary torsion')
                # Check beta_r squared directly in computed quotient coordinates.
                for col in cols:
                    composed=_linear_combination(matrices.get(degree-1,()),col,p,
                                                  len(reps.get(degree-2,())))
                    if any(composed):
                        raise ArithmeticError(f'beta_{r} squared is nonzero')
            pages.append((r,reps,matrices))
            if r==bound:
                break
            next_reps={}
            for degree in range(len(self.chain.ranks)):
                dimension=len(reps.get(degree,()))
                incoming=_Span(dimension,p)
                for v in matrices.get(degree+1,()):
                    incoming.add(v)
                cycles=_kernel(matrices[degree],p,len(reps.get(degree-1,())))
                produced=[]
                for relation in cycles:
                    if incoming.add(relation):
                        produced.append(self._lift_cycle(degree,relation,reps[degree],r,previous))
                next_reps[degree]=tuple(produced)
            for degree,chains in reps.items():
                if degree>0:
                    previous.setdefault(degree-1,[]).extend((r,x) for x in chains)
            reps=next_reps
        self._pages=pages

    def page(self,r):
        """n -> dimension of the actual lifted E_r, r >= 1."""
        self._compute()
        if type(r) is not int or r<1:
            raise ValueError('Bockstein pages begin at one')
        if r>len(self._pages):
            return self.fingerprint.page(r)
        return {n:len(x) for n,x in self._pages[r-1][1].items() if x}

    def differential(self,r):
        """n -> column-major beta_r matrix in explicit E_r basis lifts."""
        self._compute()
        if type(r) is not int or r<1:
            raise ValueError('Bockstein pages begin at one')
        if r>len(self._pages):
            return {n:() for n in range(len(self.chain.ranks))}
        return self._pages[r-1][2]

    def representatives(self,r):
        """n -> integer chains lifting every E_r basis class with d x=0 mod p**r."""
        self._compute()
        if type(r) is not int or not 1<=r<=len(self._pages):
            raise ValueError('request a computed nonzero page index')
        return self._pages[r-1][1]

    def verify(self):
        self._compute()
        if not self.fingerprint.verify():
            raise ArithmeticError('fingerprint check failed')
        for r,reps,matrices in self._pages:
            for degree,chains in reps.items():
                for chain in chains:
                    if any(v%(self.prime**r) for v in self._boundary(degree,chain)):
                        raise ArithmeticError('invalid stored integer lift')
                cols=matrices[degree]
                if len(cols)!=len(chains) or any(len(x)!=len(reps.get(degree-1,())) for x in cols):
                    raise ArithmeticError('incorrect beta matrix shape')
        return True

    def certificate(self):
        """Replay-oriented exact matrices and integer lifting witnesses."""
        self.verify()
        report={'schema':'bockstein-lift-v1','fingerprint':self.fingerprint.certificate(),
                'pages':[{'r':r,
                          'representatives':{str(n):[list(v) for v in chains]
                                             for n,chains in reps.items()},
                          'differentials':{str(n):[list(col) for col in mat]
                                           for n,mat in matrices.items()}}
                         for r,reps,matrices in self._pages]}
        report['sha256']=sha256(json.dumps(report,sort_keys=True,separators=(',',':')).encode()).hexdigest()
        return report


def replay_lifted_bockstein_certificate(receipt):
    """Validate hash and exactly recompute a lifted Bockstein receipt.

    Replays with the same Python implementation, not an independent checker.
    """
    if not isinstance(receipt,dict) or receipt.get('schema')!='bockstein-lift-v1':
        raise ValueError('unexpected Bockstein receipt schema')
    digest=receipt.get('sha256')
    body={k:v for k,v in receipt.items() if k!='sha256'}
    expected=sha256(json.dumps(body,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    if digest!=expected:
        raise ValueError('Bockstein certificate hash mismatch')
    fp=receipt.get('fingerprint')
    if not isinstance(fp,dict) or 'chain' not in fp or 'prime' not in fp:
        raise ValueError('incomplete Bockstein certificate')
    data=fp['chain']
    c=IntegralChainComplex(data['ranks'],{int(n):m for n,m in data['boundaries'].items()})
    regenerated=LiftedBockstein(c,fp['prime']).certificate()
    if regenerated!=receipt:
        raise ValueError('Bockstein certificate replay discrepancy')
    return True
