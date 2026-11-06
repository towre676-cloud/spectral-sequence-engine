"""Exact homological spectral sequences of finite filtered chain complexes.

The increasing-filtration convention is used:
  E^1[p,q] = H_{p+q}(F_p C / F_{p-1} C),
  d_r: E^r[p,q] -> E^r[p-r,q+r-1].

Everything is finite-dimensional linear algebra over a *prime* field.  This
module does not claim integral torsion, products, general convergence,
functoriality, or resolution of extension problems.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping


def _prime(p: int) -> bool:
    if isinstance(p, bool) or not isinstance(p, int) or p < 2:
        return False
    if p % 2 == 0:
        return p == 2
    d = 3
    while d * d <= p:
        if p % d == 0:
            return False
        d += 2
    return True


def _add_scaled(a, b, factor, prime):
    return tuple((x + factor * y) % prime for x, y in zip(a, b))


def _linear_combination(vectors, weights, prime, width):
    total = (0,) * width
    for vector, weight in zip(vectors, weights):
        if weight % prime:
            total = _add_scaled(total, vector, weight, prime)
    return total


class _Span:
    """Echelon span, retaining exact coordinates in *inserted* vectors."""

    def __init__(self, width: int, prime: int):
        self.width = width
        self.prime = prime
        self.vectors = []
        self.pivots = {}

    def add(self, vector):
        vector = tuple(x % self.prime for x in vector)
        if len(vector) != self.width:
            raise ValueError('vector width mismatch')
        residue = vector
        coefficients = [0] * (len(self.vectors) + 1)
        coefficients[-1] = 1
        while any(residue):
            pivot = max(i for i, x in enumerate(residue) if x)
            if pivot not in self.pivots:
                inv = pow(residue[pivot], -1, self.prime)
                self.pivots[pivot] = (
                    tuple((x * inv) % self.prime for x in residue),
                    tuple((x * inv) % self.prime for x in coefficients),
                )
                self.vectors.append(vector)
                return True
            reduced, witness = self.pivots[pivot]
            factor = (-residue[pivot]) % self.prime
            residue = _add_scaled(residue, reduced, factor, self.prime)
            coefficients = list(_add_scaled(coefficients, witness + (0,) * (len(coefficients) - len(witness)), factor, self.prime))
        return False

    def coordinates(self, vector):
        """Represent vector in the inserted basis; raise on nonmembership."""
        residue = tuple(x % self.prime for x in vector)
        if len(residue) != self.width:
            raise ValueError('vector width mismatch')
        coeff = (0,) * len(self.vectors)
        while any(residue):
            pivot = max(i for i, x in enumerate(residue) if x)
            if pivot not in self.pivots:
                raise ArithmeticError('vector is outside the claimed subspace')
            reduced, witness = self.pivots[pivot]
            factor = residue[pivot]
            residue = _add_scaled(residue, reduced, -factor, self.prime)
            coeff = _add_scaled(coeff, witness + (0,) * (len(coeff) - len(witness)), factor, self.prime)
        return coeff


def _kernel(columns, prime, output_width):
    """Basis of relations between columns of a linear transformation."""
    pivots = {}
    relations = []
    size = len(columns)
    for j, column in enumerate(columns):
        residue = tuple(x % prime for x in column)
        if len(residue) != output_width:
            raise ValueError('matrix height mismatch')
        combination = [int(i == j) for i in range(size)]
        while any(residue):
            pivot = max(i for i, x in enumerate(residue) if x)
            if pivot not in pivots:
                inv = pow(residue[pivot], -1, prime)
                pivots[pivot] = (
                    tuple((x * inv) % prime for x in residue),
                    tuple((x * inv) % prime for x in combination),
                )
                break
            reduced, witness = pivots[pivot]
            factor = -residue[pivot]
            residue = _add_scaled(residue, reduced, factor, prime)
            combination = list(_add_scaled(combination, witness, factor, prime))
        else:
            relations.append(tuple(combination))
    return tuple(relations)


def _rank(columns, prime, width):
    basis = _Span(width, prime)
    for column in columns:
        basis.add(column)
    return len(basis.vectors)


@dataclass(frozen=True)
class _PageSpace:
    cycles: tuple
    boundaries: tuple
    representatives: tuple


class FilteredChainComplex:
    """Generators (name, degree, birth, boundary coefficient dictionary).

    Every name is unique and strings are used for portable JSON certificates.
    Differential coefficients are reduced modulo the prime. The boundary of
    a degree-zero generator is necessarily zero.
    """

    def __init__(self, generators, prime=2):
        if not _prime(prime):
            raise ValueError('coefficient characteristic must be prime')
        self.prime = prime
        parsed = {}
        for record in generators:
            if len(record) != 4:
                raise ValueError('each generator must have four fields')
            name, degree, birth, boundary = record
            if not isinstance(name, str) or not name or name in parsed:
                raise ValueError('generator names must be unique, nonempty strings')
            if isinstance(degree, bool) or not isinstance(degree, int) or degree < 0:
                raise ValueError('degrees must be nonnegative integers')
            if isinstance(birth, bool) or not isinstance(birth, int) or birth < 0:
                raise ValueError('filtration levels must be nonnegative integers')
            if not isinstance(boundary, Mapping):
                raise ValueError('boundary must be a coefficient dictionary')
            if not all(isinstance(k, str) and isinstance(v, int) and not isinstance(v, bool)
                       for k, v in boundary.items()):
                raise ValueError('boundary keys and coefficients must be string and integer')
            parsed[name] = (degree, birth, {k: v % prime for k, v in boundary.items()
                                                  if v % prime})
        self.names = tuple(sorted(parsed, key=lambda k: (parsed[k][1], parsed[k][0], k)))
        self.data = parsed
        self.grades = tuple(sorted({birth for _, birth, _ in parsed.values()}))
        self.max_degree = max((degree for degree, _, _ in parsed.values()), default=-1)
        self.max_birth = max(self.grades, default=0)
        self.min_birth = min(self.grades, default=0)
        self.spread = self.max_birth - self.min_birth
        self.by_degree = {n: tuple(k for k in self.names if parsed[k][0] == n)
                          for n in range(self.max_degree + 1)}
        self.loc = {n: {name: i for i, name in enumerate(names)}
                    for n, names in self.by_degree.items()}
        for name, (degree, birth, boundary) in parsed.items():
            for target in boundary:
                if target not in parsed:
                    raise ValueError(f'{name}: unknown boundary generator {target}')
                tdegree, tbirth, _ = parsed[target]
                if tdegree != degree - 1:
                    raise ValueError(f'{name}: differential changes degree by more than one')
                if tbirth > birth:
                    raise ValueError(f'{name}: boundary raises filtration')
        self.columns = {}
        for n, names in self.by_degree.items():
            target = self.loc.get(n-1, {})
            self.columns[n] = tuple(tuple(parsed[name][2].get(t, 0)
                                          for t in self.by_degree.get(n-1, ()))
                                    for name in names)
        for n, names in self.by_degree.items():
            if n < 2:
                continue
            for name in names:
                if any(self.boundary(n-1, self.boundary(n, self.unit(n, name)))):
                    raise ValueError(f'{name}: d squared is not zero')
        self._cache = {}

    def unit(self, degree, name):
        coordinates = [0] * len(self.by_degree.get(degree, ()))
        coordinates[self.loc[degree][name]] = 1
        return tuple(coordinates)

    def boundary(self, degree, vector):
        if len(vector) != len(self.by_degree.get(degree, ())):
            raise ValueError('chain has incorrect degree or dimension')
        return _linear_combination(self.columns.get(degree, ()), vector,
                                   self.prime, len(self.by_degree.get(degree-1, ())))

    def _width(self, n):
        return len(self.by_degree.get(n, ()))

    def _lift(self, n, indices, weights):
        v = [0] * self._width(n)
        for i, value in zip(indices, weights):
            v[i] = value % self.prime
        return tuple(v)

    def _indices(self, n, predicate):
        return [j for j, name in enumerate(self.by_degree.get(n, ()))
                if predicate(self.data[name][1])]

    def _cycles(self, r, p, n):
        """Z^r_{p,n} = {x in F_p C_n : dx in F_{p-r} C_{n-1}}."""
        indices = self._indices(n, lambda birth: birth <= p)
        high = self._indices(n-1, lambda birth: birth > p-r)
        cols = [tuple(self.columns[n][i][j] for j in high) for i in indices]
        return tuple(self._lift(n, indices, combination)
                     for combination in _kernel(cols, self.prime, len(high)))

    def _boundaries(self, r, p, n):
        """F_p C_n intersect d(F_{p+r} C_{n+1})."""
        source = self._indices(n+1, lambda birth: birth <= p+r)
        high = self._indices(n, lambda birth: birth > p)
        cols = [self.columns.get(n+1, ())[i] for i in source]
        projected = [tuple(col[j] for j in high) for col in cols]
        vectors = [_linear_combination(cols, relation, self.prime, self._width(n))
                   for relation in _kernel(projected, self.prime, len(high))]
        basis = _Span(self._width(n), self.prime)
        for v in vectors:
            basis.add(v)
        return tuple(basis.vectors)

    def _space(self, r, p, n):
        if r < 1:
            raise ValueError('E_r requires r >= 1')
        key = (r, p, n)
        if key not in self._cache:
            numerator = self._cycles(r, p, n)
            denominator = (*self._cycles(r-1, p-1, n),
                           *self._boundaries(r-1, p, n))
            numerator_span = _Span(self._width(n), self.prime)
            for v in numerator:
                numerator_span.add(v)
            span = _Span(self._width(n), self.prime)
            for v in denominator:
                numerator_span.coordinates(v)  # the quotient must be well-defined
                span.add(v)
            boundaries = tuple(span.vectors)
            representatives = []
            for v in numerator:
                if span.add(v):
                    representatives.append(v)
            self._cache[key] = _PageSpace(tuple(numerator_span.vectors), boundaries,
                                          tuple(representatives))
        return self._cache[key]

    def page(self, r):
        """Nonzero (filtration p, total degree minus p) dimensions."""
        if not isinstance(r, int) or r < 0:
            raise ValueError('page index must be nonnegative')
        result = {}
        for p in self.grades:
            for n in range(self.max_degree + 1):
                if r == 0:
                    dimension = sum(self.data[k][1] == p for k in self.by_degree[n])
                else:
                    dimension = len(self._space(r, p, n).representatives)
                if dimension:
                    result[(p, n-p)] = dimension
        return result

    def differential(self, r):
        """Exact d_r matrices, as columns in the chosen quotient bases.

        Keys are nonempty source positions (p,q). Each value contains target
        position, coefficient columns (including zero columns) and rank.
        """
        if not isinstance(r, int) or r < 0:
            raise ValueError('differential index must be nonnegative')
        result = {}
        for (p, q), dimension in self.page(r).items():
            n = p + q
            if r == 0:
                sources = [self.unit(n, name) for name in self.by_degree[n]
                           if self.data[name][1] == p]
                target_names = [name for name in self.by_degree.get(n-1, ())
                                if self.data[name][1] == p]
                target_indices = [self.loc[n-1][name] for name in target_names]
                columns = [tuple(self.boundary(n, v)[j] for j in target_indices)
                           for v in sources]
                height = len(target_names)
            else:
                sources = self._space(r, p, n).representatives
                target = self._space(r, p-r, n-1)
                basis = _Span(self._width(n-1), self.prime)
                for vec in (*target.boundaries, *target.representatives):
                    basis.add(vec)
                d = len(target.boundaries)
                columns = [basis.coordinates(self.boundary(n, v))[d:]
                           for v in sources]
                height = len(target.representatives)
            if len(columns) != dimension:
                raise ArithmeticError('incorrect source dimension')
            result[(p, q)] = {'target': (p-r, q+r-1),
                              'columns': tuple(columns),
                              'rank': _rank(columns, self.prime, height)}
        return result

    def betti(self):
        result = {}
        for n in range(self.max_degree + 1):
            d_in = _rank(self.columns.get(n+1, ()), self.prime, self._width(n))
            d_out = _rank(self.columns[n], self.prime, self._width(n-1))
            b = self._width(n) - d_out - d_in
            if b < 0:
                raise ArithmeticError('invalid homology dimension')
            if b:
                result[n] = b
        return result

    def intervals(self, include_zero=False):
        """Independent ordinary persistence reduction in filtration order."""
        pos = {name: i for i, name in enumerate(self.names)}
        reduced = {}
        pairs = {}
        creators = []
        for j, name in enumerate(self.names):
            degree, _, boundary = self.data[name]
            col = {pos[k]: coefficient for k, coefficient in boundary.items()}
            while col:
                low = max(col)
                if low not in reduced:
                    inv = pow(col[low], -1, self.prime)
                    reduced[low] = {i: (v*inv) % self.prime for i, v in col.items()}
                    pairs[low] = j
                    break
                factor = -col[low]
                for i, v in reduced[low].items():
                    value = (col.get(i, 0) + factor*v) % self.prime
                    if value:
                        col[i] = value
                    else:
                        col.pop(i, None)
            if not col:
                creators.append(j)
        result = []
        for i in creators:
            name = self.names[i]
            degree, birth, _ = self.data[name]
            death = self.data[self.names[pairs[i]]][1] if i in pairs else None
            if include_zero or death is None or birth < death:
                result.append((degree, birth, death))
        return sorted(result, key=lambda a: (a[0], a[1], a[2] is None,
                                              -1 if a[2] is None else a[2]))

    def verify(self):
        """Recompute the exact page identity, square-zero and convergence."""
        pages = [self.page(r) for r in range(self.spread + 3)]
        for r in range(self.spread + 2):
            d = self.differential(r)
            current = pages[r]
            for pos, dimension in current.items():
                entry = d[pos]
                next_entry = d.get(entry['target'])
                for col in entry['columns']:
                    if next_entry is not None:
                        composed = _linear_combination(next_entry['columns'], col,
                                                     self.prime, len(next_entry['columns'][0])
                                                     if next_entry['columns'] else 0)
                        if any(composed):
                            raise ArithmeticError(f'd_{r} squared is nonzero at {pos}')
                inbound = d.get((pos[0]+r, pos[1]-r+1))
                predicted = dimension - entry['rank'] - (inbound['rank'] if inbound else 0)
                if predicted != pages[r+1].get(pos, 0):
                    raise ArithmeticError(f'E_{r+1} is not the homology of E_{r} at {pos}')
            if any(pos not in current for pos in pages[r+1]):
                raise ArithmeticError('new generator appeared on a later page')
        if pages[-1] != pages[-2]:
            raise ArithmeticError('the spectral sequence did not stabilize')
        infinity = defaultdict(int)
        for (p, q), b in pages[self.spread + 1].items():
            infinity[p+q] += b
        if dict(infinity) != self.betti():
            raise ArithmeticError('associated graded does not match total homology')
        bars = Counter(death-birth for _, birth, death in self.intervals(True)
                       if death is not None)
        for r in range(self.spread+1):
            actual = sum(v['rank'] for v in self.differential(r).values())
            if actual != bars[r]:
                raise ArithmeticError(f'd_{r} rank does not agree with barcode lengths')
        if sum(1 for _, _, death in self.intervals(True) if death is None) != sum(self.betti().values()):
            raise ArithmeticError('essential barcodes do not recover total Betti')
        return True

    def to_dict(self):
        return {'prime': self.prime, 'generators': [
            [name, self.data[name][0], self.data[name][1], self.data[name][2]]
            for name in self.names]}

    def certificate(self):
        """Deterministic replay report. Not a proof-assistant certificate."""
        self.verify()
        source = self.to_dict()
        packed = lambda v: json.dumps(v, sort_keys=True, separators=(',', ':'))
        pages = []
        for r in range(self.spread+2):
            page = self.page(r)
            differentials = self.differential(r)
            pages.append({'r': r,
                          'groups': [{'p': p, 'q': q, 'dimension': b}
                                     for (p, q), b in sorted(page.items())],
                          'differentials': [
                              {'source': [p, q], 'target': list(d['target']),
                               'columns': [list(c) for c in d['columns']], 'rank': d['rank']}
                              for (p,q), d in sorted(differentials.items())]})
        body = {'schema': 'spectral-sequence-replay-v1',
                'convention': 'increasing-filtration-homological',
                'input': source,
                'input_sha256': sha256(packed(source).encode()).hexdigest(),
                'pages': pages,
                'stable_at': self.spread+1,
                'betti': {str(k): v for k, v in self.betti().items()},
                'intervals': [list(x) for x in self.intervals(True)],
                'validated': True}
        body['report_sha256'] = sha256(packed(body).encode()).hexdigest()
        return body

    @classmethod
    def from_simplicial(cls, simplices, prime=2):
        """Oriented boundary with signs (-1)^j, for arbitrary prime fields."""
        ordered = sorted(((tuple(sorted(vertices)), birth) for vertices, birth in simplices),
                         key=lambda x: (x[1], len(x[0]), x[0]))
        names = {simplex: f's{i}' for i, (simplex, _) in enumerate(ordered)}
        if len(names) != len(ordered):
            raise ValueError('duplicate simplex')
        records = []
        for simplex, birth in ordered:
            if not simplex or len(set(simplex)) != len(simplex):
                raise ValueError('invalid simplex')
            faces = [simplex[:j] + simplex[j+1:] for j in range(len(simplex))]
            boundary = {}
            if len(simplex) > 1:
                for j, face in enumerate(faces):
                    if face not in names:
                        raise ValueError('missing face')
                    boundary[names[face]] = (-1)**j
            records.append((names[simplex], len(simplex)-1, birth, boundary))
        return cls(records, prime)


def replay_certificate(receipt):
    """Check digest and independently recompute all reported data from input."""
    if receipt.get('schema') != 'spectral-sequence-replay-v1':
        raise ValueError('unsupported certificate schema')
    body = {k: v for k, v in receipt.items() if k != 'report_sha256'}
    packed = json.dumps(body, sort_keys=True, separators=(',', ':'))
    if sha256(packed.encode()).hexdigest() != receipt.get('report_sha256'):
        raise ValueError('certificate digest mismatch')
    actual = FilteredChainComplex(**receipt['input']).certificate()
    if actual != receipt:
        raise ValueError('certificate computation mismatch')
    return True


def _example(which):
    if which == 'triangle':
        return FilteredChainComplex.from_simplicial([
            ((0,),0), ((1,),0), ((2,),0),
            ((0,1),1), ((0,2),1), ((1,2),1), ((0,1,2),2)])
    if which in ('d2', 'd3'):
        r = int(which[1:])
        return FilteredChainComplex([
            ('a',0,0,{}), ('b',1,r,{'a':1})])
    raise ValueError('unknown example')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', nargs='?', help='JSON file containing prime and generators')
    parser.add_argument('--example', choices=('triangle','d2','d3'))
    parser.add_argument('--certificate', help='write the deterministic JSON replay certificate')
    parser.add_argument('--replay', help='verify an existing certificate instead')
    args = parser.parse_args(argv)
    if args.replay:
        replay_certificate(json.loads(Path(args.replay).read_text(encoding='utf-8')))
        print('Certificate replay: PASS')
        return
    if (args.input is None) == (args.example is None):
        parser.error('select exactly one of input JSON or --example')
    engine = (_example(args.example) if args.example else
              FilteredChainComplex(**json.loads(Path(args.input).read_text(encoding='utf-8'))))
    receipt = engine.certificate()
    print(f"Prime field F_{engine.prime}; {len(engine.names)} generators")
    for item in receipt['pages']:
        groups = ', '.join(f"({x['p']},{x['q']}):{x['dimension']}" for x in item['groups'])
        diff = sum(x['rank'] for x in item['differentials'])
        print(f"E_{item['r']}: {groups or '(zero)'}; rank(d_{item['r']})={diff}")
    print(f"E_infinity at E_{receipt['stable_at']}; Betti={engine.betti()}")
    print('Exact verification: PASS')
    if args.certificate:
        Path(args.certificate).write_text(json.dumps(receipt, indent=2)+'\n', encoding='utf-8')
        print(f'Certificate: {args.certificate}')


if __name__ == '__main__':
    main()
