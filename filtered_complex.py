"""Finite simplicial filtrations over F_2.

Computes ordinary one-parameter persistence and the ranks of the E^1
associated-graded homology groups E^1_{p,q}=H_{p+q}(F_p,F_{p-1}).
Does NOT yet compute d^r, r>=1, extension data, or general sheaf cohomology.
"""
from collections import defaultdict


def _rank(columns):
    """Rank of a binary matrix represented by integer-bitset columns."""
    pivots = {}
    for original in columns:
        v = original
        while v:
            p = v.bit_length() - 1
            if p in pivots:
                v ^= pivots[p]
            else:
                pivots[p] = v
                break
    return len(pivots)


class FilteredComplex:
    def __init__(self, cells):
        """cells: iterable of (tuple of vertices, nonnegative integer birth)."""
        raw = {}
        for vertices, birth in cells:
            simplex = tuple(sorted(vertices))
            if not simplex or len(set(simplex)) != len(simplex):
                raise ValueError("cells must be nonempty simplices, without repeated vertices")
            if not isinstance(birth, int) or birth < 0:
                raise ValueError("filtration levels must be nonnegative integers")
            if simplex in raw:
                raise ValueError("duplicate simplex")
            raw[simplex] = birth
        for simplex, birth in raw.items():
            for j in range(len(simplex)):
                face = simplex[:j] + simplex[j+1:]
                if face and (face not in raw or raw[face] > birth):
                    raise ValueError("missing face or filtration order violation")
        self.cells = sorted(raw.items(), key=lambda pair: (pair[1], len(pair[0]), pair[0]))
        self.index = {cell: i for i, (cell, _) in enumerate(self.cells)}

    def _boundary(self, simplex):
        if len(simplex) == 1:
            return 0
        return sum(1 << self.index[simplex[:j] + simplex[j+1:]]
                   for j in range(len(simplex)))

    def intervals(self, include_zero=False):
        """(homology degree, birth, death|None) over F2; None means infinite."""
        reduced = {}
        pairs = {}
        creators = set()
        for j, (simplex, _) in enumerate(self.cells):
            column = self._boundary(simplex)
            while column:
                low = column.bit_length() - 1
                if low not in reduced:
                    reduced[low] = column
                    pairs[low] = j
                    break
                column ^= reduced[low]
            if column == 0:
                creators.add(j)
        result = []
        for i in sorted(creators):
            simplex, birth = self.cells[i]
            death = self.cells[pairs[i]][1] if i in pairs else None
            if include_zero or death is None or death > birth:
                result.append((len(simplex)-1, birth, death))
        return sorted(result, key=lambda x: (x[0], x[1], x[2] is None,
                                             -1 if x[2] is None else x[2]))

    def betti(self, level):
        """Betti numbers of the subcomplex born at or before level."""
        if not isinstance(level, int):
            raise ValueError("level must be integer")
        counts = defaultdict(int)
        for dim, birth, death in self.intervals():
            if birth <= level and (death is None or level < death):
                counts[dim] += 1
        return dict(sorted(counts.items()))

    def e1_ranks(self):
        """Nonzero entries (p,q): dim E1[p,q]. No higher-page claims."""
        grouped = defaultdict(list)
        for simplex, birth in self.cells:
            grouped[(birth, len(simplex)-1)].append(simplex)
        result = {}
        for p in sorted({birth for _, birth in self.cells}):
            dimensions = [d for (level, d) in grouped if level == p]
            if not dimensions:
                continue
            max_dim = max(dimensions)
            rank = {}
            for d in range(1, max_dim+1):
                lower = grouped.get((p, d-1), [])
                rows = {cell: i for i, cell in enumerate(lower)}
                columns = []
                for simplex in grouped.get((p, d), []):
                    column = 0
                    for j in range(len(simplex)):
                        face = simplex[:j] + simplex[j+1:]
                        if face in rows:
                            column ^= 1 << rows[face]
                    columns.append(column)
                rank[d] = _rank(columns)
            for d in range(max_dim+1):
                count = len(grouped.get((p, d), []))
                b = count - rank.get(d, 0) - rank.get(d+1, 0)
                if b < 0:
                    raise ArithmeticError("invalid chain complex rank")
                if b:
                    result[(p, d-p)] = b
        return result
