"""Filtered chain maps and induced spectral-sequence morphisms over prime fields.

A map is certified only in the computational sense: its generators have the
required degree and filtration and it commutes with both chain differentials.
"""

from spectral_sequence import _Span, _linear_combination


class FilteredChainMap:
    """Filtered degree-zero chain map between two FilteredChainComplex objects.

    ``images`` maps source generator names to dictionaries of target names and
    integer coefficients. Omitted generators map to zero.
    """

    def __init__(self, source, target, images):
        if source.prime != target.prime:
            raise ValueError('source and target fields must agree')
        if not isinstance(images, dict) or set(images) - set(source.names):
            raise ValueError('image declarations must name source generators')
        self.source, self.target, self.prime = source, target, source.prime
        self.images = {}
        for name in source.names:
            image = images.get(name, {})
            if not isinstance(image, dict):
                raise ValueError('images must be coefficient dictionaries')
            degree, birth, _ = source.data[name]
            vector = [0] * target._width(degree)
            for t, c in image.items():
                if t not in target.data or not isinstance(c, int) or isinstance(c, bool):
                    raise ValueError('unknown image generator or invalid coefficient')
                target_degree, target_birth, _ = target.data[t]
                if target_degree != degree or (c % self.prime and target_birth > birth):
                    raise ValueError('map must preserve degree and filtration')
                if target_degree == degree:
                    i = target.loc[degree][t]
                    vector[i] = (vector[i] + c) % self.prime
            self.images[name] = tuple(vector)
        for name in source.names:
            degree = source.data[name][0]
            left = target.boundary(degree, self.images[name])
            right = self.apply(degree-1, source.boundary(degree, source.unit(degree, name)))
            if left != right:
                raise ValueError(f'chain-map square does not commute at {name}')

    def apply(self, degree, vector):
        if len(vector) != self.source._width(degree):
            raise ValueError('input has incorrect degree or dimension')
        return _linear_combination(
            [self.images[name] for name in self.source.by_degree.get(degree, ())],
            vector, self.prime, self.target._width(degree))

    def on_page(self, r):
        """Matrices of the induced map E_r(source) to E_r(target)."""
        result = {}
        for (p, q), dimension in self.source.page(r).items():
            n = p+q
            if r == 0:
                reps = [self.source.unit(n, k) for k in self.source.by_degree[n]
                        if self.source.data[k][1] == p]
                coords = [self.target.loc[n][k] for k in self.target.by_degree.get(n, ())
                          if self.target.data[k][1] == p]
                cols = [tuple(self.apply(n, v)[i] for i in coords) for v in reps]
            else:
                reps = self.source._space(r,p,n).representatives
                target = self.target._space(r,p,n)
                basis = _Span(self.target._width(n), self.prime)
                for v in (*target.boundaries, *target.representatives):
                    basis.add(v)
                cols = [basis.coordinates(self.apply(n, v))[len(target.boundaries):]
                        for v in reps]
            if len(cols) != dimension:
                raise ArithmeticError('incorrect induced map shape')
            result[(p,q)] = tuple(cols)
        return result

    def verify_naturality(self):
        """Check d_r f_r = f_r d_r, including every stable page."""
        last = max(self.source.spread, self.target.spread) + 1
        for r in range(last+1):
            fs = self.on_page(r)
            ds = self.source.differential(r)
            dt = self.target.differential(r)
            for position, columns in fs.items():
                target_pos = ds[position]['target']
                left_map = dt.get(position)
                right_map = fs.get(target_pos, ())
                for j, image in enumerate(columns):
                    out_width = len(self.target._space(r, target_pos[0], target_pos[0]+target_pos[1]).representatives) if r else self.target.page(0).get(target_pos,0)
                    left = (_linear_combination(left_map['columns'], image, self.prime, out_width)
                            if left_map else (0,)*out_width)
                    right = _linear_combination(right_map, ds[position]['columns'][j],
                                                self.prime, out_width)
                    if left != right:
                        raise ArithmeticError(f'naturality fails for d_{r} at {position}')
        return True
