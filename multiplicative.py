"""Verified filtered differential graded associative algebras over prime fields.

Products on E_r are computed from chain representatives, not stipulated page maps.
This module uses homological grading, so d(xy)=d(x)y+(-1)^|x| x d(y).
"""
from itertools import product
from spectral_sequence import _Span, _linear_combination, _rank


class FilteredDGA:
    """Finite non-unital DGA, optionally requiring graded commutativity.

    products: dictionary mapping (source_generator, source_generator) to a
    dictionary of output generators and integer coefficients; omitted pairs
    have zero product. The finite-field underlying complex is validated first.
    """
    def __init__(self, chain, products, graded_commutative=False):
        if not isinstance(products, dict):
            raise ValueError('products must be a dictionary')
        self.chain, self.prime = chain, chain.prime
        self.products = {}
        self.graded_commutative = graded_commutative
        for pair, outputs in products.items():
            if not isinstance(pair, tuple) or len(pair) != 2 or any(x not in chain.data for x in pair):
                raise ValueError('product keys must be pairs of generators')
            if not isinstance(outputs, dict):
                raise ValueError('product values must be dictionaries')
            n = sum(chain.data[x][0] for x in pair)
            p = sum(chain.data[x][1] for x in pair)
            v = [0] * chain._width(n)
            for name, coefficient in outputs.items():
                if name not in chain.data or isinstance(coefficient, bool) or not isinstance(coefficient, int):
                    raise ValueError('invalid product coefficient or generator')
                deg, birth, _ = chain.data[name]
                if deg != n or (coefficient % self.prime and birth > p):
                    raise ValueError('product fails degree or filtration bound')
                v[chain.loc[n][name]] = (v[chain.loc[n][name]] + coefficient) % self.prime
            self.products[pair] = tuple(v)
        self.verify_chain_algebra()

    def multiply(self, left_degree, left, right_degree, right):
        c = self.chain
        if len(left) != c._width(left_degree) or len(right) != c._width(right_degree):
            raise ValueError('incorrect chain input dimensions')
        n = left_degree + right_degree
        out = [0] * c._width(n)
        for i, u in enumerate(left):
            for j, v in enumerate(right):
                if u % self.prime and v % self.prime:
                    pair = (c.by_degree[left_degree][i], c.by_degree[right_degree][j])
                    column = self.products.get(pair, (0,)*len(out))
                    for k, z in enumerate(column):
                        out[k] = (out[k] + u*v*z) % self.prime
        return tuple(out)

    def verify_chain_algebra(self):
        c = self.chain
        for a,b in product(c.names, repeat=2):
            da, db = c.data[a][0], c.data[b][0]
            va, vb = c.unit(da,a), c.unit(db,b)
            ab = self.multiply(da,va,db,vb)
            lhs = c.boundary(da+db,ab)
            first = self.multiply(da-1,c.boundary(da,va),db,vb)
            second = self.multiply(da,va,db-1,c.boundary(db,vb))
            rhs = tuple((x + (-1)**da*y)%self.prime for x,y in zip(first,second))
            if lhs != rhs:
                raise ValueError(f'Leibniz identity fails at {a}, {b}')
            if self.graded_commutative:
                ba = self.multiply(db,vb,da,va)
                if ab != tuple(((-1)**(da*db)*z)%self.prime for z in ba):
                    raise ValueError(f'graded commutativity fails at {a}, {b}')
        for a,b,d in product(c.names,repeat=3):
            da,db,dd = (c.data[x][0] for x in (a,b,d))
            va,vb,vd = (c.unit(n,x) for n,x in ((da,a),(db,b),(dd,d)))
            left = self.multiply(da+db,self.multiply(da,va,db,vb),dd,vd)
            right = self.multiply(da,va,db+dd,self.multiply(db,vb,dd,vd))
            if left != right:
                raise ValueError(f'associativity fails at {a}, {b}, {d}')
        return True

    def _reps(self,r,p,n):
        c = self.chain
        if r:
            return c._space(r,p,n).representatives
        return tuple(c.unit(n,name) for name in c.by_degree.get(n,()) if c.data[name][1]==p)

    def _coordinates(self,r,p,n,v):
        c = self.chain
        if r:
            target = c._space(r,p,n)
            s = _Span(c._width(n),self.prime)
            for z in (*target.boundaries,*target.representatives):
                s.add(z)
            return s.coordinates(v)[len(target.boundaries):]
        return tuple(v[c.loc[n][name]] for name in c.by_degree.get(n,()) if c.data[name][1]==p)

    def page_product(self,r,p,n,s,m):
        """Columns for E_r(p,n) tensor E_r(s,m) -> E_r(p+s,n+m).

        The pairs are ordered lexicographically by the first and second basis
        indices. Returns a tuple of columns, including zero columns.
        """
        if not isinstance(r,int) or r<0:
            raise ValueError('page must be a nonnegative integer')
        return tuple(self._coordinates(r,p+s,n+m,self.multiply(n,a,m,b))
                     for a in self._reps(r,p,n) for b in self._reps(r,s,m))

    def verify_page_well_defined(self, through=None):
        """Check that changing a class by a page denominator changes no product.

        This directly checks the subtle multiplicative descent condition on
        representatives, rather than relying solely on chain-level Leibniz.
        """
        c=self.chain
        if through is None:
            through=c.spread+1
        for r in range(1,through+1):
            positions=[(p,p+q) for p,q in c.page(r)]
            for p,n in positions:
                a_space=c._space(r,p,n)
                for s,m in positions:
                    b_space=c._space(r,s,m)
                    target=c._space(r,p+s,n+m)
                    denom=_Span(c._width(n+m),self.prime)
                    for v in target.boundaries:
                        denom.add(v)
                    for a in a_space.boundaries:
                        for b in b_space.representatives:
                            denom.coordinates(self.multiply(n,a,m,b))
                    for a in a_space.representatives:
                        for b in b_space.boundaries:
                            denom.coordinates(self.multiply(n,a,m,b))
        return True

    def verify_page_leibniz(self,through=None):
        """Check product Leibniz on every chosen page, including d_0."""
        c=self.chain
        if through is None:
            through=c.spread+1
        for r in range(through+1):
            positions=[(p,p+q) for p,q in c.page(r)]
            for p,n in positions:
                for s,m in positions:
                    for a in self._reps(r,p,n):
                        for b in self._reps(r,s,m):
                            ab=self.multiply(n,a,m,b)
                            left=self._coordinates(r,p+s-r,n+m-1,c.boundary(n+m,ab))
                            da=c.boundary(n,a); db=c.boundary(m,b)
                            rhs1=self.multiply(n-1,da,m,b)
                            rhs2=self.multiply(n,a,m-1,db)
                            rhs=tuple((x+(-1)**n*y)%self.prime for x,y in zip(rhs1,rhs2))
                            right=self._coordinates(r,p+s-r,n+m-1,rhs)
                            if left != right:
                                raise ArithmeticError(f'page Leibniz fails on E_{r}')
        return True
