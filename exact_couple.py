"""Un-graded finite-field exact-couple validator and first derived differential.

This is a finite-dimensional algebraic exact couple D -i-> D -j-> E -k-> D.
It checks exactness at all three vertices and computes d_1=j*k and E_2.
Graded shifts and derived-couple iteration are not implemented here.
"""
from spectral_sequence import _Span, _kernel, _rank, _linear_combination, _prime


class ExactCouple:
    def __init__(self, prime, dim_d, dim_e, i, j, k):
        if not _prime(prime):
            raise ValueError('prime field required')
        if not all(isinstance(v,int) and not isinstance(v,bool) and v>=0 for v in (dim_d,dim_e)):
            raise ValueError('dimensions must be nonnegative integers')
        self.prime=prime
        self.dim_d=dim_d
        self.dim_e=dim_e
        self.i=self._validate(i,dim_d,dim_d)
        self.j=self._validate(j,dim_d,dim_e)
        self.k=self._validate(k,dim_e,dim_d)
        self.verify()

    def _validate(self, columns, width, height):
        if len(columns)!=width or any(len(c)!=height for c in columns):
            raise ValueError('matrix shape does not match couple dimensions')
        if any(not isinstance(v,int) or isinstance(v,bool) for c in columns for v in c):
            raise ValueError('integer matrix entries required')
        return tuple(tuple(x%self.prime for x in c) for c in columns)

    def _compose(self, first, second, height):
        return tuple(_linear_combination(first,col,self.prime,height) for col in second)

    def verify(self):
        for a,b,height in ((self.j,self.i,self.dim_e),(self.k,self.j,self.dim_d),(self.i,self.k,self.dim_d)):
            if any(any(z) for z in self._compose(a,b,height)):
                raise ValueError('exact-couple adjacent composites must vanish')
        # exactness by dimension, given vanishing composite: rank(previous)+rank(next)=dim
        rk_i=_rank(self.i,self.prime,self.dim_d)
        rk_j=_rank(self.j,self.prime,self.dim_e)
        rk_k=_rank(self.k,self.prime,self.dim_d)
        if rk_i+rk_j!=self.dim_d or rk_j+rk_k!=self.dim_e or rk_k+rk_i!=self.dim_d:
            raise ValueError('im(previous) must equal ker(next) at each vertex')
        return True

    def differential(self):
        """d_1 : E -> E as column-major matrix."""
        return self._compose(self.j,self.k,self.dim_e)

    def e2_dimension(self):
        d=self.differential()
        if any(any(x) for x in self._compose(d,d,self.dim_e)):
            raise ArithmeticError('exact-couple derived differential is not square-zero')
        return self.dim_e-2*_rank(d,self.prime,self.dim_e)

    def derive(self):
        """Construct the derived exact couple on D'=im(i), E'=ker(d)/im(d).

        Basis choices are deterministic with the declared source-column order.
        A nonzero class is represented by a cycle in the original E.
        This construction is valid for finite ungraded exact couples.
        """
        p = self.prime
        dim_d, dim_e = self.dim_d, self.dim_e
        d = self.differential()

        # D' has a basis of actual images i(x), with chosen preimages x.
        ds = _Span(dim_d,p)
        preimages = []
        for index,col in enumerate(self.i):
            if ds.add(col):
                preimages.append(tuple(int(i==index) for i in range(dim_d)))

        # E' is cycles modulo boundaries; retain classes as actual E chains.
        cycles = _kernel(d,p,dim_e)
        es = _Span(dim_e,p)
        for v in d:
            es.add(v)
        boundary_count = len(es.vectors)
        reps = []
        for v in cycles:
            if es.add(v):
                reps.append(v)
        def eclass(v):
            return es.coordinates(v)[boundary_count:]

        # i' is the old i restricted to im i.
        i1 = [ds.coordinates(self._compose(self.i,(v,),dim_d)[0])
              for v in ds.vectors]
        # j'(i x) = [j(x)]; the ambiguity lies in im(jk)=im(d).
        j1 = [eclass(self._compose(self.j,(preimage,),dim_e)[0])
              for preimage in preimages]
        # k'([e])=k(e) in ker j = im i, because d(e)=0.
        k1 = [ds.coordinates(self._compose(self.k,(v,),dim_d)[0])
              for v in reps]
        return ExactCouple(p,len(ds.vectors),len(reps),i1,j1,k1)

    def derived_pages(self, max_steps=None):
        """Report the E dimensions and differential ranks until D vanishes/stabilizes.

        Finite exact couples with i invertible can have persistent D and zero E;
        we stop once the derived couple's data repeat, or use max_steps.
        """
        if max_steps is None:
            max_steps = self.dim_d + 2
        if not isinstance(max_steps,int) or max_steps<0:
            raise ValueError('max_steps must be nonnegative')
        pages = []
        current = self
        for r in range(1,max_steps+2):
            rank = _rank(current.differential(),current.prime,current.dim_e)
            pages.append({'r':r,'D_dimension':current.dim_d,
                          'E_dimension':current.dim_e,'d_rank':rank})
            if current.dim_e == 0 or (r>1 and pages[-1]['D_dimension']==pages[-2]['D_dimension'] and pages[-1]['E_dimension']==pages[-2]['E_dimension'] and rank==0):
                break
            if r<=max_steps:
                current=current.derive()
        return pages
