"""Singly graded exact couples and their derived graded differentials.

D_n --i--> D_n --j--> E_n --k--> D_(n-1), with exactness
im(i)=ker(j), im(j)=ker(k), im(k)=ker(i) at each relevant degree.
The E differential d_r always lowers homological degree by one.

This is a *singly graded algebraic* exact couple. Bigraded shifts,
filtration indices, and connecting homomorphisms from a specific fibration
must be supplied by a further geometric construction.
"""
from collections import defaultdict

from exact_couple import ExactCouple
from spectral_sequence import _Span, _kernel, _rank


def _homogeneous_degree(vector, labels):
    seen={labels[i] for i,x in enumerate(vector) if x}
    if len(seen)!=1:
        raise ArithmeticError('an exact-couple derived basis vector is not homogeneous')
    return seen.pop()


class GradedExactCouple:
    """Finite-field exact couple with ordered homogeneous D,E basis labels.

    Matrices i:D->D, j:D->E, k:E->D use column-major tuples.
    i,j preserve homological degree; k lowers it by one. All three
    exactness identities are checked by the underlying exact-couple engine.
    """
    def __init__(self, prime, d_degrees, e_degrees, i, j, k):
        self.d_degrees=tuple(d_degrees)
        self.e_degrees=tuple(e_degrees)
        if any(type(n) is not int for n in self.d_degrees+self.e_degrees):
            raise ValueError('homological degrees must be integers')
        self.couple=ExactCouple(prime,len(self.d_degrees),len(self.e_degrees),i,j,k)
        self.prime=prime
        self.verify()

    def verify(self):
        c=self.couple
        for columns,source,target,shift,name in (
            (c.i,self.d_degrees,self.d_degrees,0,'i'),
            (c.j,self.d_degrees,self.e_degrees,0,'j'),
            (c.k,self.e_degrees,self.d_degrees,-1,'k'),
        ):
            for j,col in enumerate(columns):
                for i,x in enumerate(col):
                    if x and target[i] != source[j] + shift:
                        raise ValueError(f'{name} has wrong degree at matrix entry ({i},{j})')
        return c.verify()

    def differential(self):
        """Homogeneous degree -1 matrix j*k on E, with global basis order."""
        return self.couple.differential()

    def differential_by_degree(self):
        d=self.differential()
        out={}
        for n in sorted(set(self.e_degrees)):
            sources=[i for i,x in enumerate(self.e_degrees) if x==n]
            targets=[i for i,x in enumerate(self.e_degrees) if x==n-1]
            columns=tuple(tuple(d[j][i] for i in targets) for j in sources)
            out[n]={'target_degree':n-1,'columns':columns,
                    'rank':_rank(columns,self.prime,len(targets))}
        return out

    def page_dimensions(self):
        out=defaultdict(int)
        for n in self.e_degrees:
            out[n]+=1
        return dict(sorted(out.items()))

    def derive(self):
        """Construct graded derived couple via images and cycle/quotient reps."""
        c=self.couple
        derived=c.derive()
        ds=_Span(c.dim_d,self.prime)
        d_degrees=[]
        for col in c.i:
            if ds.add(col):
                d_degrees.append(_homogeneous_degree(col,self.d_degrees))
        d=c.differential()
        cycles=_kernel(d,self.prime,c.dim_e)
        es=_Span(c.dim_e,self.prime)
        for col in d:
            es.add(col)
        e_degrees=[]
        for vec in cycles:
            if es.add(vec):
                e_degrees.append(_homogeneous_degree(vec,self.e_degrees))
        if len(d_degrees)!=derived.dim_d or len(e_degrees)!=derived.dim_e:
            raise ArithmeticError('grading reconstruction mismatches derived quotient ranks')
        return GradedExactCouple(self.prime,d_degrees,e_degrees,
                                 derived.i,derived.j,derived.k)

    def derived_pages(self,max_steps=None):
        if max_steps is None:
            max_steps=len(self.d_degrees)+2
        if type(max_steps) is not int or max_steps<0:
            raise ValueError('max_steps must be a nonnegative integer')
        state=self
        history=[]
        for r in range(1,max_steps+2):
            data=state.differential_by_degree()
            history.append({'r':r, 'D_dimensions':dict(sorted((n,state.d_degrees.count(n))
                                           for n in set(state.d_degrees))),
                            'E_dimensions':state.page_dimensions(),
                            'd_ranks':{n:v['rank'] for n,v in data.items() if v['rank']}})
            if state.couple.dim_e == 0:
                break
            if r <= max_steps:
                following=state.derive()
                old=state.page_dimensions()
                expected={n:dim - data.get(n,{}).get('rank',0) - data.get(n+1,{}).get('rank',0)
                          for n,dim in old.items()}
                expected={n:dim for n,dim in expected.items() if dim}
                if following.page_dimensions()!=expected:
                    raise ArithmeticError('graded E_(r+1) does not match homology of E_r')
                state=following
        return history
