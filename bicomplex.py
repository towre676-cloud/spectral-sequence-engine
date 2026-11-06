"""First-quadrant double complexes and their two exact filtration spectral sequences.

Horizontal d_h:(p,q)->(p-1,q), vertical d_v:(p,q)->(p,q-1).
We require d_h d_v + d_v d_h = 0. The total differential is d_h+d_v;
column filtration by p has E_1[p,q]=H_q(C[p,*],d_v).

A double complex is an algebraic input, not automatically a chain model
for a geometric fibration or the Serre spectral sequence of that fibration.
"""
from __future__ import annotations

from collections import defaultdict
from typing import Mapping

from spectral_sequence import FilteredChainComplex, _prime, _rank, _Span, _kernel


def tensor_double_complex(left: FilteredChainComplex, right: FilteredChainComplex):
    """Tensor-product bicomplex of two finite prime-field chain complexes.

    Horizontal d_h=d_left ⊗ 1; vertical d_v=(-1)^deg(left) 1 ⊗ d_right.
    The total differential is the standard Koszul-signed tensor differential.
    Existing filtration births on the inputs do not enter this construction:
    the new column filtration is their *homological degree* on the left.
    """
    if not isinstance(left,FilteredChainComplex) or not isinstance(right,FilteredChainComplex):
        raise ValueError('both tensor factors must be FilteredChainComplex instances')
    if left.prime!=right.prime:
        raise ValueError('tensor factors require the same prime field')
    names={(a,b):repr((a,b)) for a in left.names for b in right.names}
    cells=[]
    for a in left.names:
        p,_,da=left.data[a]
        for b in right.names:
            q,_,db=right.data[b]
            horizontal={names[(target,b)]:v for target,v in da.items()}
            vertical={names[(a,target)]:(-1)**p*v for target,v in db.items()}
            cells.append((names[(a,b)],p,q,vertical,horizontal))
    return DoubleComplex(cells,prime=left.prime)


class DoubleComplex:
    """Generators (name, p, q, vertical_boundary, horizontal_boundary)."""
    def __init__(self, generators, prime=2):
        if not _prime(prime):
            raise ValueError('coefficient characteristic must be prime')
        self.prime = prime
        self.data = {}
        for entry in generators:
            if len(entry) != 5:
                raise ValueError('each generator must have five fields')
            name, p, q, dv, dh = entry
            if not isinstance(name, str) or not name or name in self.data:
                raise ValueError('generator names must be unique nonempty strings')
            if any(type(i) is not int or i < 0 for i in (p, q)):
                raise ValueError('bidegrees must be nonnegative integers')
            if any(not isinstance(m, Mapping) or any(not isinstance(k,str) or type(v) is not int
                                                     for k,v in m.items()) for m in (dv,dh)):
                raise ValueError('boundary must be a string-to-integer mapping')
            self.data[name] = (p, q,
                               {k: v % prime for k,v in dv.items() if v % prime},
                               {k: v % prime for k,v in dh.items() if v % prime})
        self.names = tuple(sorted(self.data, key=lambda name: (self.data[name][0],self.data[name][1],name)))
        for name, (p,q,dv,dh) in self.data.items():
            for target in dv:
                if target not in self.data or self.data[target][:2] != (p,q-1):
                    raise ValueError(f'{name}: vertical boundary has wrong bidegree or unknown generator')
            for target in dh:
                if target not in self.data or self.data[target][:2] != (p-1,q):
                    raise ValueError(f'{name}: horizontal boundary has wrong bidegree or unknown generator')
        def compose(first, second):
            out = defaultdict(int)
            for intermediate, x in second.items():
                for target,y in self.data[intermediate][first].items():
                    out[target] = (out[target] + x*y) % prime
            return {k:v for k,v in out.items() if v}
        for name, (_,_,dv,dh) in self.data.items():
            vv = compose(2,dv)
            hh = compose(3,dh)
            hv = compose(3,dv)
            vh = compose(2,dh)
            if vv or hh or any((hv.get(k,0)+vh.get(k,0))%prime
                                for k in hv.keys()|vh.keys()):
                raise ValueError(f'{name}: double-complex square/anticommutation identities fail')
        self._column = None
        self._row = None

    def _total(self, horizontal=True):
        chain=[]
        for name in self.names:
            p,q,dv,dh = self.data[name]
            boundary=defaultdict(int)
            for mapping in (dv,dh):
                for target,x in mapping.items():
                    boundary[target] = (boundary[target]+x)%self.prime
            chain.append((name, p+q, p if horizontal else q,
                          {k:v for k,v in boundary.items() if v}))
        return FilteredChainComplex(chain,prime=self.prime)

    def column_spectral_sequence(self):
        """Filter by horizontal degree p: E_1 is vertical homology."""
        if self._column is None:
            self._column=self._total(horizontal=True)
        return self._column

    def row_spectral_sequence(self):
        """Filter by vertical degree q: E_1 is horizontal homology."""
        if self._row is None:
            self._row=self._total(horizontal=False)
        return self._row

    def _directional_e1(self, direction):
        """Independent associated-graded rank computation, without page engine."""
        groups=defaultdict(list)
        for name, (p,q,*_) in self.data.items():
            if direction=='vertical':
                groups[(p,q)].append(name)
            else:
                groups[(q,p)].append(name)
        result={}
        for grade,degree in groups:
            names=groups[(grade,degree)]
            lower=groups.get((grade,degree-1),[])
            upper=groups.get((grade,degree+1),[])
            component=2 if direction=='vertical' else 3
            outgoing=[tuple(self.data[name][component].get(target,0) for target in lower)
                      for name in names]
            incoming=[tuple(self.data[name][component].get(target,0) for target in names)
                      for name in upper]
            dim=len(names)-_rank(outgoing,self.prime,len(lower))-_rank(incoming,self.prime,len(names))
            if dim < 0:
                raise ArithmeticError('negative directional homology dimension')
            if dim:
                result[(grade,degree)] = dim
        return result

    def _directional_e2(self, direction):
        """Independent E2: homology of the induced cross-axis d1 on E1.

        This computes vertical homology, induced horizontal maps, then their
        homology, without consulting FilteredChainComplex.page() or d_r().
        """
        groups=defaultdict(list)
        for name,(p,q,*_) in self.data.items():
            grade,degree = (p,q) if direction=='vertical' else (q,p)
            groups[(grade,degree)].append(name)
        for names in groups.values():
            names.sort()
        along=2 if direction=='vertical' else 3
        across=3 if direction=='vertical' else 2
        reps={}
        bounds={}
        for (grade,degree),names in groups.items():
            lower=groups.get((grade,degree-1),())
            upper=groups.get((grade,degree+1),())
            incoming=_Span(len(names),self.prime)
            for name in upper:
                incoming.add(tuple(self.data[name][along].get(target,0) for target in names))
            boundary_basis=tuple(incoming.vectors)
            out_cols=[tuple(self.data[name][along].get(target,0) for target in lower)
                      for name in names]
            cycles=_kernel(out_cols,self.prime,len(lower))
            chosen=[]
            for cycle in cycles:
                if incoming.add(cycle):
                    chosen.append(cycle)
            bounds[(grade,degree)]=boundary_basis
            reps[(grade,degree)]=tuple(chosen)
        d1={}
        for (grade,degree),sources in reps.items():
            source_names=groups[(grade,degree)]
            target=(grade-1,degree)
            target_names=groups.get(target,())
            target_index={name:i for i,name in enumerate(target_names)}
            ambient=_Span(len(target_names),self.prime)
            for vec in (*bounds.get(target,()),*reps.get(target,())):
                ambient.add(vec)
            offset=len(bounds.get(target,()))
            columns=[]
            for source in sources:
                image=[0]*len(target_names)
                for coefficient,name in zip(source,source_names):
                    for cell,v in self.data[name][across].items():
                        image[target_index[cell]]=(image[target_index[cell]]+coefficient*v)%self.prime
                columns.append(ambient.coordinates(tuple(image))[offset:])
            d1[(grade,degree)]=tuple(columns)
        result={}
        for (grade,degree),basis in reps.items():
            dim=len(basis)
            outgoing=_rank(d1[(grade,degree)],self.prime,len(reps.get((grade-1,degree),())))
            incoming=_rank(d1.get((grade+1,degree),()),self.prime,dim)
            remaining=dim-outgoing-incoming
            if remaining < 0:
                raise ArithmeticError('cross-axis differential has impossible rank')
            if remaining:
                result[(grade,degree)]=remaining
        return result

    def verify(self):
        """Check two independently assembled E_1 pages and common abutment."""
        col=self.column_spectral_sequence()
        row=self.row_spectral_sequence()
        vertical=self._directional_e1('vertical')
        horizontal=self._directional_e1('horizontal')
        if col.page(1) != vertical:
            raise ArithmeticError('vertical homology does not match E_1')
        if row.page(1) != horizontal:
            raise ArithmeticError('horizontal homology does not match E_1')
        if col.page(2) != self._directional_e2('vertical'):
            raise ArithmeticError('column E_2 disagrees with independent double-complex E_2')
        if row.page(2) != self._directional_e2('horizontal'):
            raise ArithmeticError('row E_2 disagrees with independent double-complex E_2')
        if not col.verify() or not row.verify():
            raise ArithmeticError('spectral-sequence verification failed')
        if col.betti() != row.betti():
            raise ArithmeticError('two filtrations must give the same total homology')
        return True
