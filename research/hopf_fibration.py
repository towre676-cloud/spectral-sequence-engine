"""Published simplicial Hopf map, filtered by the base 2-skeleton.

K. V. Madahar & K. S. Sarkaria, Geometriae Dedicata 2000,
via the Macaulay2 SpectralSequences documentation:
https://www.macaulay2.com/doc/Macaulay2-1.25.11/share/doc/Macaulay2/SpectralSequences/html/___Computing_spthe_sp__Serre_sp__Spectral_sp__Sequence_spassociated_spto_spa_sp__Hopf_sp__Fibration.html

S^1 -> S^3 -> S^2. Vertices a_i,b_i,c_i,d_i project to
vertices a,b,c,d of the tetrahedral boundary. We enumerate all faces
of the 36 listed tetrahedra; filtration p of each simplex is the
dimension of its projected base simplex. This is the filtration from
preimages of the base skeleta, not a hand-engineered page table.
"""
from itertools import combinations
from collections import Counter

from spectral_sequence import FilteredChainComplex
from simplicial_maps import SimplicialComplex, SimplicialMap

TETRAHEDRA = """
a_0 b_0 b_1 c_1
a_0 b_0 c_0 c_1
a_0 a_1 b_1 c_1
b_0 b_1 c_1 d_1
b_0 c_0 c_1 d_2
a_0 a_1 c_1 d_2
a_0 c_0 c_1 d_2
b_0 c_1 d_1 d_2
b_1 c_1 c_2 a_2
b_1 c_1 a_1 a_2
b_1 b_2 c_2 a_2
c_1 c_2 a_2 d_1
c_1 a_1 a_2 d_2
b_1 b_2 a_2 d_2
b_1 a_1 a_2 d_2
c_1 a_2 d_1 d_2
c_2 a_2 a_0 b_0
c_2 a_2 b_2 b_0
c_2 c_0 a_0 b_0
a_2 a_0 b_0 d_1
a_2 b_2 b_0 d_2
c_2 c_0 b_0 d_2
c_2 b_2 b_0 d_2
a_2 b_0 d_1 d_2
a_0 b_0 b_1 d_1
a_0 b_1 d_0 d_1
b_1 c_1 c_2 d_1
b_1 c_2 d_0 d_1
a_0 a_2 c_2 d_1
a_0 c_2 d_0 d_1
a_0 b_1 d_0 d_2
a_0 a_1 b_1 d_2
b_1 c_2 d_0 d_2
b_1 b_2 c_2 d_2
a_0 c_2 d_0 d_2
a_0 c_0 c_2 d_2
"""


def simplicial_hopf_filtration():
    """Return all faces with their induced base-skeletal filtration births."""
    tetrahedra=[tuple(row.split()) for row in TETRAHEDRA.splitlines() if row.strip()]
    if len(tetrahedra)!=36 or len(set(tuple(sorted(t)) for t in tetrahedra))!=36:
        raise ArithmeticError('expected 36 distinct published tetrahedra')
    if any(len(t)!=4 or len(set(t))!=4 for t in tetrahedra):
        raise ArithmeticError('incorrect tetrahedral facet')
    simplices=set()
    for tet in tetrahedra:
        for k in range(1,5):
            simplices.update(combinations(sorted(tet),k))
    return [(s,len({v[0] for v in s})-1) for s in sorted(simplices)]


def hopf_simplicial_map():
    """Feed the published 36 facets and vertex projection to the generic engine.

    The separate simplicial_hopf_filtration() above remains an independent
    combinatorial reference for the expected birth labels of all 168 faces.
    """
    tetrahedra=[tuple(row.split()) for row in TETRAHEDRA.splitlines() if row.strip()]
    if len(tetrahedra)!=36 or len(set(frozenset(t) for t in tetrahedra))!=36:
        raise ArithmeticError('published triangulation must have 36 distinct facets')
    total=SimplicialComplex(tetrahedra)
    base=SimplicialComplex([tuple(v) for v in ('abc','abd','acd','bcd')])
    projection=SimplicialMap(total,base,{v:v[0] for v in total.vertices})
    if set(projection.induced_filtration())!=set(simplicial_hopf_filtration()):
        raise ArithmeticError('automatically induced and published Hopf filtrations disagree')
    return projection


def hopf_complex(prime=2):
    """Hopf chain complex constructed by the generic simplicial-map adapter."""
    return hopf_simplicial_map().chain_complex(prime)


def report(prime=2):
    c=hopf_complex(prime)
    e0=c.page(0)
    e1=c.page(1)
    e2=c.page(2)
    ranks={r:sum(v['rank'] for v in c.differential(r).values()) for r in (0,1,2)}
    return {'prime':prime,'n_simplices':len(c.names),'e0':e0,'e1':e1,'e2':e2,
            'd_ranks':ranks,'e3':c.page(3),'betti':c.betti(),'verified':c.verify()}


if __name__=='__main__':
    for p in (2,3,5):
        print(report(p))
