"""Reproducible non-interactive geometric spectral-sequence comparison checks.

Run: python -m research.geometric_benchmarks
All assertions are exact; report is not a formal topology theorem.
"""
from bockstein import LiftedBockstein
from simplicial_maps import _example


def main():
    hopf=_example('hopf')
    product=_example('product-sphere')
    assert hopf.domain.f_vector()==product.domain.f_vector()=={0:12,1:48,2:72,3:36}
    for prime in (2,3,5,7):
        h,p=hopf.chain_complex(prime), product.chain_complex(prime)
        assert h.page(2)==p.page(2)=={(0,0):1,(0,1):1,(2,0):1,(2,1):1}
        assert h.differential(2)[(2,0)]['rank']==1
        assert p.differential(2)[(2,0)]['rank']==0
        assert h.betti()=={0:1,3:1}
        assert p.betti()=={0:1,1:1,2:1,3:1}
        assert h.verify() and p.verify()
    print('Hopf vs trivial product: 168 simplices each, identical E2, d2 ranks 1 vs 0; PASS')

    torus=_example('torus')
    klein=_example('klein')
    assert torus.domain.integral_chain_complex().homology()[1]=={'free_rank':2,'torsion':()}
    assert klein.domain.integral_chain_complex().homology()[1]=={'free_rank':1,'torsion':(2,)}
    for prime in (2,3,5,7):
        if prime==2:
            assert klein.chain_complex(prime).betti()=={0:1,1:2,2:1}
        else:
            assert klein.chain_complex(prime).betti()=={0:1,1:1}
        fiber=klein.fiber_over_vertex(0)
        # The independent local-system calculation is in unit tests using
        # the original pre-gluing fiber automorphism.
        assert fiber.chain_complex(prime)[0].betti()=={0:1,1:1}
        assert klein.verify(prime)
    b=LiftedBockstein(klein.domain.integral_chain_complex(),2)
    assert b.verify() and any(any(x for x in col) for col in b.differential(1)[2])
    print('Torus vs Klein: integral Z/2 torsion and nonzero mod-2 Bockstein; PASS')

    twisted=_example('twisted-sphere-torus')
    assert twisted.domain.integral_chain_complex().homology()[2]=={'free_rank':0,'torsion':(2,)}
    assert twisted.chain_complex(3).betti()=={0:1,1:1}
    assert twisted.chain_complex(2).betti()=={0:1,1:1,2:1,3:1}
    print('Twisted S2-bundle over S1: H2(Z)=Z/2 and odd/even field distinction; PASS')
    print('Geometric benchmark suite: PASS')


if __name__=='__main__':
    main()
