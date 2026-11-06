import itertools
import random
import unittest
from integral import IntegralChainComplex, _diagonalize, _integer_inverse

class IntegralHomologyTests(unittest.TestCase):
    def test_moore_spaces_and_free_groups(self):
        for n in (2,3,4,6,12):
            c=IntegralChainComplex([1,1,1],{2:[[n]]})
            self.assertEqual(c.homology(),{0:{'free_rank':1,'torsion':()},
                              1:{'free_rank':0,'torsion':(n,)},2:{'free_rank':0,'torsion':()}})
        c=IntegralChainComplex([1,2,0],{})
        self.assertEqual(c.homology()[1],{'free_rank':2,'torsion':()})

    def test_normalize_coprime_and_nondividing_torsion(self):
        c=IntegralChainComplex([0,2,2],{2:[[2,0],[0,3]]})
        self.assertEqual(c.homology()[1]['torsion'],(6,))
        c=IntegralChainComplex([0,2,2],{2:[[4,0],[0,6]]})
        self.assertEqual(c.homology()[1]['torsion'],(2,12))

    def test_nontrivial_kernel_basis_and_boundary(self):
        # d1=[2,4] and d2=[-2,1]^T => H1=0 and H0=Z/2.
        c=IntegralChainComplex([1,2,1],{1:[[2,4]],2:[[-2],[1]]})
        self.assertEqual(c.homology()[0],{'free_rank':0,'torsion':(2,)})
        self.assertEqual(c.homology()[1],{'free_rank':0,'torsion':()})

    def test_invalid_complex(self):
        with self.assertRaises(ValueError):
            IntegralChainComplex([1,1,1],{1:[[1]],2:[[1]]})
        with self.assertRaises(ValueError):
            IntegralChainComplex([1,1],{1:[[1,2]]})

    def test_diagonalization_unimodular(self):
        rng=random.Random(812)
        for rows,cols in itertools.product(range(5),repeat=2):
            for trial in range(16):
                a=[[rng.randrange(-5,6) for _ in range(cols)] for _ in range(rows)]
                diag,v,rank=_diagonalize(a,cols)
                inv=_integer_inverse(v)
                self.assertEqual([[sum(v[i][k]*inv[k][j] for k in range(cols)) for j in range(cols)] for i in range(cols)],
                                 [[int(i==j) for j in range(cols)] for i in range(cols)])
                self.assertEqual(sum(bool(diag[i][i]) for i in range(min(rows,cols))),rank)

if __name__=='__main__': unittest.main()
