import unittest
from spectral_sequence import FilteredChainComplex
from multiplicative import FilteredDGA
from exact_couple import ExactCouple


class MultiplicativeTests(unittest.TestCase):
    def exterior(self,prime):
        # Exterior algebra on two closed degree-one classes and their degree-two product.
        c=FilteredChainComplex([('one',0,0,{}),('x',1,1,{}),('y',1,1,{}),('xy',2,2,{})],prime)
        prod={('one',v):{v:1} for v in c.names}
        prod.update({(v,'one'):{v:1} for v in c.names})
        prod[('x','y')]={'xy':1}; prod[('y','x')]={'xy':-1}
        return FilteredDGA(c,prod,graded_commutative=True)

    def test_exterior_algebras_odd_characteristic(self):
        for prime in (2,3,5,7):
            dga=self.exterior(prime)
            self.assertEqual(dga.page_product(1,1,1,1,1),((0,), (1,),((-1)%prime,), (0,)))
            self.assertTrue(dga.verify_page_leibniz())
            self.assertTrue(dga.verify_page_well_defined())
            self.assertTrue(dga.chain.verify())

    def test_nontrivial_differential_and_product(self):
        # Unital algebra F_p[u]/u² with |u|=0, plus contractible ideal a<-b.
        # one*a=a, one*b=b; other pair-products vanish.
        for prime in (2,3,5):
            c=FilteredChainComplex([('one',0,0,{}),('a',0,0,{}),('b',1,2,{'a':1})],prime)
            prod={('one',v):{v:1} for v in c.names}
            prod.update({(v,'one'):{v:1} for v in c.names})
            dga=FilteredDGA(c,prod,graded_commutative=True)
            self.assertTrue(dga.verify_page_leibniz())
            self.assertTrue(dga.verify_page_well_defined())
            self.assertEqual(c.differential(2)[(2,-1)]['rank'],1)

    def test_dga_rejects_bad_products(self):
        c=FilteredChainComplex([('a',0,0,{}),('b',1,2,{'a':1})],3)
        with self.assertRaises(ValueError):
            FilteredDGA(c,{('a','b'):{'a':1}})
        with self.assertRaises(ValueError):
            FilteredDGA(c,{('b','a'):{'b':1}})  # fails Leibniz
        with self.assertRaises(ValueError):
            FilteredDGA(c,{('a','a'):{'b':1}})


class ExactCoupleTests(unittest.TestCase):
    def test_exact_couple_with_nonzero_d1(self):
        for prime in (2,3,5):
            c=ExactCouple(prime,1,2,[(0,)],[(1,0)],[(0,),(1,)])
            self.assertEqual(c.differential(),((0,0),(1,0)))
            self.assertEqual(c.e2_dimension(),0)

    def test_hidden_second_differential_from_derived_couple(self):
        # i(b)=a, j(b)=x, k(y)=a. d1=jk=0 but d2 sends [y] -> [x].
        for prime in (2,3,5):
            first=ExactCouple(prime,2,2,[(0,0),(1,0)],
                              [(0,0),(1,0)],[(0,0),(1,0)])
            self.assertEqual(first.differential(),((0,0),(0,0)))
            second=first.derive()
            self.assertEqual((second.dim_d,second.dim_e),(1,2))
            self.assertEqual(second.differential(),((0,0),(1,0)))
            self.assertEqual(second.e2_dimension(),0)
            third=second.derive()
            self.assertEqual(third.dim_e,0)
            self.assertEqual([x['d_rank'] for x in first.derived_pages()], [0,1,0])

    def test_zero_couple(self):
        c=ExactCouple(3,0,0,[],[],[])
        self.assertEqual(c.e2_dimension(),0)

    def test_nonexact_rejected(self):
        with self.assertRaises(ValueError):
            ExactCouple(2,1,1,[(0,)],[(0,)],[(0,)])
        with self.assertRaises(ValueError):
            ExactCouple(2,1,2,[(0,)],[(1,0)],[(1,),(0,)])
        with self.assertRaises(ValueError):
            ExactCouple(4,0,0,[],[],[])


if __name__ == '__main__':
    unittest.main()
