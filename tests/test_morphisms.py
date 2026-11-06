import unittest
from morphisms import FilteredChainMap
from spectral_sequence import FilteredChainComplex


class FilteredMapTests(unittest.TestCase):
    def fixture(self, prime=3):
        return FilteredChainComplex([
            ('a',0,0,{}), ('b',0,1,{}),
            ('u',1,1,{'a':1,'b':1}), ('v',1,2,{'b':1})],prime)

    def test_identity_zero_and_scaling(self):
        for prime in (2,3,5):
            c = self.fixture(prime)
            for scalar in (0,1,2):
                f = FilteredChainMap(c,c,{
                    name: {name: scalar} for name in c.names})
                self.assertTrue(f.verify_naturality())
                for r in range(4):
                    for pos,cols in f.on_page(r).items():
                        if scalar == 0:
                            self.assertTrue(all(not any(x) for x in cols))

    def test_inclusion_and_filtration_constraint(self):
        c = self.fixture()
        sub = FilteredChainComplex([('a',0,0,{})],3)
        f = FilteredChainMap(sub,c,{'a':{'a':1}})
        self.assertTrue(f.verify_naturality())
        with self.assertRaises(ValueError):
            FilteredChainMap(sub,c,{'a':{'b':1}})

    def test_bad_chain_map_and_bad_degree(self):
        c = self.fixture()
        with self.assertRaises(ValueError):
            FilteredChainMap(c,c,{'v':{'v':1}})  # d(v) maps to zero
        with self.assertRaises(ValueError):
            FilteredChainMap(c,c,{'a':{'u':1}})
        with self.assertRaises(ValueError):
            FilteredChainMap(c,c,{'nonexistent':{}})
        other = self.fixture(2)
        with self.assertRaises(ValueError):
            FilteredChainMap(c,other,{})

    def test_map_between_differential_lift_models(self):
        c = self.fixture(5)
        # A coordinate change of u and v is not a chain map unless it
        # respects their boundaries; map a -> a, b -> b, u -> u, v -> v.
        mapped = FilteredChainMap(c,c,{name:{name:1} for name in c.names})
        self.assertTrue(mapped.verify_naturality())
        self.assertEqual(mapped.on_page(2)[(2,-1)], ((1,),))


if __name__ == '__main__':
    unittest.main()
