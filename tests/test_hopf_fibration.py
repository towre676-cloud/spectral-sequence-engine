"""Published Hopf-map triangulation, independent Macaulay2 page checks."""
import unittest

from research.hopf_fibration import simplicial_hopf_filtration, hopf_complex


class HopfFibrationTests(unittest.TestCase):
    def test_published_facets_and_filtration_counts(self):
        cells=simplicial_hopf_filtration()
        self.assertEqual(len(cells),168)
        self.assertEqual(len({v for simplex,_ in cells for v in simplex}),12)
        self.assertEqual(max(len(simplex) for simplex,_ in cells),4)
        for simplex,p in cells:
            self.assertEqual(p,len({v[0] for v in simplex})-1)
            self.assertLessEqual(p,2)
        for prime in (2,3,5):
            chain=hopf_complex(prime)
            self.assertEqual(chain.page(0),{
                (0,0):12,(0,1):12,(1,0):36,(1,1):36,(2,0):36,(2,1):36})

    def test_exact_published_E1_E2_E3_and_hopf_transgression(self):
        # Independent published reference: Macaulay2 SpectralSequences,
        # Hopf fibration worked example, over QQ (rank comparison).
        for prime in (2,3,5,7):
            chain=hopf_complex(prime)
            self.assertEqual(chain.page(1),{
                (0,0):4,(0,1):4,(1,0):6,(1,1):6,(2,0):4,(2,1):4})
            self.assertEqual(chain.page(2),{
                (0,0):1,(0,1):1,(2,0):1,(2,1):1})
            d2=chain.differential(2)
            self.assertEqual(d2[(2,0)]['target'],(0,1))
            self.assertEqual(d2[(2,0)]['rank'],1)
            self.assertEqual(chain.page(3),{(0,0):1,(2,1):1})
            self.assertEqual(chain.betti(),{0:1,3:1})
            self.assertEqual({r:sum(v['rank'] for v in chain.differential(r).values())
                             for r in range(3)}, {0:70,1:12,2:1})
            self.assertTrue(chain.verify())


if __name__=='__main__':
    unittest.main()
