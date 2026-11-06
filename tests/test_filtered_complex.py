import unittest
from filtered_complex import FilteredComplex


def triangle():
    return FilteredComplex([
        ((0,), 0), ((1,), 0), ((2,), 0),
        ((0,1), 1), ((0,2), 1), ((1,2), 1),
        ((0,1,2), 2)
    ])


class FilteredComplexTests(unittest.TestCase):
    def test_boundaries_and_betti(self):
        c = triangle()
        self.assertEqual(c.betti(0), {0: 3})
        self.assertEqual(c.betti(1), {0: 1, 1: 1})
        self.assertEqual(c.betti(2), {0: 1})

    def test_intervals(self):
        c = triangle()
        self.assertEqual(c.intervals().count((0, 0, 1)), 2)
        self.assertIn((0, 0, None), c.intervals())
        self.assertIn((1, 1, 2), c.intervals())

    def test_e1_relative_groups(self):
        e1 = triangle().e1_ranks()
        self.assertEqual(e1, {(0, 0): 3, (1, 0): 3, (2, 0): 1})

    def test_invalid_filtration(self):
        with self.assertRaises(ValueError):
            FilteredComplex([((0,1), 0), ((0,), 0)])
        with self.assertRaises(ValueError):
            FilteredComplex([((0,), 1), ((1,), 0), ((0,1), 0)])

    def test_single_vertex_and_circle(self):
        self.assertEqual(FilteredComplex([(('v',), 0)]).intervals(),
                         [(0, 0, None)])
        graph = FilteredComplex([((0,), 0), ((1,), 0), ((2,), 0),
                                 ((0,1), 1), ((1,2), 1), ((0,2), 1)])
        self.assertEqual(graph.betti(2), {0: 1, 1: 1})

    def test_exact_compatibility_bridge(self):
        c = triangle()
        self.assertEqual(c.page(1), c.e1_ranks())
        self.assertEqual(c.page(2), {(0, 0): 1})
        self.assertEqual(sum(x['rank'] for x in c.differential(1).values()), 3)


if __name__ == '__main__':
    unittest.main()
