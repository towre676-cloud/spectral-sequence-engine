"""Exact examples and deterministic generated regression checks (stdlib only)."""
import copy
import itertools
import json
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from spectral_sequence import (FilteredChainComplex, _Span, _example,
                               _kernel, _rank, replay_certificate)


class LinearAlgebraTests(unittest.TestCase):
    def test_field_reduction_and_coordinates(self):
        for prime in (2, 3, 5, 7):
            s = _Span(3, prime)
            self.assertTrue(s.add((2, 1, 0)))
            self.assertTrue(s.add((0, 1, 1)))
            self.assertFalse(s.add((4, 2, 0)))
            self.assertEqual(s.coordinates((2, 2, 1)), (1, 1))
            self.assertEqual(_rank([(2, 1), (4, 2)], prime, 2), 1)
            self.assertEqual(len(_kernel([(1, 0), (1, 0)], prime, 2)), 1)
            with self.assertRaises(ArithmeticError):
                s.coordinates((1, 0, 0))

    def test_empty_linear_algebra(self):
        for prime in (2, 3):
            s = _Span(0, prime)
            self.assertFalse(s.add(()))
            self.assertEqual(s.coordinates(()), ())
            self.assertEqual(_kernel([(), ()], prime, 0), ((1, 0), (0, 1)))


class SpectralSequenceTests(unittest.TestCase):
    def test_triangle_matches_legacy_persistence_and_e1(self):
        from filtered_complex import FilteredComplex
        simplices = [((0,), 0), ((1,), 0), ((2,), 0),
                     ((0,1), 1), ((0,2), 1), ((1,2), 1),
                     ((0,1,2), 2)]
        old = FilteredComplex(simplices)
        new = FilteredChainComplex.from_simplicial(simplices)
        self.assertEqual(old.e1_ranks(), new.page(1))
        self.assertEqual(old.intervals(), new.intervals())
        self.assertEqual(new.page(1), {(0,0): 3, (1,0): 3, (2,0): 1})
        self.assertEqual(sum(x['rank'] for x in new.differential(1).values()), 3)
        self.assertEqual(new.page(2), {(0,0): 1})
        self.assertEqual(new.betti(), {0: 1})
        self.assertEqual(new.page(3), new.page(4))
        self.assertTrue(new.verify())

    def test_same_filtration_d0(self):
        for prime in (2,3,5):
            c = FilteredChainComplex([
                ('a',0,0,{}), ('b',1,0,{'a':1})], prime)
            self.assertEqual(sum(d['rank'] for d in c.differential(0).values()), 1)
            self.assertEqual(c.page(1), {})
            self.assertEqual(c.intervals(), [])
            self.assertEqual(c.intervals(True), [(0, 0, 0)])
            self.assertTrue(c.verify())

    def test_higher_differential_d2_d3(self):
        for prime in (2, 3, 5, 7):
            for length in (2, 3):
                c = FilteredChainComplex([
                    ('a',0,0,{}), ('b',1,length,{'a':1})], prime)
                for r in range(1, length):
                    self.assertEqual(sum(v['rank'] for v in c.differential(r).values()), 0)
                self.assertEqual(sum(v['rank'] for v in c.differential(length).values()), 1)
                self.assertEqual(c.page(length+1), {})
                self.assertTrue(c.verify())

    def test_hidden_d2_after_lifting(self):
        # At filtration 1, d_0(u)=b; d_1(v) is zero in E_1.
        # Nevertheless v-u has boundary -a and supports a nonzero d_2.
        for prime in (2,3,5):
            c = FilteredChainComplex([
                ('a',0,0,{}), ('b',0,1,{}),
                ('u',1,1,{'a':1,'b':1}), ('v',1,2,{'b':1})], prime)
            self.assertEqual(c.page(1), {(0,0):1,(2,-1):1})
            self.assertEqual(c.differential(1)[(2,-1)]['rank'], 0)
            self.assertEqual(c.differential(2)[(2,-1)]['rank'], 1)
            self.assertTrue(c.verify())

    def test_two_dimensional_simplices_odd_characteristic(self):
        simplices = [(x,0) for x in [(0,), (1,), (2,),
                    (0,1), (0,2), (1,2), (0,1,2)]]
        for prime in (2,3,5,7):
            c = FilteredChainComplex.from_simplicial(simplices, prime)
            self.assertEqual(c.betti(), {0:1})
            self.assertEqual(sum(x['rank'] for x in c.differential(0).values()), 3)
            self.assertEqual(c.page(1), {(0,0):1})
            self.assertTrue(c.verify())

    def test_disconnected_and_empty(self):
        c = FilteredChainComplex([],3)
        self.assertEqual(c.page(1), {})
        self.assertEqual(c.betti(), {})
        self.assertTrue(c.verify())
        c = FilteredChainComplex([('a',0,0,{}), ('b',0,2,{})],5)
        self.assertEqual(c.page(1), {(0,0):1,(2,-2):1})
        self.assertEqual(c.betti(), {0:2})
        self.assertTrue(c.verify())

    def test_input_validation(self):
        examples = [
            [('a',0,0,{}),('a',1,0,{})],
            [('a',0,0,{'b':1})],
            [('a',0,0,{}),('b',1,0,{'a':1}),('c',2,0,{'b':1})],
            [('a',0,1,{}),('b',1,0,{'a':1})],
            [('a',0,-1,{})],
            [('a',0,0,{'a':1})],
            [('a',False,0,{})],
            [('a',0,0,{'ghost':1})],
        ]
        for generators in examples:
            with self.subTest(generators=generators):
                with self.assertRaises(ValueError):
                    FilteredChainComplex(generators)
        for prime in (0,1,4,8,9,-2,False,2.5):
            with self.assertRaises(ValueError):
                FilteredChainComplex([],prime)
        with self.assertRaises(ValueError):
            FilteredChainComplex.from_simplicial([((0,1),0)])

    def test_deterministic_replay_and_tamper(self):
        c = _example('triangle')
        receipt = c.certificate()
        self.assertTrue(replay_certificate(receipt))
        self.assertEqual(receipt, c.certificate())
        changed = copy.deepcopy(receipt)
        changed['pages'][1]['groups'][0]['dimension'] += 1
        with self.assertRaises(ValueError):
            replay_certificate(changed)
        changed = copy.deepcopy(receipt)
        changed['schema'] = 'unknown'
        with self.assertRaises(ValueError):
            replay_certificate(changed)
        with tempfile.TemporaryDirectory() as tmp:
            f = Path(tmp)/'proof.json'
            f.write_text(json.dumps(receipt))
            result = subprocess.run([sys.executable, 'spectral_sequence.py',
                                      '--replay', str(f)], capture_output=True,
                                    text=True, check=True)
            self.assertIn('PASS', result.stdout)

    def test_cli_and_file_schema(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            f = root/'complex.json'; out = root/'receipt.json'
            f.write_text(json.dumps({'prime': 3, 'generators': [
                ['x', 0, 0, {}], ['y', 1, 2, {'x': 1}]]}))
            result = subprocess.run([sys.executable, 'spectral_sequence.py',
                                     str(f), '--certificate', str(out)],
                                    capture_output=True, text=True, check=True)
            self.assertIn('rank(d_2)=1', result.stdout)
            self.assertTrue(replay_certificate(json.loads(out.read_text())))

    def test_seeded_random_simplicial_filtrations(self):
        # Differential/page reduction and unrelated barcode reduction must agree.
        for seed in range(80):
            rng = random.Random(seed)
            n = rng.randrange(1, 7)
            births = {(v,): rng.randrange(4) for v in range(n)}
            for d in (1, 2, 3):
                for simplex in itertools.combinations(range(n), d+1):
                    faces = [simplex[:j]+simplex[j+1:] for j in range(d+1)]
                    if all(f in births for f in faces) and rng.randrange(3) == 0:
                        births[simplex] = max(births[f] for f in faces)+rng.randrange(3)
            for prime in (2, 3, 5):
                with self.subTest(seed=seed, prime=prime):
                    c = FilteredChainComplex.from_simplicial(births.items(), prime)
                    self.assertTrue(c.verify())
                    self.assertEqual(c.certificate()['stable_at'], c.spread+1)


if __name__ == '__main__':
    unittest.main()
