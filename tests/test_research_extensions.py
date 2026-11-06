import random
import unittest
from bicomplex import DoubleComplex, tensor_double_complex
from bockstein import BocksteinSpectralSequence, LiftedBockstein, replay_lifted_bockstein_certificate
from integral import IntegralChainComplex
from graded_couple import GradedExactCouple
from spectral_sequence import FilteredChainComplex


def staircase(length, prime):
    """Vertical cancellation chain with first nonzero d_length."""
    if length < 2:
        raise ValueError('minimum length is two')
    data=[('x',length,0,{}, {'a1':1}),
          ('y',0,length-1,{}, {})]
    for j in range(1,length):
        p=length-j
        q=j-1
        data.append((f'a{j}',p,q,{}, {}))
        target='y' if j==length-1 else f'a{j+1}'
        data.append((f'b{j}',p,q+1,{f'a{j}':-1},{target:1}))
    return DoubleComplex(data,prime)


class DoubleComplexTests(unittest.TestCase):
    def test_hidden_d2_and_d3_for_multiple_primes(self):
        for prime in (2,3,5):
            for r in (2,3,4,5):
                b=staircase(r,prime)
                self.assertTrue(b.verify())
                col=b.column_spectral_sequence()
                self.assertEqual(col.page(1),{(0,r-1):1,(r,0):1})
                self.assertTrue(all(sum(e['rank'] for e in col.differential(i).values())==0
                                    for i in range(1,r)))
                self.assertEqual(sum(e['rank'] for e in col.differential(r).values()),1)
                self.assertEqual(col.page(r+1),{})
                self.assertEqual(b.row_spectral_sequence().page(1),{})

    def test_vertical_homology_and_both_abutments(self):
        b=DoubleComplex([
            ('u',1,1,{}, {}), ('v',0,2,{}, {}),
            ('w',1,0,{}, {}), ('z',0,1,{}, {}),
        ],7)
        self.assertTrue(b.verify())
        self.assertEqual(b.column_spectral_sequence().betti(),{1:2,2:2})
        self.assertEqual(b.row_spectral_sequence().betti(),{1:2,2:2})

    def test_rejects_wrong_bidegree_and_noncommuting_square(self):
        with self.assertRaises(ValueError):
            DoubleComplex([('x',0,0,{}, {'missing':1})])
        with self.assertRaises(ValueError):
            DoubleComplex([
                ('x',1,1,{'a':1},{'b':1}),
                ('a',1,0,{}, {'t':1}),
                ('b',0,1,{'t':1},{}),
                ('t',0,0,{}, {})],3)
        with self.assertRaises(ValueError):
            DoubleComplex([('x',0,0,{'x':1},{})],2)

    def test_tensor_product_kunneth_and_signed_differential(self):
        for prime in (2,3,5,7):
            left=FilteredChainComplex([('a',0,0,{}),('b',0,1,{}),
                                       ('x',1,2,{'a':1,'b':-1}),
                                       ('y',1,3,{'a':1,'b':-1})],prime)
            right=FilteredChainComplex([('u',0,0,{}),('v',0,1,{}),
                                        ('z',1,2,{'u':1,'v':-1})],prime)
            bic=tensor_double_complex(left,right)
            self.assertTrue(bic.verify())
            # H(left): H0=1,H1=1; H(right): H0=1, H1=0.
            self.assertEqual(bic.column_spectral_sequence().betti(),{0:1,1:1})
            self.assertEqual(bic.column_spectral_sequence().page(2),
                             {(0,0):1,(1,0):1})
            self.assertEqual(bic.row_spectral_sequence().page(2),
                             {(0,0):1,(0,1):1})

    def test_tensor_product_rejects_mixed_characteristic(self):
        a=FilteredChainComplex([('a',0,0,{})],2)
        b=FilteredChainComplex([('b',0,0,{})],3)
        with self.assertRaises(ValueError): tensor_double_complex(a,b)

    def test_empty_and_nonzero_anticommutator(self):
        self.assertTrue(DoubleComplex([],2).verify())
        b=DoubleComplex([
            ('x',1,1,{'a':1},{'b':1}),
            ('a',1,0,{}, {'t':1}),
            ('b',0,1,{'t':-1},{}),
            ('t',0,0,{}, {})],3)
        self.assertTrue(b.verify())


class BocksteinTests(unittest.TestCase):
    def test_p_power_torsion_detected_on_exact_page(self):
        for p in (2,3,5):
            for exponent in (1,2,3,4):
                c=IntegralChainComplex([1,1],{1:[[p**exponent]]})
                seq=BocksteinSpectralSequence(c,p)
                self.assertTrue(seq.verify())
                self.assertEqual(seq.page(1),{0:1,1:1})
                self.assertEqual(seq.page(exponent),{0:1,1:1})
                self.assertEqual(seq.differential_ranks(exponent),{1:1})
                self.assertEqual(seq.page(exponent+1),{})
                self.assertEqual(seq.torsion_from_differentials()[0],(p**exponent,))

    def test_away_from_prime_and_free_part(self):
        c=IntegralChainComplex([2,2],{1:[[6,0],[0,0]]})
        for p in (2,3,5):
            s=BocksteinSpectralSequence(c,p)
            self.assertTrue(s.verify())
            self.assertEqual(s.page(100),{0:1,1:1})
            self.assertEqual(s.differential_ranks(1),{1:1} if p in (2,3) else {})

    def test_mixed_torsion_across_degrees(self):
        c=IntegralChainComplex([3,3,2], {1:[[4,0,0],[0,9,0],[0,0,0]],
                                         2:[[0,0],[0,0],[8,25]]})
        for p in (2,3,5):
            s=BocksteinSpectralSequence(c,p)
            self.assertTrue(s.verify())
            for r in range(1,s.max_exponent+2):
                self.assertEqual(s.page(r+1),
                    {n:x for n,x in ((n,s.page(r).get(n,0)-s.differential_ranks(r).get(n,0)
                                      -s.differential_ranks(r).get(n+1,0))
                                      for n in range(3)) if x})

    def test_random_diagonal_presentations(self):
        random.seed(1106)
        for _ in range(200):
            widths=random.randrange(1,6)
            factors=[random.choice([0,1,2,3,4,6,8,9,16,25,27,30,50,125])
                     for _ in range(widths)]
            c=IntegralChainComplex([widths,widths],
                                   {1:[[factors[i] if i==j else 0 for j in range(widths)]
                                       for i in range(widths)]})
            for p in (2,3,5):
                seq=BocksteinSpectralSequence(c,p)
                self.assertTrue(seq.verify())

    def test_lifted_higher_bockstein_differentials(self):
        for p in (2,3,5):
            for r in (1,2,3,4):
                c=IntegralChainComplex([1,1],{1:[[p**r]]})
                lift=LiftedBockstein(c,p)
                self.assertTrue(lift.verify())
                for i in range(1,r):
                    self.assertEqual(lift.differential(i)[1],((0,),))
                self.assertEqual(lift.differential(r)[1],((1,),))
                self.assertEqual(lift.page(r+1),{})

    def test_lifts_on_non_diagonal_presentations(self):
        rng=random.Random(24246)
        for i in range(180):
            dim=rng.randrange(1,5)
            factor=[rng.choice([0,1,2,3,4,6,8,9,16,25,27,30,125])
                    for _ in range(dim)]
            a=[[factor[k] if k==j else 0 for j in range(dim)] for k in range(dim)]
            # Unimodular row/column changes preserve the homology but can
            # force nontrivial integer lift corrections and cancellations.
            for _ in range(4):
                x,y=rng.sample(range(dim),2) if dim>1 else (0,0)
                if x==y: continue
                v=rng.choice([-2,-1,1,2])
                if rng.choice([False,True]):
                    for j in range(dim): a[x][j]+=v*a[y][j]
                else:
                    for k in range(dim): a[k][x]+=v*a[k][y]
            c=IntegralChainComplex([dim,dim],{1:a})
            for p in (2,3,5):
                lift=LiftedBockstein(c,p)
                self.assertTrue(lift.verify())
                self.assertEqual(lift.page(1),lift.fingerprint.page(1))
                self.assertEqual(lift.page(100),lift.fingerprint.page(100))

    def test_random_multidegree_chain_complexes_and_unimodular_changes(self):
        rng=random.Random(909273)
        for trial in range(100):
            generators=[[] for _ in range(4)]
            pairings=[]
            for _ in range(rng.randrange(1,8)):
                degree=rng.choice((1,2,3))
                source=len(generators[degree]); generators[degree].append(('s',trial))
                target=len(generators[degree-1]); generators[degree-1].append(('t',trial))
                coefficient=rng.choice((1,2,3,4,6,8,9,16,18,25,27,30,75))
                pairings.append((degree,source,target,coefficient))
            for degree in range(4):
                for _ in range(rng.randrange(3)):
                    generators[degree].append(('free',trial))
            ranks=list(map(len,generators))
            boundaries={n:[[0]*ranks[n] for _ in range(ranks[n-1])]
                        for n in range(1,4)}
            for degree,source,target,coefficient in pairings:
                boundaries[degree][target][source]=coefficient
            # Change integer bases, and propagate the inverse transformation
            # through adjacent maps so d**2 stays exactly zero.
            for degree in range(4):
                for _ in range(7):
                    if ranks[degree]<2: continue
                    a,b=rng.sample(range(ranks[degree]),2)
                    factor=rng.choice((-2,-1,1,2))
                    if degree>0:
                        for row in boundaries[degree]:
                            row[a]+=factor*row[b]
                    if degree<3:
                        upper=boundaries[degree+1]
                        for j in range(ranks[degree+1]):
                            upper[b][j]-=factor*upper[a][j]
            chain=IntegralChainComplex(ranks,boundaries)
            for prime in (2,3,5):
                lifted=LiftedBockstein(chain,prime)
                self.assertTrue(lifted.verify())
                fp=lifted.fingerprint
                self.assertEqual(lifted.page(fp.max_exponent+1),fp.page(fp.max_exponent+1))

    def test_lift_certificate_and_tamper_detection(self):
        c=IntegralChainComplex([2,2],{1:[[9,3],[0,0]]})
        lift=LiftedBockstein(c,3)
        report=lift.certificate()
        self.assertTrue(replay_lifted_bockstein_certificate(report))
        altered=dict(report)
        altered['sha256']='0'*64
        with self.assertRaises(ValueError):
            replay_lifted_bockstein_certificate(altered)
        altered=dict(report)
        altered['pages']=list(report['pages'])
        altered['pages'][0]=dict(altered['pages'][0])
        altered['pages'][0]['r']=23
        import hashlib, json
        altered['sha256']=hashlib.sha256(json.dumps({k:v for k,v in altered.items() if k!='sha256'},
                           sort_keys=True,separators=(',',':')).encode()).hexdigest()
        with self.assertRaises(ValueError):
            replay_lifted_bockstein_certificate(altered)

    def test_two_adjacent_boundary_maps_with_nontrivial_kernel(self):
        # C_2 -> C_1 -> C_0: incoming map hits the kernel of d1.
        c=IntegralChainComplex([1,3,2],
                 {1:[[0,6,0]],2:[[4,0],[0,0],[0,9]]})
        for p in (2,3,5):
            lift=LiftedBockstein(c,p)
            self.assertTrue(lift.verify())
            self.assertEqual(lift.page(5),lift.fingerprint.page(5))

    def test_invalid_prime_and_page(self):
        c=IntegralChainComplex([1],{})
        with self.assertRaises(ValueError):
            BocksteinSpectralSequence(c,4)
        s=BocksteinSpectralSequence(c,2)
        for bad in (0,-1,True):
            with self.assertRaises(ValueError):
                s.page(bad)


class GradedExactCoupleTests(unittest.TestCase):
    def test_hidden_second_differential_preserves_grading(self):
        for p in (2,3,5):
            c=GradedExactCouple(p, (0,0), (0,1),
                                ((0,0),(1,0)),
                                ((0,0),(1,0)),
                                ((0,0),(1,0)))
            self.assertEqual(c.differential_by_degree()[1]['rank'],0)
            self.assertEqual(c.derive().differential_by_degree()[1]['rank'],1)
            self.assertEqual(c.derive().derive().page_dimensions(),{})
            self.assertEqual([x['d_ranks'] for x in c.derived_pages()], [{},{1:1},{}])

    def test_first_differential_and_degrees(self):
        c=GradedExactCouple(3,(0,),(0,1),
                            ((0,),),((1,0),),((0,),(1,)))
        self.assertEqual(c.differential_by_degree()[1]['rank'],1)
        self.assertEqual(c.derive().page_dimensions(),{})

    def test_mixed_degree_stable_blocks(self):
        c=GradedExactCouple(2,(0,0,5,5),(0,1,5,6),
                            ((0,0,0,0),(1,0,0,0),(0,0,0,0),(0,0,1,0)),
                            ((0,0,0,0),(1,0,0,0),(0,0,0,0),(0,0,1,0)),
                            ((0,0,0,0),(1,0,0,0),(0,0,0,0),(0,0,1,0)))
        self.assertEqual(c.derive().page_dimensions(),{0:1,1:1,5:1,6:1})
        self.assertEqual(c.derive().differential_by_degree()[1]['rank'],1)
        self.assertEqual(c.derive().differential_by_degree()[6]['rank'],1)

    def test_rejects_wrong_degrees(self):
        with self.assertRaises(ValueError):
            GradedExactCouple(3,(0,0),(0,0),
                              ((0,0),(1,0)),((0,0),(1,0)),((0,0),(1,0)))
        with self.assertRaises(ValueError):
            GradedExactCouple(3,(0,False),(0,1),
                              ((0,0),(1,0)),((0,0),(1,0)),((0,0),(1,0)))


if __name__=='__main__':
    unittest.main()
