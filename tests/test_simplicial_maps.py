"""Reusable geometric filtration, mapping-torus and Serre benchmark regressions."""
from collections import Counter
from bockstein import LiftedBockstein
from hashlib import sha256
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest

from simplicial_maps import (
    SimplicialComplex, SimplicialMap, _example, _packed,
    simplicial_product, mapping_torus_projection, monodromy_on_homology,
    circle_local_system_e2, replay_simplicial_map_certificate,
)


def cycle():
    return SimplicialComplex([('a','b'),('b','c'),('a','c')])


def sphere():
    return SimplicialComplex([tuple(v) for v in ('abc','abd','acd','bcd')])


def page_rank(chain,r):
    return sum(x['rank'] for x in chain.differential(r).values())


class ComplexConstructionTests(unittest.TestCase):
    def test_closure_and_counts(self):
        k=SimplicialComplex([('x','z','y')])
        self.assertEqual(k.f_vector(), {0:3,1:3,2:1})
        self.assertEqual(k.facets, (('x','y','z'),))
        self.assertEqual(SimplicialComplex([('a',),('a','b')]).facets,(('a','b'),))
        self.assertEqual(SimplicialComplex([]).f_vector(),{})

    def test_valid_mixed_vertices_and_nested_product_labels(self):
        k=SimplicialComplex([(1,'s',(2,'q'))])
        self.assertEqual(k.f_vector(),{0:3,1:3,2:1})
        m=SimplicialMap(k,SimplicialComplex([('t',)]),{v:'t' for v in k.vertices})
        self.assertEqual(set(m.births.values()),{0})
        self.assertTrue(m.verify(3))
        self.assertEqual(m.certificate(3),SimplicialMap.from_data(m.source_data(3)).certificate(3))

    def test_invalid_facets_and_vertices(self):
        for facets in [[()],[(1,1)],[(True,)],[(1.5,)],[([1],)],[(0,1),(1,0)]]:
            with self.subTest(facets=repr(facets)):
                with self.assertRaises(ValueError):SimplicialComplex(facets)

    def test_closed_subcomplex_validation(self):
        k=SimplicialComplex([(0,1,2)])
        with self.assertRaises(ValueError):k.subcomplex([(0,1)])
        with self.assertRaises(ValueError):k.subcomplex([(0,),(3,)])
        self.assertEqual(k.subcomplex([(0,),(1,)]).f_vector(),{0:2})

    def test_invalid_simplicial_vertex_maps(self):
        k=SimplicialComplex([('a','b')])
        b=SimplicialComplex([(1,),(2,)])
        with self.assertRaises(ValueError):SimplicialMap(k,b,{'a':1,'b':2})
        with self.assertRaises(ValueError):SimplicialMap(k,b,{'a':1})
        with self.assertRaises(ValueError):SimplicialMap(k,b,{'a':1,'b':3})
        with self.assertRaises(ValueError):SimplicialMap(k,b,{'a':1,'b':1,'extra':2})
        with self.assertRaises(ValueError):SimplicialMap.from_data({'domain_facets':[['a'],['b']],
            'base_facets':[[0]],'vertex_map':[['a',0],['a',0]]})

    def test_empty_domain_and_base(self):
        m=SimplicialMap(SimplicialComplex([]),SimplicialComplex([]),{})
        self.assertEqual(m.chain_complex().betti(),{})
        self.assertEqual(m.chain_complex().page(1),{})
        self.assertTrue(m.verify())
        self.assertEqual(m.preimage_skeleton(-1).f_vector(),{})
        replay_simplicial_map_certificate(json.loads(json.dumps(m.certificate())))


class GeometricMapTests(unittest.TestCase):
    def test_hopf_from_facets_reproduces_published_pages(self):
        m=_example('hopf')
        self.assertEqual(m.domain.f_vector(),{0:12,1:48,2:72,3:36})
        self.assertEqual(m.base.f_vector(),{0:4,1:6,2:4})
        self.assertEqual(Counter(m.births.values()),{0:24,1:72,2:72})
        self.assertEqual(m.domain.dimension,3)
        for p in (2,3,5,7):
            with self.subTest(prime=p):
                c=m.chain_complex(p)
                self.assertEqual(c.page(2),{(0,0):1,(0,1):1,(2,0):1,(2,1):1})
                self.assertEqual(c.differential(2)[(2,0)]['rank'],1)
                self.assertEqual(c.page(3),{(0,0):1,(2,1):1})
                self.assertEqual(c.betti(),{0:1,3:1})
                self.assertTrue(m.verify(p))

    def test_product_sphere_same_e2_as_hopf_but_d2_zero(self):
        m=simplicial_product(sphere(),cycle())
        h=_example('hopf')
        self.assertEqual(len(m.domain.simplices),len(h.domain.simplices))
        self.assertEqual(m.domain.dimension,h.domain.dimension)
        for p in (2,3,5):
            with self.subTest(prime=p):
                c=m.chain_complex(p)
                self.assertEqual(c.page(2),h.chain_complex(p).page(2))
                self.assertEqual(c.differential(2)[(2,0)]['rank'],0)
                self.assertEqual(c.betti(),{0:1,1:1,2:1,3:1})
                self.assertTrue(m.verify(p))
                for v in m.base.vertices:
                    self.assertEqual(m.fiber_over_vertex(v).chain_complex(p)[0].betti(),{0:1,1:1})

    def test_hopf_fibers_and_base_subcomplex_preimages(self):
        m=_example('hopf')
        for v in m.base.vertices:
            fiber=m.fiber_over_vertex(v)
            self.assertEqual(fiber.f_vector(),{0:3,1:3})
            self.assertEqual(fiber.chain_complex(3)[0].betti(),{0:1,1:1})
        self.assertEqual(m.preimage_skeleton(-1).f_vector(),{})
        self.assertEqual(len(m.preimage_skeleton(0).simplices),24)
        self.assertEqual(len(m.preimage_skeleton(1).simplices),96)
        self.assertEqual(len(m.preimage_skeleton(2).simplices),168)
        self.assertEqual(m.preimage_subcomplex(SimplicialComplex([('a',)])).simplices,
                         m.fiber_over_vertex('a').simplices)
        with self.assertRaises(ValueError):m.fiber_over_vertex('unknown')
        with self.assertRaises(ValueError):m.preimage_skeleton('zero')
        with self.assertRaises(ValueError):m.preimage_subcomplex(SimplicialComplex([('x',)]))

    def test_square_product_triangulation_and_naturality(self):
        segment=SimplicialComplex([(0,1)])
        m=simplicial_product(segment,segment)
        self.assertEqual(m.domain.f_vector(),{0:4,1:5,2:2})
        for p in (2,3,5):
            self.assertEqual(m.chain_complex(p).betti(),{0:1})
            self.assertTrue(m.chain_map(p).verify_naturality())
            self.assertTrue(m.verify(p))

    def test_degenerate_simplices_map_to_zero_and_signs(self):
        k=SimplicialComplex([('a','b','c')])
        base=SimplicialComplex([('u','v','w')])
        m=SimplicialMap(k,base,{'a':'v','b':'u','c':'w'})
        f=m.chain_map(5)
        source_lookup=m.domain.chain_complex(5,births=m.births)[1]
        self.assertEqual(sum(f.images[source_lookup[('a','b','c')]]),4)  # orientation sign -1
        self.assertTrue(f.verify_naturality())
        collapse=SimplicialMap(k,base,{'a':'u','b':'u','c':'v'})
        collapsed_map=collapse.chain_map(3)
        sc=collapse.domain.chain_complex(3,births=collapse.births)[1]
        self.assertFalse(any(collapsed_map.images[sc[('a','b','c')]]))
        self.assertFalse(any(collapsed_map.images[sc[('a','b')]]))
        self.assertTrue(collapse.verify(3))

    def test_nonsurjective_map_is_not_necessarily_serre(self):
        k=SimplicialComplex([('v',)])
        b=cycle()
        m=SimplicialMap(k,b,{'v':'a'})
        self.assertEqual(m.chain_complex(3).page(2),{(0,0):1})
        self.assertEqual(b.chain_complex(3)[0].betti(),{0:1,1:1})
        self.assertTrue(m.verify(3))

    def test_label_and_facets_permutation_does_not_change_pages(self):
        rng=random.Random(12084)
        original=_example('hopf')
        pages={r:original.chain_complex(3).page(r) for r in range(4)}
        ranks={r:page_rank(original.chain_complex(3),r) for r in range(4)}
        for repeat in range(3):
            facets=[list(face) for face in original.domain.facets]
            rng.shuffle(facets)
            for face in facets:rng.shuffle(face)
            bfacets=[list(face) for face in original.base.facets]
            rng.shuffle(bfacets)
            test=SimplicialMap(SimplicialComplex(facets),SimplicialComplex(bfacets),original.vertex_map)
            c=test.chain_complex(3)
            self.assertEqual({r:c.page(r) for r in range(4)},pages)
            self.assertEqual({r:page_rank(c,r) for r in range(4)},ranks)


class MonodromyTests(unittest.TestCase):
    def test_trivial_and_klein_tori(self):
        fiber=cycle()
        identity={x:x for x in fiber.vertices}
        reflect={'a':'b','b':'a','c':'c'}
        torus=mapping_torus_projection(fiber,identity)
        klein=mapping_torus_projection(fiber,reflect)
        self.assertEqual(torus.domain.f_vector(),{0:9,1:27,2:18})
        self.assertEqual(klein.domain.f_vector(),{0:9,1:27,2:18})
        for p in (2,3,5,7):
            with self.subTest(prime=p):
                self.assertEqual(torus.chain_complex(p).betti(),{0:1,1:2,2:1})
                self.assertEqual(klein.chain_complex(p).betti(),
                                 {0:1,1:2,2:1} if p==2 else {0:1,1:1})
                self.assertEqual(torus.chain_complex(p).page(2),
                                 circle_local_system_e2(fiber,identity,p))
                self.assertEqual(klein.chain_complex(p).page(2),
                                 circle_local_system_e2(fiber,reflect,p))
                self.assertEqual(monodromy_on_homology(fiber,reflect,p)[1],
                                 (((-1)%p,),))
                self.assertTrue(klein.verify(p))

    def test_transposition_of_discrete_two_point_fiber_is_connected_cover(self):
        fiber=SimplicialComplex([(0,),(1,)])
        swap={0:1,1:0}
        trivial=mapping_torus_projection(fiber,{0:0,1:1},segments=4)
        double=mapping_torus_projection(fiber,swap,segments=4)
        for p in (2,3,5):
            self.assertEqual(trivial.chain_complex(p).betti(),{0:2,1:2})
            self.assertEqual(double.chain_complex(p).betti(),{0:1,1:1})
            self.assertEqual(double.chain_complex(p).page(2),
                             circle_local_system_e2(fiber,swap,p))
            self.assertTrue(double.verify(p))

    def test_higher_dimensional_local_system_and_integral_torsion(self):
        untwisted=_example('sphere-torus')
        twisted=_example('twisted-sphere-torus')
        fiber=sphere()
        reflection={'a':'b','b':'a','c':'c','d':'d'}
        self.assertEqual(twisted.domain.f_vector(),{0:12,1:48,2:72,3:36})
        self.assertEqual(untwisted.domain.integral_chain_complex().homology()[2],
                         {'free_rank':1,'torsion':()})
        self.assertEqual(twisted.domain.integral_chain_complex().homology()[2],
                         {'free_rank':0,'torsion':(2,)})
        for p in (2,3,5):
            with self.subTest(prime=p):
                expected=({0:1,1:1,2:1,3:1} if p==2 else {0:1,1:1})
                self.assertEqual(twisted.chain_complex(p).betti(),expected)
                self.assertEqual(twisted.chain_complex(p).page(2),
                                 circle_local_system_e2(fiber,reflection,p))
                self.assertEqual(monodromy_on_homology(fiber,reflection,p)[2],
                                 (((-1)%p,),))
                self.assertTrue(twisted.verify(p))

    def test_klein_integral_torsion_and_first_bockstein(self):
        torus=_example('torus').domain.integral_chain_complex()
        klein=_example('klein').domain.integral_chain_complex()
        self.assertEqual(torus.homology()[1],{'free_rank':2,'torsion':()})
        self.assertEqual(klein.homology()[1],{'free_rank':1,'torsion':(2,)})
        self.assertEqual(klein.homology()[2],{'free_rank':0,'torsion':()})
        for model,expected in ((torus,0),(klein,1)):
            b=LiftedBockstein(model,prime=2)
            rank=sum(sum(any(x for x in col) for col in matrix)
                     for matrix in b.differential(1).values())
            self.assertEqual(rank,expected)
            self.assertTrue(b.verify())
        # In this fixture the nonzero map goes H_2(F2) -> H_1(F2).
        differential=LiftedBockstein(klein,prime=2).differential(1)
        self.assertTrue(any(any(x for x in col) for col in differential[2]))

    def test_mapping_torus_rejects_invalid_gluings(self):
        fiber=cycle()
        with self.assertRaises(ValueError):mapping_torus_projection(fiber,{'a':'a'},3)
        with self.assertRaises(ValueError):mapping_torus_projection(fiber,
            {'a':'a','b':'a','c':'c'},3)
        with self.assertRaises(ValueError):mapping_torus_projection(fiber,
            {'a':'a','b':'b','c':'c'},2)
        path=SimplicialComplex([(0,1),(1,2)])
        with self.assertRaises(ValueError):mapping_torus_projection(path,{0:1,1:0,2:2})


class ReceiptAndFuzzTests(unittest.TestCase):
    def test_human_json_roundtrip_cli(self):
        data={'schema':'simplicial-map-input-plain-v1','prime':3,
              'domain_facets':[['s','t']],'base_facets':[['a','b']],
              'vertex_map':[['s','a'],['t','b']]}
        m=SimplicialMap.from_data(data)
        self.assertTrue(m.verify(3))
        with tempfile.TemporaryDirectory() as tmp:
            inp=Path(tmp)/'input.json'
            receipt=Path(tmp)/'receipt.json'
            inp.write_text(json.dumps(data),encoding='utf-8')
            run=subprocess.run([sys.executable,'-m','simplicial_maps',str(inp),
                               '--certificate',str(receipt)],capture_output=True,text=True)
            self.assertEqual(run.returncode,0,run.stdout+run.stderr)
            self.assertIn('F_3',run.stdout)
            rerun=subprocess.run([sys.executable,'-m','simplicial_maps',
                                  '--replay',str(receipt)],capture_output=True,text=True)
            self.assertEqual(rerun.returncode,0,rerun.stdout+rerun.stderr)
            self.assertIn('PASS',rerun.stdout)

    def test_certificates_and_digest_tamper_detection(self):
        for name in ('torus','klein','hopf'):
            m=_example(name)
            receipt=json.loads(json.dumps(m.certificate(3)))
            self.assertTrue(replay_simplicial_map_certificate(receipt))
            changed=json.loads(json.dumps(receipt))
            changed['pages'][1]['groups'][0][2]+=1
            with self.assertRaises(ValueError):replay_simplicial_map_certificate(changed)
            # A digest is an integrity check, not a signature. Even if someone
            # recomputes it after tampering, recomputation must catch changed data.
            changed['report_sha256']=sha256(_packed({k:v for k,v in changed.items()
                                                       if k!='report_sha256'}).encode()).hexdigest()
            with self.assertRaises(ValueError):replay_simplicial_map_certificate(changed)

    def test_random_valid_simplicial_maps(self):
        rng=random.Random(100932)
        bases=[SimplicialComplex([(0,1),(1,2),(0,2)]),
               SimplicialComplex([(0,1,2)]),
               SimplicialComplex([(0,1),(2,3)]),
               SimplicialComplex([(0,), (1,), (2,)])]
        for trial in range(100):
            base=rng.choice(bases)
            n=rng.randrange(1,7)
            vertex_map={f'v{i}':rng.choice(base.vertices) for i in range(n)}
            vertices=list(vertex_map)
            faces=[(v,) for v in vertices]
            for k in range(2,min(n,4)+1):
                for sub in __import__('itertools').combinations(vertices,k):
                    image=tuple(sorted({vertex_map[v] for v in sub}))
                    if image in base.simplex_set and rng.random()<0.3:
                        faces.append(sub)
            m=SimplicialMap(SimplicialComplex(faces),base,vertex_map)
            p=rng.choice((2,3,5))
            with self.subTest(trial=trial,prime=p):
                c=m.chain_complex(p)
                self.assertTrue(m.verify(p))
                self.assertEqual(m.relative_e1_page(p),c.page(1))
                self.assertEqual(Counter(b for _,b in m.induced_filtration()),Counter(m.births.values()))
                self.assertEqual(c.betti(),m.domain.chain_complex(p)[0].betti())

    def test_seeded_audit_receipt_is_reproducible(self):
        from research.random_geometry_audit import audit
        report=audit(15,20261106)
        self.assertEqual(report,audit(15,20261106))
        self.assertEqual(report['verified'],15)
        self.assertEqual(sum(report['by_prime'].values()),15)
        self.assertEqual(len(report['input_sha256']),64)
        with self.assertRaises(ValueError):audit(-1)

    def test_input_rejections(self):
        with self.assertRaises(ValueError):SimplicialMap.from_data({'schema':'foo'})
        with self.assertRaises(ValueError):SimplicialMap.from_data({'schema':'simplicial-map-input-v1',
            'domain_facets':[[['bool',True]]],'base_facets':[],'vertex_map':[]})
        with self.assertRaises(ValueError):_example('hopf').certificate(4)


if __name__=='__main__':unittest.main()
