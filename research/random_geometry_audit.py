"""Seeded adversarial-ish regression for finite geometric map filtrations.

Build small valid simplicial maps across several base triangulations and prime
fields, verify the entire spectral sequence, independent relative E1 ranks,
chain-map naturality, stabilization, persistence, and total homology.

Reproduce: python -m research.random_geometry_audit --trials 1000 --seed 20261106
A deterministic randomized smoke test, NOT mathematical formal verification.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import random

from simplicial_maps import SimplicialComplex, SimplicialMap


def audit(trials=1000, seed=20261106):
    if type(trials) is not int or trials<0 or type(seed) is not int:
        raise ValueError('trials must be a nonnegative integer and seed an integer')
    rng=random.Random(seed)
    bases=[SimplicialComplex([(0,1),(1,2),(0,2)]),
           SimplicialComplex([(0,1,2)]),
           SimplicialComplex([(0,1),(2,3)]),
           SimplicialComplex([(0,1,2),(0,1,3),(0,2,3),(1,2,3)]),
           SimplicialComplex([(0,),(1,),(2,)])]
    tally=Counter()
    digest=sha256()
    max_simplices=0
    for t in range(trials):
        base=rng.choice(bases)
        n=rng.randrange(1,8)
        vertex_map={f'x{i}':rng.choice(base.vertices) for i in range(n)}
        facets=[(v,) for v in vertex_map]
        vertices=list(vertex_map)
        for k in range(2,min(n,4)+1):
            for face in combinations(vertices,k):
                image=tuple(sorted({vertex_map[v] for v in face}))
                if image in base.simplex_set and rng.random()<0.28:
                    facets.append(face)
        projection=SimplicialMap(SimplicialComplex(facets),base,vertex_map)
        prime=rng.choice((2,3,5,7))
        digest.update(json.dumps(projection.source_data(prime),sort_keys=True,
                                 separators=(',',':')).encode()+b'\n')
        if not projection.verify(prime):
            raise AssertionError(f'geometry verification failed at trial {t}')
        if projection.chain_complex(prime).betti()!=projection.domain.chain_complex(prime)[0].betti():
            raise AssertionError(f'geometric total homology mismatch at trial {t}')
        tally[prime]+=1
        max_simplices=max(max_simplices,len(projection.domain.simplices))
    return {'schema':'random-geometric-map-audit-v1','seed':seed,'trials':trials,
            'verified':trials,'by_prime':{str(p):tally[p] for p in (2,3,5,7)},
            'max_domain_simplices':max_simplices,'input_sha256':digest.hexdigest()}


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--trials',type=int,default=1000)
    parser.add_argument('--seed',type=int,default=20261106)
    parser.add_argument('--output',help='optional deterministic summary JSON output')
    args=parser.parse_args(argv)
    receipt=audit(args.trials,args.seed)
    text=json.dumps(receipt,indent=2,sort_keys=True)+'\n'
    if args.output:
        Path(args.output).write_text(text,encoding='utf-8')
    print(text,end='')


if __name__=='__main__':
    main()
