import copy
from itertools import combinations
import json
from math import gcd
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from integral import IntegralChainComplex, integral_mapping_cone
from integral_receipt_checker import verify_class_receipt


def determinant(matrix):
    if not matrix:
        return 1
    return sum(
        (-1) ** j * value * determinant([row[:j] + row[j + 1:] for row in matrix[1:]])
        for j, value in enumerate(matrix[0])
    )


def divisors(matrix):
    """Nonzero determinantal divisors, including the empty minor."""
    result = [1]
    for k in range(1, min(len(matrix), len(matrix[0])) + 1):
        divisor = gcd(*(
            determinant([[matrix[i][j] for j in cols] for i in rows])
            for rows in combinations(range(len(matrix)), k)
            for cols in combinations(range(len(matrix[0])), k)
        ))
        if not divisor:
            break
        result.append(divisor)
    return result


class IndependentReceiptTests(unittest.TestCase):
    def check_both(self, complex_, degree, chain, expected):
        receipt = json.loads(json.dumps(complex_.class_certificate(degree, chain)))
        status, order = expected
        self.assertEqual((receipt["status"], receipt.get("order")), (status, order))
        with patch("integral._diagonalize", side_effect=AssertionError("reduction called")):
            self.assertTrue(complex_.verify_class_certificate(degree, chain, receipt))
            self.assertTrue(verify_class_receipt(
                complex_.ranks, complex_.boundaries, degree, chain, receipt,
            ))
        return receipt

    def test_cone_and_empty_groups(self):
        circle = IntegralChainComplex([0, 1], {})
        cone = integral_mapping_cone(circle, circle, {1: [[4]]})
        for z, expected in [(1, ("torsion", 4)), (2, ("torsion", 2)), (4, ("boundary", None))]:
            self.check_both(cone, 1, [z], expected)
        c = IntegralChainComplex([0, 0, 3, 0, 2], {})
        for degree, chain, expected in [
            (1, [], ("boundary", None)), (2, [1, 0, 0], ("infinite", None)),
            (4, [0, 0], ("boundary", None)), (4, [1, 2], ("infinite", None)),
        ]:
            receipt = self.check_both(c, degree, chain, expected)
            self.assertTrue(verify_class_receipt(c.ranks, {}, degree, chain, receipt))

    def test_nondiagonal_complex_against_determinantal_divisors(self):
        rng = random.Random(20261115)
        # Integral coordinates on ker([2, 4, 0, 0, 0]).
        for factors in [(2, 4, 12, 0), (4, 9, 0, 0), (1, 2, 8, 16), (0, 0, 0, 0)]:
            boundary = [[a if i == j else 0 for j in range(4)] for i, a in enumerate(factors)]
            for _ in range(8):
                i, j = rng.sample(range(4), 2)
                coefficient = rng.choice([-2, -1, 1, 2])
                boundary[i] = [x + coefficient * y for x, y in zip(boundary[i], boundary[j])]
                for row in boundary:
                    row[j] += coefficient * row[i]
            incoming = [[-2 * value for value in boundary[0]], *boundary]
            c = IntegralChainComplex([1, 5, 4], {1: [[2, 4, 0, 0, 0]], 2: incoming})
            d = divisors(boundary)
            smith = [b // a for a, b in zip(d, d[1:])]
            self.assertEqual(
                c.homology()[1],
                {"free_rank": 5 - len(d), "torsion": tuple(k for k in smith if k > 1)},
            )
            self.check_both(c, 1, [1, 0, 0, 0, 0], ("noncycle", None))
            self.check_both(c, 1, [row[0] for row in incoming], ("boundary", None))
            for _ in range(24):
                z = [rng.randint(-9, 9) for _ in range(4)]
                augmented = divisors([row + [value] for row, value in zip(boundary, z)])
                if len(augmented) > len(d):
                    expected = ("infinite", None)
                else:
                    order = d[-1] // augmented[-1]
                    expected = ("boundary", None) if order == 1 else ("torsion", order)
                self.check_both(c, 1, [-2 * z[0], *z], expected)

    def test_hidden_extension_contrasting_total_orders(self):
        # Both filtrations give Z/2 at (0,1) and (1,0) on E1.
        for coupling, order, factors in [(0, 2, (2, 2)), (1, 4, (4,))]:
            c = IntegralChainComplex([0, 2, 2], {2: [[2, coupling], [0, 2]]})
            self.assertEqual(c.homology()[1]["torsion"], factors)
            self.check_both(c, 1, [0, 1], ("torsion", order))

    def test_torsion_forgery_and_malformed_fields(self):
        c = IntegralChainComplex([0, 2, 2], {2: [[4, 0], [0, 9]]})
        z = [2, 1]
        receipt = self.check_both(c, 1, z, ("torsion", 18))
        for changes in [
            {"kind": "integral-homology-class/2"}, {"status": []}, {"order": True},
            {"order": 18.0}, {"lift": [9, 3]}, {"lift": [True, 2]}, {"lift": [9]},
            {"lift": [9, 2, 0]}, {"order": 36, "lift": [18, 4]}, {"witnesses": []},
            {"witnesses": None}, {"witnesses": [False]}, {"extra": 0},
            {"input_sha256": "0" * 64},
        ]:
            with self.subTest(changes=changes):
                self.assertFalse(verify_class_receipt(c.ranks, c.boundaries, 1, z, receipt | changes))
        for changes in [
            {"modulus": 0}, {"modulus": -4}, {"modulus": 1}, {"modulus": True},
            {"modulus": 2}, {"modulus": 8}, {"functional": [1, 1]},
            {"functional": [True, 0]}, {"functional": [1]}, {"unknown": 1},
        ]:
            with self.subTest(changes=changes):
                forged = copy.deepcopy(receipt)
                forged["witnesses"][0].update(changes)
                self.assertFalse(verify_class_receipt(c.ranks, c.boundaries, 1, z, forged))
                self.assertFalse(c.verify_class_certificate(1, z, forged))
        # Witnesses need not be canonical. Twice the first row works modulo eight.
        alternative = copy.deepcopy(receipt)
        alternative["witnesses"][0] = {"modulus": 8, "functional": [2, 0]}
        self.assertTrue(verify_class_receipt(c.ranks, c.boundaries, 1, z, alternative))
        self.assertTrue(c.verify_class_certificate(1, z, alternative))
        self.assertFalse(verify_class_receipt(c.ranks, c.boundaries, 1, [6, 1], receipt))
        self.assertFalse(verify_class_receipt(c.ranks, c.boundaries, 2, z, receipt))
        self.assertFalse(verify_class_receipt(c.ranks, {2: [[8, 0], [0, 9]]}, 1, z, receipt))

    def test_other_forged_statuses_and_invalid_complex(self):
        c = IntegralChainComplex([1, 3, 1], {1: [[2, 4, 0]], 2: [[-8], [4], [0]]})
        for z, changes in [
            ([1, 0, 0], {"value": 4}),
            ([1, 0, 0], {"row": -1}),
            ([-8, 4, 0], {"lift": [0]}),
            ([-4, 2, 1], {"functional": [1, 0, 0]}),
            ([-4, 2, 1], {"functional": [2, 4, 0]}),
        ]:
            forged = c.class_certificate(1, z) | changes
            self.assertFalse(verify_class_receipt(c.ranks, c.boundaries, 1, z, forged))
            self.assertFalse(c.verify_class_certificate(1, z, forged))
        for ranks, boundaries, degree, chain in [
            ([1, 1, 1], {1: [[1]], 2: [[1]]}, 0, [1]),
            ([1, 1], {1: [[]]}, 0, [1]),
            ([True], {}, 0, [1]),
            ([1, 1], {True: [[0]]}, 0, [1]),
            ([1], {}, True, [1]),
            ([1], {}, 0, [True]),
        ]:
            with self.assertRaises(ValueError):
                verify_class_receipt(ranks, boundaries, degree, chain, {})

    def test_standalone_cli_and_strict_json(self):
        checker = Path(__file__).resolve().parents[1] / "integral_receipt_checker.py"
        c = IntegralChainComplex([0, 1, 1], {2: [[4]]})
        valid = json.dumps(c.class_certificate(1, [2]))
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            (directory / checker.name).write_bytes(checker.read_bytes())
            (directory / "input.json").write_text(
                '{"ranks":[0,1,1],"boundaries":{"2":[[4]]},"degree":1,"chain":[2]}',
                encoding="utf-8",
            )
            for serialized, expected_code in [
                (valid, 0), (valid.replace('"order": 2', '"order": 4'), 1),
                ('{"status":"torsion",' + valid[1:], 2),
                (valid.replace('"order": 2', '"order": 2.0'), 2),
                (valid.replace('"order": 2', '"order": NaN'), 2),
            ]:
                with self.subTest(code=expected_code, serialized=serialized):
                    (directory / "receipt.json").write_text(serialized, encoding="utf-8")
                    result = subprocess.run(
                        [sys.executable, "-I", checker.name, "input.json", "receipt.json"],
                        cwd=directory, capture_output=True, text=True,
                    )
                    self.assertEqual(result.returncode, expected_code, result.stderr)
