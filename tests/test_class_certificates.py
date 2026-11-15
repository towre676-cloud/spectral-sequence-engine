import copy
import json
import unittest
from unittest.mock import patch

from bockstein import LiftedBockstein
from integral import IntegralChainComplex, integral_mapping_cone
from spectral_sequence import FilteredChainComplex


class IntegralClassReceiptTests(unittest.TestCase):
    def test_filtered_degree_four_mapping_cone(self):
        circle = IntegralChainComplex([0, 1], {})
        cone = integral_mapping_cone(circle, circle, {1: [[4]]})
        self.assertEqual(cone.homology()[1], {"free_rank": 0, "torsion": (4,)})
        expected = {1: 4, 2: 2, 3: 4, 4: 1}
        for value, order in expected.items():
            with self.subTest(value=value):
                receipt = json.loads(json.dumps(cone.class_certificate(1, [value])))
                self.assertEqual(receipt.get("order", 1), order)
                with patch("integral._diagonalize", side_effect=AssertionError):
                    self.assertTrue(cone.verify_class_certificate(1, [value], receipt))
        for prime in (2, 3, 5, 7):
            field_cone = FilteredChainComplex([
                ("a", 1, 0, {}), ("b", 2, 1, {"a": 4})], prime=prime)
            self.assertEqual(field_cone.page(1), {(0, 1): 1, (1, 1): 1})
            self.assertEqual(field_cone.page(2), field_cone.page(1) if prime == 2 else {})
            # At p=2 the cycle 2a is zero as a vector; at odd p it bounds.
            reduced = (2 % prime,)
            if prime == 2:
                self.assertEqual(reduced, (0,))
            else:
                self.assertEqual(field_cone.boundary(2, (2 * pow(4, -1, prime) % prime,)), reduced)
        bockstein = LiftedBockstein(cone, 2)
        self.assertEqual(bockstein.differential(1)[2], ((0,),))
        self.assertEqual(bockstein.differential(2)[2], ((1,),))
        self.assertEqual(bockstein.page(3), {})

    def test_tampering_and_input_binding(self):
        c = IntegralChainComplex([0, 1, 1], {2: [[4]]})
        receipt = c.class_certificate(1, [2])
        forged = copy.deepcopy(receipt)
        forged.update(order=4, lift=[2])
        self.assertFalse(c.verify_class_certificate(1, [2], forged))
        forged = copy.deepcopy(receipt)
        forged["witnesses"][0]["modulus"] = 2
        self.assertFalse(c.verify_class_certificate(1, [2], forged))
        self.assertFalse(c.verify_class_certificate(1, [6], receipt))
        other = IntegralChainComplex([0, 1, 1], {2: [[8]]})
        self.assertFalse(other.verify_class_certificate(1, [2], receipt))

    def test_nondiagonal_complex_all_outcomes(self):
        c = IntegralChainComplex([1, 3, 2], {1: [[2, 4, 0]], 2: [[-8, 0], [4, 0], [0, 0]]})
        for z, status in [([1, 0, 0], "noncycle"), ([-8, 4, 0], "boundary"),
                          ([-4, 2, 0], "torsion"), ([-4, 2, 1], "infinite")]:
            with self.subTest(z=z):
                receipt = c.class_certificate(1, z)
                self.assertEqual(receipt["status"], status)
                self.assertTrue(c.verify_class_certificate(1, z, receipt))

    def test_multiple_moduli_and_empty_groups(self):
        c = IntegralChainComplex([0, 2, 2], {2: [[4, 0], [0, 9]]})
        receipt = c.class_certificate(1, [2, 1])
        self.assertEqual(receipt["order"], 18)
        self.assertEqual(len(receipt["witnesses"]), 2)
        self.assertTrue(c.verify_class_certificate(1, [2, 1], receipt))
        for c, degree, z in [(IntegralChainComplex([0, 2], {}), 0, []),
                             (IntegralChainComplex([1], {}), 0, [1])]:
            self.assertTrue(c.verify_class_certificate(degree, z, c.class_certificate(degree, z)))

    def test_mapping_cone_sign_and_chain_map_validation(self):
        c = IntegralChainComplex([1, 1], {1: [[3]]})
        cone = integral_mapping_cone(c, c, {0: [[2]], 1: [[2]]})
        self.assertEqual(cone.boundaries, {1: [[3, 2]], 2: [[2], [-3]]})
        with self.assertRaisesRegex(ValueError, "!= 0"):
            integral_mapping_cone(c, c, {0: [[2]], 1: [[1]]})


if __name__ == "__main__":
    unittest.main()
