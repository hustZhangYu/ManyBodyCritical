import math
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from basis import nlabel
import JumpOperator
import MbcH


def dense_reference(L, n):
    number = np.zeros((2**L, 2**L))
    for site in range(L):
        number += np.kron(np.eye(2**site),
                          np.kron([[1, 0], [0, 0]], np.eye(2**(L-site-1))))
    return [index for index in range(2**L) if number[index, index] == n]


class BasisTests(unittest.TestCase):
    def test_every_small_sector_matches_kronecker_order(self):
        for L in range(7):
            for n in range(-1, L + 2):
                with self.subTest(L=L, n=n):
                    self.assertEqual(nlabel(L, n), dense_reference(L, n))

    def test_larger_sector_without_dense_allocation(self):
        labels = nlabel(14, 7)
        self.assertEqual(len(labels), math.comb(14, 7))
        self.assertEqual(labels, sorted(labels))

    def test_both_existing_modules_reexport_the_same_function(self):
        self.assertIs(JumpOperator.nlabel, nlabel)
        self.assertIs(MbcH.nlabel, nlabel)


if __name__ == "__main__":
    unittest.main()
