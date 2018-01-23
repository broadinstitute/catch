"""Tests for longest_common_substring module.
"""

import unittest

from catch.utils import longest_common_substring as lcf

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class TestLCSWithKMismatches(unittest.TestCase):
    """Tests the k_lcf function.
    """

    def test_different(self):
        a = 'ABC'
        b = 'DEF'
        self.assertEqual(lcf.k_lcf(a, b, 0), (0, 0, 0))
        self.assertEqual(lcf.k_lcf(a, b, 1), (1, 0, 2))
        self.assertEqual(lcf.k_lcf(a, b, 2), (2, 0, 1))
        self.assertEqual(lcf.k_lcf(a, b, 3), (3, 0, 0))
        self.assertEqual(lcf.k_lcf(a, b, 4), (3, 0, 0))

    def test_equal(self):
        a = 'ABCDEFGHIJ'
        b = 'ABCDEFGHIJ'
        for k in range(0, 15):
            self.assertEqual(lcf.k_lcf(a, b, k), (10, 0, 0))

    def test_almost_equal1(self):
        a = 'ABCDEFGHIJ'
        b = 'ABCDEFGKIJ'
        self.assertEqual(lcf.k_lcf(a, b, 0), (7, 0, 0))
        for k in range(1, 15):
            self.assertEqual(lcf.k_lcf(a, b, k), (10, 0, 0))

    def test_almost_equal2(self):
        a = 'ABCABCDEFGHIJKLM'
        b = 'DEFGABCDEFGHIJKLN'
        self.assertEqual(lcf.k_lcf(a, b, 0), (12, 3, 4))
        self.assertEqual(lcf.k_lcf(a, b, 1), (13, 2, 3))
        self.assertEqual(lcf.k_lcf(a, b, 2), (14, 1, 2))

    def test_complex1(self):
        a = 'ABCDEFLONGESTCOMMONSUBSTRINGQRSTU'
        b = 'QRSLONGESTCOMMONTUBVTRINGVYX'
        self.assertEqual(lcf.k_lcf(a, b, 0), (13, 6, 3))
        self.assertEqual(lcf.k_lcf(a, b, 1), (16, 6, 3))
        self.assertEqual(lcf.k_lcf(a, b, 2), (22, 6, 3))
        self.assertEqual(lcf.k_lcf(a, b, 3), (23, 5, 2))
        self.assertEqual(lcf.k_lcf(a, b, 4), (24, 4, 1))
        self.assertEqual(lcf.k_lcf(a, b, 5), (25, 3, 0))
        self.assertEqual(lcf.k_lcf(a, b, 6), (26, 3, 0))
        self.assertEqual(lcf.k_lcf(a, b, 7), (27, 3, 0))

    def test_complex2(self):
        a = 'AGTCGCTGCCTCGTGCACATTG'
        b = 'GTATAATGTCGCAGCGTCGGCC'
        self.assertEqual(lcf.k_lcf(a, b, 0), (5, 1, 7))
        self.assertEqual(lcf.k_lcf(a, b, 1), (8, 1, 7))
        self.assertEqual(lcf.k_lcf(a, b, 2), (12, 1, 7))
        self.assertEqual(lcf.k_lcf(a, b, 3), (13, 0, 6))
        self.assertEqual(lcf.k_lcf(a, b, 4), (15, 1, 7))
        # flip above
        a, b = b, a
        self.assertEqual(lcf.k_lcf(a, b, 0), (5, 7, 1))
        self.assertEqual(lcf.k_lcf(a, b, 1), (8, 7, 1))
        self.assertEqual(lcf.k_lcf(a, b, 2), (12, 7, 1))
        self.assertEqual(lcf.k_lcf(a, b, 3), (13, 6, 0))
        self.assertEqual(lcf.k_lcf(a, b, 4), (15, 7, 1))


class TestLCSAroundAnchorWithKMismatches(unittest.TestCase):
    """Tests the k_lcf_around_anchor function.
    """

    def test_just_anchor(self):
        a = 'ABCDEFGHIJKLM'
        b = 'XYZDEFSTUVWXY'
        self.assertEqual(lcf.k_lcf_around_anchor(a, b, 3, 6, 0), (3, 3))

    def test_anchor_and_around(self):
        a = 'ABCDEFGHIJKLM'
        b = 'XYCDEFSTUVWXY'
        self.assertEqual(lcf.k_lcf_around_anchor(a, b, 3, 6, 0), (4, 2))
        a = 'ABCDEFGHIJKLM'
        b = 'XYCDEFGTUVWXY'
        self.assertEqual(lcf.k_lcf_around_anchor(a, b, 3, 6, 0), (5, 2))
        a = 'ABCDEFGHIJKLM'
        b = 'XYCDEFGHUVWXY'
        self.assertEqual(lcf.k_lcf_around_anchor(a, b, 3, 6, 0), (6, 2))
        a = 'ABCDEFGHIJKLM'
        b = 'XBCDEFSTUVWXY'
        self.assertEqual(lcf.k_lcf_around_anchor(a, b, 3, 6, 0), (5, 1))
        a = 'ABCDEFGHIJKLM'
        b = 'ABCDEFSTUVWXY'
        self.assertEqual(lcf.k_lcf_around_anchor(a, b, 3, 6, 0), (6, 0))
        a = 'ABCDEFGHIJKLM'
        b = 'XYCDEFGHIJKLM'
        self.assertEqual(lcf.k_lcf_around_anchor(a, b, 3, 6, 0), (11, 2))

    def test_anchor_start(self):
        a = 'ABCDEFGHIJKLM'
        b = 'ABCDEFSTUVWXY'
        self.assertEqual(lcf.k_lcf_around_anchor(a, b, 0, 3, 0), (6, 0))

    def test_anchor_end(self):
        a = 'ABCDEFGHIJKLM'
        b = 'XYCDEFSTUJKLM'
        self.assertEqual(lcf.k_lcf_around_anchor(a, b, 10, 13, 0), (4, 9))

    def test_more_mismatches_than_needed(self):
        a = 'ABCDEFGHIJKLM'
        b = 'ABCDEFGHIJKLM'
        self.assertEqual(lcf.k_lcf_around_anchor(a, b, 3, 6, 0), (13, 0))
        self.assertEqual(lcf.k_lcf_around_anchor(a, b, 3, 6, 1), (13, 0))
        self.assertEqual(lcf.k_lcf_around_anchor(a, b, 3, 6, 2), (13, 0))
        a = 'ABCDEFGHIJKLM'
        b = 'AXCDEFGHIJKLM'
        self.assertEqual(lcf.k_lcf_around_anchor(a, b, 3, 6, 0), (11, 2))
        self.assertEqual(lcf.k_lcf_around_anchor(a, b, 3, 6, 1), (13, 0))
        self.assertEqual(lcf.k_lcf_around_anchor(a, b, 3, 6, 2), (13, 0))

    def test_with_one_mismatch(self):
        a = 'ABCDEFGHIJKLM'
        b = 'XBZDEFSTUVWXY'
        self.assertEqual(lcf.k_lcf_around_anchor(a, b, 3, 6, 0), (3, 3))
        a = 'ABCDEFGHIJKLM'
        b = 'XBZDEFSTUVWXY'
        self.assertEqual(lcf.k_lcf_around_anchor(a, b, 3, 6, 1), (5, 1))
        a = 'ABCDEFGHIJKLM'
        b = 'XBZDEFSHIVWXY'
        self.assertEqual(lcf.k_lcf_around_anchor(a, b, 3, 6, 1), (6, 3))
        a = 'ABCDEFGHIJKLM'
        b = 'XYZDEFSTUVWXY'
        self.assertEqual(lcf.k_lcf_around_anchor(a, b, 3, 6, 1), (4, 3))
        a = 'ABCDEFGHIJKLM'
        b = 'XYZDEFGTUVWXY'
        self.assertEqual(lcf.k_lcf_around_anchor(a, b, 3, 6, 1), (5, 3))

    def test_with_two_mismatches(self):
        a = 'ABCDEFGHIJKLM'
        b = 'XBZDEFSTUVWXY'
        self.assertEqual(lcf.k_lcf_around_anchor(a, b, 3, 6, 2), (6, 1))
        a = 'ABCDEFGHIJKLM'
        b = 'ABZDEFSTUVWXY'
        self.assertEqual(lcf.k_lcf_around_anchor(a, b, 3, 6, 2), (7, 0))
        a = 'ABCDEFGHIJKLM'
        b = 'XBZDEFGTUJKLM'
        self.assertEqual(lcf.k_lcf_around_anchor(a, b, 3, 6, 2), (10, 3))
        a = 'ABCDEFGHIJKLM'
        b = 'XBZDEFSHUVWXY'
        self.assertEqual(lcf.k_lcf_around_anchor(a, b, 3, 6, 2), (7, 1))

    def test_with_three_mismatches(self):
        a = 'ABCDEFGHIJKLM'
        b = 'ABZDEFSTIJWXY'
        self.assertEqual(lcf.k_lcf_around_anchor(a, b, 3, 6, 3), (10, 0))
