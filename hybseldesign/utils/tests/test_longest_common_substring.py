"""Tests for longest_common_substring module.
"""

# Author: Hayden Metsky <hayden@mit.edu>

import unittest

from hybseldesign.utils import longest_common_substring as lcf


"""Tests the longest common substring with k mismatches
implementation.
"""
class TestLongestCommonSubstringWithKMismatches(unittest.TestCase):

  def test_different(self):
    a = 'ABC'
    b = 'DEF'
    self.assertEqual(lcf.k_lcf(a,b,0), (0,0,0))
    self.assertEqual(lcf.k_lcf(a,b,1), (1,0,2))
    self.assertEqual(lcf.k_lcf(a,b,2), (2,0,1))
    self.assertEqual(lcf.k_lcf(a,b,3), (3,0,0))
    self.assertEqual(lcf.k_lcf(a,b,4), (3,0,0))

  def test_equal(self):
    a = 'ABCDEFGHIJ'
    b = 'ABCDEFGHIJ'
    for k in xrange(0, 15):
      self.assertEqual(lcf.k_lcf(a,b,k), (10,0,0))

  def test_almost_equal1(self):
    a = 'ABCDEFGHIJ'
    b = 'ABCDEFGKIJ'
    self.assertEqual(lcf.k_lcf(a,b,0), (7,0,0))
    for k in xrange(1, 15):
      self.assertEqual(lcf.k_lcf(a,b,k), (10,0,0))

  def test_almost_equal2(self):
    a = 'ABCABCDEFGHIJKLM'
    b = 'DEFGABCDEFGHIJKLN'
    self.assertEqual(lcf.k_lcf(a,b,0), (12,3,4))
    self.assertEqual(lcf.k_lcf(a,b,1), (13,2,3))
    self.assertEqual(lcf.k_lcf(a,b,2), (14,1,2))

  def test_complex1(self):
    a = 'ABCDEFLONGESTCOMMONSUBSTRINGQRSTU'
    b = 'QRSLONGESTCOMMONTUBVTRINGVYX'
    self.assertEqual(lcf.k_lcf(a,b,0), (13,6,3))
    self.assertEqual(lcf.k_lcf(a,b,1), (16,6,3))
    self.assertEqual(lcf.k_lcf(a,b,2), (22,6,3))
    self.assertEqual(lcf.k_lcf(a,b,3), (23,5,2))
    self.assertEqual(lcf.k_lcf(a,b,4), (24,4,1))
    self.assertEqual(lcf.k_lcf(a,b,5), (25,3,0))
    self.assertEqual(lcf.k_lcf(a,b,6), (26,3,0))
    self.assertEqual(lcf.k_lcf(a,b,7), (27,3,0))

  def test_complex2(self):
    a = 'AGTCGCTGCCTCGTGCACATTG'
    b = 'GTATAATGTCGCAGCGTCGGCC'
    self.assertEqual(lcf.k_lcf(a,b,0), (5,1,7))
    self.assertEqual(lcf.k_lcf(a,b,1), (8,1,7))
    self.assertEqual(lcf.k_lcf(a,b,2), (12,1,7))
    self.assertEqual(lcf.k_lcf(a,b,3), (13,0,6))
    self.assertEqual(lcf.k_lcf(a,b,4), (15,1,7))
    # flip above
    a, b = b, a
    self.assertEqual(lcf.k_lcf(a,b,0), (5,7,1))
    self.assertEqual(lcf.k_lcf(a,b,1), (8,7,1))
    self.assertEqual(lcf.k_lcf(a,b,2), (12,7,1))
    self.assertEqual(lcf.k_lcf(a,b,3), (13,6,0))
    self.assertEqual(lcf.k_lcf(a,b,4), (15,7,1))

