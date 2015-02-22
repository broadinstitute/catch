"""Tests for probe module.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

import unittest
import numpy as np

from hybseldesign import probe


"""Tests methods in the Probe class.
"""
class TestProbe(unittest.TestCase):

  def setUp(self):
    self.a = probe.Probe.from_str('ATCGTCGCGGATCG')
    self.b = probe.Probe.from_str('ATCCTCGCGTATNG')
    self.c = probe.Probe.from_str('ATCGTCGCGGATC')
    self.d = probe.Probe.from_str('GATCGTCGCGGATC')
    self.e = probe.Probe.from_str('GGATTGTCGGGGAT')
    self.f = probe.Probe.from_str('GTCGCGGAACGGGG')
    self.g = probe.Probe.from_str('GTCGCTGATCGATC')

  """Test that probe parses the string correctly.
  """
  def test_parse_str(self):
    np.testing.assert_array_equal(self.a.seq,
        np.array(['A','T','C','G','T','C','G','C','G',
                  'G','A','T','C','G']))

  """Test mismatches method.
  """
  def test_mismatches(self):
    self.assertEqual(self.a.mismatches(self.a), 0)
    self.assertEqual(self.a.mismatches(self.b), 3)
    self.assertEqual(self.b.mismatches(self.a), 3)

  """Test mismatches_at_offset method.
  """
  def test_mismatches_at_offset(self):
    self.assertEqual(self.a.mismatches_at_offset(self.d, -1), 0)
    self.assertEqual(self.a.mismatches_at_offset(self.e, -2), 2)
    self.assertEqual(self.a.mismatches_at_offset(self.f, 3), 1)
    self.assertRaises(ValueError, self.a.mismatches_at_offset,
        self.c, 1)
    self.assertRaises(ValueError, self.a.mismatches_at_offset,
        self.b, 15)

  """Test min_mismatches_within_shift method.
  """
  def test_min_mismatches_within_shift(self):
    self.assertEqual(self.a.min_mismatches_within_shift(self.g, 5), 1)
    self.assertEqual(self.g.min_mismatches_within_shift(self.a, 5), 1)
    self.assertEqual(self.a.min_mismatches_within_shift(self.g, 2), 8)
    self.assertEqual(self.g.min_mismatches_within_shift(self.a, 2), 8)
    self.assertEqual(self.a.min_mismatches_within_shift(self.b, 0), 3)
    self.assertEqual(self.b.min_mismatches_within_shift(self.a, 0), 3)
    self.assertEqual(self.a.min_mismatches_within_shift(self.b, 2), 3)
    self.assertEqual(self.b.min_mismatches_within_shift(self.a, 2), 3)

  """Test share_some_kmers method.
  """
  def test_share_some_kmers_nonmemoized(self):
    np.random.seed(1)
    args = { 'k': 5, 'num_kmers_to_test': 10,
             'memoize_kmer_hashes': False }
    a = probe.Probe.from_str('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
    b = probe.Probe.from_str('ZYXWVUTSRQPONMLKJIHGFEDCBA')
    c = probe.Probe.from_str('ABCXDEFGHIJKLMNOPQRATUVWYZ')
    ab, ba, ac, ca = 0, 0, 0, 0
    for i in xrange(100):
      if a.shares_some_kmers(b, **args):
        ab += 1
      if b.shares_some_kmers(a, **args):
        ba += 1
      if a.shares_some_kmers(c, **args):
        ac += 1
      if c.shares_some_kmers(a, **args):
        ca += 1
    self.assertLess(ab, 10)
    self.assertLess(ba, 10)
    self.assertGreater(ac, 90)
    self.assertGreater(ca, 90)

  """Test share_some_kmers method.
  """
  def test_share_some_kmers_memoized(self):
    np.random.seed(1)
    args = { 'k': 5, 'num_kmers_to_test': 10,
             'memoize_kmer_hashes': True }
    a = probe.Probe.from_str('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
    b = probe.Probe.from_str('ZYXWVUTSRQPONMLKJIHGFEDCBA')
    c = probe.Probe.from_str('ABCXDEFGHIJKLMNOPQRATUVWYZ')
    ab, ba, ac, ca = 0, 0, 0, 0
    for i in xrange(100):
      if a.shares_some_kmers(b, **args):
        ab += 1
      if b.shares_some_kmers(a, **args):
        ba += 1
      if a.shares_some_kmers(c, **args):
        ac += 1
      if c.shares_some_kmers(a, **args):
        ca += 1
    self.assertLess(ab, 10)
    self.assertLess(ba, 10)
    self.assertGreater(ac, 90)
    self.assertGreater(ca, 90)

