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
             'memoize_kmers': False }
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
             'memoize_kmers': True }
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

  """Test construct_kmers method.
  """
  def test_construct_kmers(self):
    a = probe.Probe.from_str('ABCDEFGHI')
    self.assertEqual(a.construct_kmers(4),
        set(['ABCD','BCDE','CDEF','DEFG','EFGH','FGHI']))


"""Tests construct_kmer_probe_map function.
"""
class TestConstructKmerProbeMap(unittest.TestCase):

  def make_random_probe(self, length):
    bases = ['A','T','C','G']
    s = "".join(np.random.choice(bases, size=length, replace=True))
    return probe.Probe.from_str(s)

  """Make 50 random probes. From them, construct the kmer probe map
  with k=15 and with 10 kmers per probe. For each probe, check that
  at least 8 of its kmers can be found in this map (because kmers
  are selected from the probes with replacement, not all 10 may be
  present, and indeed not even 8 may be present).
  """
  def test_random(self):
    np.random.seed(1)
    k = 15
    num_kmers_per_probe = 10
    probes = [self.make_random_probe(100) for _ in xrange(50)]
    kmer_map = probe.construct_kmer_probe_map(probes, k=k,
        num_kmers_per_probe=num_kmers_per_probe)
    for p in probes:
      num_found = 0
      for kmer in p.construct_kmers(k):
        if p in kmer_map[kmer]:
          num_found += 1
      self.assertGreaterEqual(num_found, 0.8*num_kmers_per_probe)

  def test_shared_kmer(self):
    np.random.seed(1)
    a = probe.Probe.from_str('ABCDEFG')
    b = probe.Probe.from_str('XYZDEFH')
    probes = [a, b]
    # Use a high num_kmers_per_probe to ensure all possible
    # kmers are selected to be put into the map
    kmer_map = probe.construct_kmer_probe_map(probes, k=3,
        num_kmers_per_probe=50)
    self.assertTrue(a in kmer_map['DEF'])
    self.assertTrue(b in kmer_map['DEF'])
    self.assertTrue(a in kmer_map['ABC'])
    self.assertFalse(b in kmer_map['ABC'])
    self.assertFalse(a in kmer_map['XYZ'])
    self.assertTrue(b in kmer_map['XYZ'])
    self.assertTrue(a in kmer_map['EFG'])
    self.assertFalse(b in kmer_map['EFG'])
    self.assertFalse(a in kmer_map['EFH'])
    self.assertTrue(b in kmer_map['EFH'])


"""Tests the probe_covers_sequence_by_longest_common_substring
function.
"""
class TestProbeCoversSequenceByLongestCommonSubstring(unittest.TestCase):

  def setUp(self):
    self.seq = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  def test_match_no_mismatches(self):
    f = probe.probe_covers_sequence_by_longest_common_substring(0, 6)
    match = f(probe.Probe.from_str('GHIJKL'), self.seq)
    self.assertTrue(match != None)
    start, end = match
    self.assertEqual(start, 6)
    self.assertEqual(end, 12)
    match = f(probe.Probe.from_str('FGHIJKLM'), self.seq)
    self.assertTrue(match != None)
    start, end = match
    self.assertEqual(start, 5)
    self.assertEqual(end, 13)

  def test_match_with_mismatches(self):
    f = probe.probe_covers_sequence_by_longest_common_substring(1, 6)
    match = f(probe.Probe.from_str('GHIXKL'), self.seq)
    self.assertTrue(match != None)
    start, end = match
    self.assertEqual(start, 6)
    self.assertEqual(end, 12)
    match = f(probe.Probe.from_str('GHIJKX'), self.seq)
    self.assertTrue(match != None)
    start, end = match
    self.assertEqual(start, 6)
    self.assertEqual(end, 12)
    match = f(probe.Probe.from_str('FGHIJKXM'), self.seq)
    self.assertTrue(match != None)
    start, end = match
    self.assertEqual(start, 5)
    self.assertEqual(end, 13)

  def test_no_match_no_mismatches(self):
    f = probe.probe_covers_sequence_by_longest_common_substring(0, 6)
    match = f(probe.Probe.from_str('GHIXKL'), self.seq)
    self.assertTrue(match == None)
    match = f(probe.Probe.from_str('GHIJK'), self.seq)
    self.assertTrue(match == None)
    
  def test_no_match_with_mismatches(self):
    f = probe.probe_covers_sequence_by_longest_common_substring(1, 6)
    match = f(probe.Probe.from_str('GHIXKX'), self.seq)
    self.assertTrue(match == None)
    f = probe.probe_covers_sequence_by_longest_common_substring(1, 6)
    match = f(probe.Probe.from_str('GHIJ'), self.seq)
    self.assertTrue(match == None)


"""Tests the find_probe_covers_in_sequence function.
"""
class TestFindProbeCoversInSequence(unittest.TestCase):

  pass
