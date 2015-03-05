"""Tests for probe module.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

import unittest
import numpy as np
from collections import defaultdict

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
        if kmer in kmer_map and p in kmer_map[kmer]:
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

  def test_positions(self):
    np.random.seed(1)
    a = probe.Probe.from_str('ABCDEFGABC')
    b = probe.Probe.from_str('XYZDEFHGHI')
    probes = [a, b]
    # Use a high num_kmers_per_probe to ensure all possible
    # kmers are selected to be put into the map
    kmer_map = probe.construct_kmer_probe_map(probes, k=3,
        num_kmers_per_probe=50, include_positions=True)
    self.assertItemsEqual(kmer_map['DEF'], [(a,3), (b,3)])
    self.assertItemsEqual(kmer_map['ABC'], [(a,0), (a,7)])
    self.assertItemsEqual(kmer_map['XYZ'], [(b,0)])
    self.assertItemsEqual(kmer_map['EFG'], [(a,4)])
    self.assertItemsEqual(kmer_map['EFH'], [(b,4)])


"""Tests probe_covers_sequence_by_longest_common_substring
function.
"""
class TestProbeCoversSequenceByLongestCommonSubstring(unittest.TestCase):

  def setUp(self):
    self.seq = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  def test_match_no_mismatches(self):
    f = probe.probe_covers_sequence_by_longest_common_substring(0, 6)
    match = f('ZZZABCGHIJKLXYZ', self.seq, 6, 9)
    self.assertTrue(match != None)
    start, end = match
    self.assertEqual(start, 6)
    self.assertEqual(end, 12)
    match = f('ZZZZAFGHIJKLMDEF', self.seq, 6, 9)
    self.assertTrue(match != None)
    start, end = match
    self.assertEqual(start, 5)
    self.assertEqual(end, 13)

  def test_match_with_mismatches(self):
    f = probe.probe_covers_sequence_by_longest_common_substring(1, 6)
    match = f('ZZZGHIGHIXKLDEF', self.seq, 6, 9)
    self.assertTrue(match != None)
    start, end = match
    self.assertEqual(start, 6)
    self.assertEqual(end, 12)
    match = f('ZZZZZZGHIJKXSWZ', self.seq, 6, 9)
    self.assertTrue(match != None)
    start, end = match
    self.assertEqual(start, 6)
    self.assertEqual(end, 12)
    match = f('ZZAGTFGHIJKXM', self.seq, 6, 9)
    self.assertTrue(match != None)
    start, end = match
    self.assertEqual(start, 5)
    self.assertEqual(end, 13)

  def test_no_match_no_mismatches(self):
    f = probe.probe_covers_sequence_by_longest_common_substring(0, 6)
    match = f('ZZZABCGHIXKLXYZ', self.seq, 6, 9)
    self.assertTrue(match == None)
    match = f('ZZZZZAGHIJKBC', self.seq, 6, 9)
    self.assertTrue(match == None)
    
  def test_no_match_with_mismatches(self):
    f = probe.probe_covers_sequence_by_longest_common_substring(1, 6)
    match = f('ZZZABCGHIXKXXYZ', self.seq, 6, 9)
    self.assertTrue(match == None)
    f = probe.probe_covers_sequence_by_longest_common_substring(1, 6)
    match = f('ZZZZZAGHIJYZ', self.seq, 6, 9)
    self.assertTrue(match == None)


"""Tests find_probe_covers_in_sequence function.
"""
class TestFindProbeCoversInSequence(unittest.TestCase):

  """Tests with short sequence, short probes, and small k
  where each probe appears zero or one times.
  """
  def test_one_or_no_occurrence(self):
    np.random.seed(1)
    sequence = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    a = probe.Probe.from_str('GHIJKL')
    b = probe.Probe.from_str('STUVWX')
    c = probe.Probe.from_str('ACEFHJ')
    probes = [a, b, c]
    k = 2
    num_kmers_per_probe = 10
    kmer_probe_map = probe.construct_kmer_probe_map(probes, k,
        num_kmers_per_probe, include_positions=True)
    f = probe.probe_covers_sequence_by_longest_common_substring(0, 6)
    found = probe.find_probe_covers_in_sequence(sequence,
              kmer_probe_map, k, f)
    self.assertItemsEqual(found[a], [(6,12)])
    self.assertItemsEqual(found[b], [(18,24)])
    self.assertFalse(c in found)

  """Tests with short sequence, short probes, and small k
  where one probe appears twice.
  """
  def test_two_occurrences(self):
    np.random.seed(1)
    sequence = 'ABCDEFGHIJKLMNOPCDEFGHQRSTU'
    a = probe.Probe.from_str('CDEFGH')
    b = probe.Probe.from_str('GHIJKL')
    c = probe.Probe.from_str('STUVWX')
    probes = [a, b, c]
    k = 2
    num_kmers_per_probe = 10
    kmer_probe_map = probe.construct_kmer_probe_map(probes, k,
        num_kmers_per_probe, include_positions=True)
    f = probe.probe_covers_sequence_by_longest_common_substring(0, 6)
    found = probe.find_probe_covers_in_sequence(sequence,
              kmer_probe_map, k, f)
    self.assertItemsEqual(found[a], [(2,8), (16,22)])
    self.assertItemsEqual(found[b], [(6,12)])
    self.assertFalse(c in found)

  """Tests with short sequence, short probes, and small k
  where probes contain more than what they cover.
  """
  def test_more_than_cover(self):
    np.random.seed(1)
    sequence = 'ABCDEFGHIJKLMNOPQR'+('Z'*100)+'STUVWXYZ'
    a = probe.Probe.from_str('XYZCDEFGHIJKABCSTUVWXABC')
    b = probe.Probe.from_str('PQRSGHIJKLMNXYZ')
    c = probe.Probe.from_str('ABCFGHIJKLZAZAZAGHIJKL')
    probes = [a, b, c]
    k = 4
    num_kmers_per_probe = 100
    kmer_probe_map = probe.construct_kmer_probe_map(probes, k,
        num_kmers_per_probe, include_positions=True)
    f = probe.probe_covers_sequence_by_longest_common_substring(0, 6)
    found = probe.find_probe_covers_in_sequence(sequence,
              kmer_probe_map, k, f)
    self.assertItemsEqual(found[a], [(2,11), (118,124)])
    self.assertItemsEqual(found[b], [(6,14)])
    self.assertItemsEqual(found[c], [(5,12)])

  """Tests with short sequence, short probes, and small k
  where the sequence and probes have repetitive sequences, so that
  one probe can cover a lot of the sequence.
  """
  def test_repetitive(self):
    np.random.seed(1)
    sequence = 'ABCAAAAAAAAAAXYZXYZXYZXYZAAAAAAAAAAAAAXYZ'
    a = probe.Probe.from_str('NAAAAAAN')
    probes = [a]
    k = 3
    num_kmers_per_probe = 20
    kmer_probe_map = probe.construct_kmer_probe_map(probes, k,
        num_kmers_per_probe, include_positions=True)
    f = probe.probe_covers_sequence_by_longest_common_substring(0, 6)
    found = probe.find_probe_covers_in_sequence(sequence,
              kmer_probe_map, k, f)
    self.assertItemsEqual(found[a], [(3,13), (25,38)])

  def test_random_small_genome(self):
    self.run_random(100, 15000, 25000, 300)

  def test_random_large_genome(self):
    self.run_random(2, 1500000, 2500000, 30000)

  """Runs, n times, a test in which a long sequence is randomly
  generated and probes are generated from that sequence. Makes 100 bp
  probes with k=15 and num_kmers_per_probe=10 and creates the probes
  with the intention of determining coverage with a longest common
  substring. The genome size is randomly chosen between genome_min
  and genome_max, and the number of probes generated is num_probes.
  """
  def run_random(self, n, genome_min, genome_max, num_probes):
    np.random.seed(1)
    for n in xrange(n):
      # Make a random sequence
      seq_length = np.random.randint(genome_min, genome_max)
      sequence = "".join(np.random.choice(['A','T','C','G'],
          size=seq_length, replace=True))
      desired_probe_cover_ranges = defaultdict(list)
      # Make num_probes random probes
      probes = []
      for m in xrange(num_probes):
        probe_length = 100
        subseq_start = np.random.randint(0, seq_length-probe_length)
        subseq_end = subseq_start + probe_length
        cover_length = np.random.randint(80, 100)
        cover_start = subseq_start + \
            np.random.randint(0, probe_length-cover_length+1)
        cover_end = min(seq_length, cover_start + cover_length)
        probe_str_cover = sequence[cover_start:cover_end]
        # Add random bases before and after what the probe should
        # cover
        probe_str_start = \
            "".join(np.random.choice(['A','T','C','G'],
            size=cover_start-subseq_start, replace=True))
        probe_str_end = \
            "".join(np.random.choice(['A','T','C','G'],
            size=subseq_end-cover_end, replace=True))
        probe_str = probe_str_start + probe_str_cover + probe_str_end
        # Add 0, 1, 2, or 3 random mismatches
        for k in xrange(np.random.randint(0, 4)):
          pos = np.random.randint(0, probe_length)
          base_choices = [b for b in ['A','T','C','G'] if b != probe_str[pos]]
          probe_str = probe_str[:pos] + \
              "".join(np.random.choice(base_choices, size=1)) + \
              probe_str[(pos+1):]
        p = probe.Probe.from_str(probe_str)
        desired_probe_cover_ranges[p].append((cover_start, cover_end))
        probes += [p]
      kmer_probe_map = probe.construct_kmer_probe_map(probes,
          k=15, num_kmers_per_probe=10, include_positions=True)
      f = probe.probe_covers_sequence_by_longest_common_substring(3, 80)
      found = probe.find_probe_covers_in_sequence(sequence,
          kmer_probe_map, k=15,
          cover_range_for_probe_in_subsequence_fn=f)
      # Check that this didn't find any extraneous probes and that
      # it found at least 95% of the original (it may miss some
      # due to false negatives in the approach)
      self.assertLessEqual(len(found), len(probes))
      self.assertGreaterEqual(len(found), 0.95*len(probes))
      # Check that each desired probe was found correctly
      for p, cover_ranges in desired_probe_cover_ranges.iteritems():
        if p not in found:
          continue
        found_cover_ranges = found[p]
        # This probe most likely was found once, but could have
        # been missed (due to false negatives in the approach) and
        # may have been found more than once due to chance (but
        # probably not too much more!)
        self.assertTrue(len(found_cover_ranges) in [1,2])
        # The cover ranges should have been captured, and the ones
        # found may extend past what was desired by a small amount due
        # to allowing mismatches and chance
        # Because of mismatches possibly added to the end of the
        # desired cover range, what was recaptured may not always
        # encompass the entire cover range, so allow some small tolerance
        for desired_cv in cover_ranges:
          found_desired_cv = False
          for found_cv in found_cover_ranges:
            left_diff = desired_cv[0] - found_cv[0]
            right_diff = found_cv[1] - desired_cv[1]
            if left_diff >= -5 and left_diff < 15:
              if right_diff >= -5 and right_diff < 15:
                found_desired_cv = True
                break
          self.assertTrue(found_desired_cv)

