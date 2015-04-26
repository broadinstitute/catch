"""Tests for probe module.
"""

from collections import defaultdict
import unittest

import numpy as np

from hybseldesign import probe

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class TestProbe(unittest.TestCase):
    """Tests methods in the Probe class.
    """

    def setUp(self):
        self.a = probe.Probe.from_str('ATCGTCGCGGATCG')
        self.b = probe.Probe.from_str('ATCCTCGCGTATNG')
        self.c = probe.Probe.from_str('ATCGTCGCGGATC')
        self.d = probe.Probe.from_str('GATCGTCGCGGATC')
        self.e = probe.Probe.from_str('GGATTGTCGGGGAT')
        self.f = probe.Probe.from_str('GTCGCGGAACGGGG')
        self.g = probe.Probe.from_str('GTCGCTGATCGATC')

    def make_random_probe(self, length):
        bases = ['A', 'T', 'C', 'G']
        s = "".join(np.random.choice(bases, size=length, replace=True))
        return probe.Probe.from_str(s)

    def test_parse_str(self):
        """Test that probe parses the string correctly.
        """
        np.testing.assert_array_equal(self.a.seq,
                                      np.array(['A', 'T', 'C', 'G', 'T', 'C',
                                                'G', 'C', 'G', 'G', 'A', 'T',
                                                'C', 'G']))

    def test_mismatches(self):
        """Test mismatches method.
        """
        self.assertEqual(self.a.mismatches(self.a), 0)
        self.assertEqual(self.a.mismatches(self.b), 3)
        self.assertEqual(self.b.mismatches(self.a), 3)

    def test_mismatches_at_offset(self):
        """Test mismatches_at_offset method.
        """
        self.assertEqual(self.a.mismatches_at_offset(self.d, -1), 0)
        self.assertEqual(self.a.mismatches_at_offset(self.e, -2), 2)
        self.assertEqual(self.a.mismatches_at_offset(self.f, 3), 1)
        self.assertRaises(ValueError, self.a.mismatches_at_offset, self.c, 1)
        self.assertRaises(ValueError, self.a.mismatches_at_offset, self.b, 15)

    def test_min_mismatches_within_shift(self):
        """Test min_mismatches_within_shift method.
        """
        self.assertEqual(self.a.min_mismatches_within_shift(self.g, 5), 1)
        self.assertEqual(self.g.min_mismatches_within_shift(self.a, 5), 1)
        self.assertEqual(self.a.min_mismatches_within_shift(self.g, 2), 8)
        self.assertEqual(self.g.min_mismatches_within_shift(self.a, 2), 8)
        self.assertEqual(self.a.min_mismatches_within_shift(self.b, 0), 3)
        self.assertEqual(self.b.min_mismatches_within_shift(self.a, 0), 3)
        self.assertEqual(self.a.min_mismatches_within_shift(self.b, 2), 3)
        self.assertEqual(self.b.min_mismatches_within_shift(self.a, 2), 3)

    def test_reverse_complement(self):
        """Test reverse_complement method.
        """
        a_rc = self.a.reverse_complement()
        a_rc_desired = probe.Probe.from_str('CGATCCGCGACGAT')
        self.assertEqual(a_rc, a_rc_desired)

    def test_with_prepended_str(self):
        """Test with_prepended_str method.
        """
        a_prepended = self.a.with_prepended_str('TATA')
        a_prepended_desired = probe.Probe.from_str('TATAATCGTCGCGGATCG')
        self.assertEqual(a_prepended, a_prepended_desired)

    def test_with_appended_str(self):
        """Test with_appended_str method.
        """
        a_appended = self.a.with_appended_str('TATA')
        a_appended_desired = probe.Probe.from_str('ATCGTCGCGGATCGTATA')
        self.assertEqual(a_appended, a_appended_desired)

    def test_identifier(self):
        """Test identifier method.

        This randomly produces 100 probes and checks that their identifiers
        are all unique. They are not guaranteed to be, but certainly
        should be.
        """
        np.random.seed(1)
        probes = [self.make_random_probe(100) for _ in xrange(100)]
        identifiers = set([p.identifier() for p in probes])
        self.assertEqual(len(identifiers), 100)

    def test_share_some_kmers_nonmemoized(self):
        """Test share_some_kmers method.
        """
        np.random.seed(1)
        args = {'k': 5, 'num_kmers_to_test': 10, 'memoize_kmers': False}
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

    def test_share_some_kmers_memoized(self):
        """Test share_some_kmers method.
        """
        np.random.seed(1)
        args = {'k': 5, 'num_kmers_to_test': 10, 'memoize_kmers': True}
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

    def test_construct_kmers(self):
        """Test construct_kmers method.
        """
        a = probe.Probe.from_str('ABCDEFGHI')
        self.assertEqual(a.construct_kmers(4),
                         ['ABCD', 'BCDE', 'CDEF', 'DEFG', 'EFGH', 'FGHI'])


class TestConstructRandKmerProbeMap(unittest.TestCase):
    """Tests _construct_rand_kmer_probe_map function.
    """

    def make_random_probe(self, length):
        bases = ['A', 'T', 'C', 'G']
        s = "".join(np.random.choice(bases, size=length, replace=True))
        return probe.Probe.from_str(s)

    def test_random(self):
        """Make 50 random probes. From them, construct the kmer probe map
        with k=15 and with 10 kmers per probe. For each probe, check that
        at least 8 of its kmers can be found in this map (because kmers
        are selected from the probes with replacement, not all 10 may be
        present, and indeed not even 8 may be present).
        """
        np.random.seed(1)
        k = 15
        num_kmers_per_probe = 10
        probes = [self.make_random_probe(100) for _ in xrange(50)]
        kmer_map = probe._construct_rand_kmer_probe_map(
            probes,
            k=k,
            num_kmers_per_probe=num_kmers_per_probe)
        for p in probes:
            num_found = 0
            for kmer in p.construct_kmers(k):
                if kmer in kmer_map and p in kmer_map[kmer]:
                    num_found += 1
            self.assertGreaterEqual(num_found, 0.8 * num_kmers_per_probe)

    def test_shared_kmer(self):
        np.random.seed(1)
        a = probe.Probe.from_str('ABCDEFG')
        b = probe.Probe.from_str('XYZDEFH')
        probes = [a, b]
        # Use a high num_kmers_per_probe to ensure all possible
        # kmers are selected to be put into the map
        kmer_map = probe._construct_rand_kmer_probe_map(probes,
                                                        k=3,
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
        kmer_map = probe._construct_rand_kmer_probe_map(probes,
                                                        k=3,
                                                        num_kmers_per_probe=50,
                                                        include_positions=True)
        self.assertItemsEqual(kmer_map['DEF'], [(a, 3), (b, 3)])
        self.assertItemsEqual(kmer_map['ABC'], [(a, 0), (a, 7)])
        self.assertItemsEqual(kmer_map['XYZ'], [(b, 0)])
        self.assertItemsEqual(kmer_map['EFG'], [(a, 4)])
        self.assertItemsEqual(kmer_map['EFH'], [(b, 4)])


class TestConstructPigeonholedKmerProbeMap(unittest.TestCase):
    """Tests _construct_pigeonholed_kmer_probe_map function.
    """

    def test_no_mismatches(self):
        a = probe.Probe.from_str('ABCDEFGHIJ')
        b = probe.Probe.from_str('ZYXWVUTSRQ')
        probes = [a, b]
        kmer_map = probe._construct_pigeonholed_kmer_probe_map(
            probes, 0, min_k=5)
        # k-mers equal to the full length of the probe should be
        # chosen
        self.assertTrue(a in kmer_map[a.seq_str])
        self.assertTrue(b in kmer_map[b.seq_str])
        self.assertFalse(a in kmer_map[b.seq_str])
        self.assertFalse(b in kmer_map[a.seq_str])

    def test_too_small_k(self):
        a = probe.Probe.from_str('ABCDEFGHIJ')
        b = probe.Probe.from_str('ZYXWVUTSRQ')
        probes = [a, b]
        with self.assertRaises(probe.PigeonholeRequiresTooSmallKmerSizeError):
            # Should pick k=5, but requires k=6
            probe._construct_pigeonholed_kmer_probe_map(
                probes, 1, min_k=6)
        with self.assertRaises(probe.PigeonholeRequiresTooSmallKmerSizeError):
            # Should pick k=2, but requires k=3
            probe._construct_pigeonholed_kmer_probe_map(
                probes, 3, min_k=3)

    def test_one_mismatch(self):
        a = probe.Probe.from_str('ABCDEFGHIJ')
        b = probe.Probe.from_str('ZYXWVUTSRQ')
        probes = [a, b]
        kmer_map = probe._construct_pigeonholed_kmer_probe_map(
            probes, 1, min_k=2)
        # Should pick k=5
        self.assertEquals(len(kmer_map), 4)
        self.assertItemsEqual(kmer_map['ABCDE'], [a])
        self.assertItemsEqual(kmer_map['FGHIJ'], [a])
        self.assertItemsEqual(kmer_map['ZYXWV'], [b])
        self.assertItemsEqual(kmer_map['UTSRQ'], [b])

    def test_shared_kmer(self):
        a = probe.Probe.from_str('ABCDEFGHIJ')
        b = probe.Probe.from_str('ZYXWVABCDE')
        probes = [a, b]
        kmer_map = probe._construct_pigeonholed_kmer_probe_map(
            probes, 1, min_k=2)
        # Should pick k=5
        self.assertEquals(len(kmer_map), 3)
        self.assertItemsEqual(kmer_map['ABCDE'], [a, b])
        self.assertItemsEqual(kmer_map['FGHIJ'], [a])
        self.assertItemsEqual(kmer_map['ZYXWV'], [b])

    def test_positions(self):
        a = probe.Probe.from_str('ABCDEFGH')
        b = probe.Probe.from_str('ZYXWVUAB')
        probes = [a, b]
        kmer_map = probe._construct_pigeonholed_kmer_probe_map(
            probes, 3, min_k=2, include_positions=True)
        # Should pick k=2
        self.assertEquals(len(kmer_map), 7)
        self.assertItemsEqual(kmer_map['AB'], [(a, 0), (b, 6)])
        self.assertItemsEqual(kmer_map['CD'], [(a, 2)])
        self.assertItemsEqual(kmer_map['EF'], [(a, 4)])
        self.assertItemsEqual(kmer_map['GH'], [(a, 6)])
        self.assertItemsEqual(kmer_map['ZY'], [(b, 0)])
        self.assertItemsEqual(kmer_map['XW'], [(b, 2)])
        self.assertItemsEqual(kmer_map['VU'], [(b, 4)])


class TestProbeCoversSequenceByLongestCommonSubstring(unittest.TestCase):
    """Tests probe_covers_sequence_by_longest_common_substring function.
    """

    def setUp(self):
        self.seq = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    def test_match_no_mismatches(self):
        f = probe.probe_covers_sequence_by_longest_common_substring(0, 6)
        match = f('ZZZABCGHIJKLXYZ', self.seq, 6, 9)
        self.assertTrue(match is not None)
        start, end = match
        self.assertEqual(start, 6)
        self.assertEqual(end, 12)
        match = f('ZZZZAFGHIJKLMDEF', self.seq, 6, 9)
        self.assertTrue(match is not None)
        start, end = match
        self.assertEqual(start, 5)
        self.assertEqual(end, 13)

    def test_match_with_mismatches(self):
        f = probe.probe_covers_sequence_by_longest_common_substring(1, 6)
        match = f('ZZZGHIGHIXKLDEF', self.seq, 6, 9)
        self.assertTrue(match is not None)
        start, end = match
        self.assertEqual(start, 6)
        self.assertEqual(end, 12)
        match = f('ZZZZZZGHIJKXSWZ', self.seq, 6, 9)
        self.assertTrue(match is not None)
        start, end = match
        self.assertEqual(start, 6)
        self.assertEqual(end, 12)
        match = f('ZZAGTFGHIJKXM', self.seq, 6, 9)
        self.assertTrue(match is not None)
        start, end = match
        self.assertEqual(start, 5)
        self.assertEqual(end, 13)

    def test_no_match_no_mismatches(self):
        f = probe.probe_covers_sequence_by_longest_common_substring(0, 6)
        match = f('ZZZABCGHIXKLXYZ', self.seq, 6, 9)
        self.assertTrue(match is None)
        match = f('ZZZZZAGHIJKBC', self.seq, 6, 9)
        self.assertTrue(match is None)

    def test_no_match_with_mismatches(self):
        f = probe.probe_covers_sequence_by_longest_common_substring(1, 6)
        match = f('ZZZABCGHIXKXXYZ', self.seq, 6, 9)
        self.assertTrue(match is None)
        f = probe.probe_covers_sequence_by_longest_common_substring(1, 6)
        match = f('ZZZZZAGHIJYZ', self.seq, 6, 9)
        self.assertTrue(match is None)


class TestFindProbeCoversInSequence(unittest.TestCase):
    """Tests find_probe_covers_in_sequence function.
    """

    def test_one_or_no_occurrence(self):
        """Tests with short sequence and short probes
        where each probe appears zero or one times.
        """
        np.random.seed(1)
        sequence = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        a = probe.Probe.from_str('GHIJKL')
        b = probe.Probe.from_str('STUVWX')
        c = probe.Probe.from_str('ACEFHJ')
        probes = [a, b, c]
        kmer_map = probe.construct_kmer_probe_map_to_find_probe_covers(
            probes, 0, 6, min_k=6)
        f = probe.probe_covers_sequence_by_longest_common_substring(0, 6)
        found = probe.find_probe_covers_in_sequence(sequence, kmer_map, f)
        self.assertItemsEqual(found[a], [(6, 12)])
        self.assertItemsEqual(found[b], [(18, 24)])
        self.assertFalse(c in found)

    def test_two_occurrences(self):
        """Tests with short sequence and short probes
        where one probe appears twice.
        """
        np.random.seed(1)
        sequence = 'ABCDEFGHIJKLMNOPCDEFGHQRSTU'
        a = probe.Probe.from_str('CDEFGH')
        b = probe.Probe.from_str('GHIJKL')
        c = probe.Probe.from_str('STUVWX')
        probes = [a, b, c]
        kmer_map = probe.construct_kmer_probe_map_to_find_probe_covers(
            probes, 0, 6, min_k=6)
        f = probe.probe_covers_sequence_by_longest_common_substring(0, 6)
        found = probe.find_probe_covers_in_sequence(sequence, kmer_map, f)
        self.assertItemsEqual(found[a], [(2, 8), (16, 22)])
        self.assertItemsEqual(found[b], [(6, 12)])
        self.assertFalse(c in found)

    def test_more_than_cover(self):
        """Tests with short sequence and short probes
        where probes contain more than what they cover.
        """
        np.random.seed(1)
        sequence = 'ABCDEFGHIJKLMNOPQR' + ('Z' * 100) + 'STUVWXYZ'
        a = probe.Probe.from_str('XYZCDEFGHIJKABCSTUVWXABC')
        b = probe.Probe.from_str('PQRSGHIJKLMNXYZ')
        c = probe.Probe.from_str('ABCFGHIJKLZAZAZAGHIJKL')
        probes = [a, b, c]
        # This should default to the random approach, so set k (rather than
        # min_k)
        kmer_map = probe.construct_kmer_probe_map_to_find_probe_covers(
            probes, 0, 6, k=6)
        f = probe.probe_covers_sequence_by_longest_common_substring(0, 6)
        found = probe.find_probe_covers_in_sequence(sequence, kmer_map, f)
        self.assertItemsEqual(found[a], [(2, 11), (118, 124)])
        self.assertItemsEqual(found[b], [(6, 14)])
        self.assertItemsEqual(found[c], [(5, 12)])

    def test_repetitive(self):
        """Tests with short sequence and short probes
        where the sequence and probes have repetitive sequences, so that
        one probe can cover a lot of the sequence.
        """
        np.random.seed(1)
        sequence = 'ABCAAAAAAAAAAXYZXYZXYZXYZAAAAAAAAAAAAAXYZ'
        a = probe.Probe.from_str('NAAAAAAN')
        probes = [a]
        # This should default to the random approach, so set k (rather than
        # min_k)
        kmer_map = probe.construct_kmer_probe_map_to_find_probe_covers(
            probes, 0, 6, k=6)
        f = probe.probe_covers_sequence_by_longest_common_substring(0, 6)
        found = probe.find_probe_covers_in_sequence(sequence, kmer_map, f)
        self.assertItemsEqual(found[a], [(3, 13), (25, 38)])

    def test_pigeonhole_with_mismatch(self):
        """Tests with short sequence and short probes
        where the call to construct_kmer_probe_map_to_find_probe_covers tries
        the pigeonhole approach.
        """
        np.random.seed(1)
        sequence = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        a = probe.Probe.from_str('GHIJXL')
        b = probe.Probe.from_str('BTUVWX')
        c = probe.Probe.from_str('ACEFHJ')
        probes = [a, b, c]

        kmer_map = probe.construct_kmer_probe_map_to_find_probe_covers(
            probes, 1, 6, min_k=3, k=4)
        # This should try the pigeonhole approach, which should choose k=3
        for kmer in kmer_map.keys():
            self.assertEquals(len(kmer), 3)
        f = probe.probe_covers_sequence_by_longest_common_substring(1, 6)
        found = probe.find_probe_covers_in_sequence(sequence, kmer_map, f)
        self.assertItemsEqual(found[a], [(6, 12)])
        self.assertItemsEqual(found[b], [(18, 24)])
        self.assertFalse(c in found)

        kmer_map = probe.construct_kmer_probe_map_to_find_probe_covers(
            probes, 1, 6, min_k=4, k=4)
        # This should try the pigeonhole approach and fail because it
        # chooses k=3, but min_k=4. So it should then try the random
        # approach with k=4.
        for kmer in kmer_map.keys():
            self.assertEquals(len(kmer), 4)
        f = probe.probe_covers_sequence_by_longest_common_substring(1, 6)
        found = probe.find_probe_covers_in_sequence(sequence, kmer_map, f)
        self.assertItemsEqual(found[a], [(6, 12)])
        self.assertItemsEqual(found[b], [(18, 24)])
        self.assertFalse(c in found)

    def test_random_small_genome(self):
        self.run_random(100, 15000, 25000, 300)

    def test_random_large_genome(self):
        self.run_random(2, 1500000, 2500000, 30000)

    def run_random(self, n, genome_min, genome_max, num_probes):
        """Run tests with a randomly generated sequence.

        Repeatedly runs tests in which a sequence is randomly generated,
        probes are generated from that sequence, and then the probes are
        looked up in the sequence.

        Creates the probes with the intention of determining coverage with
        a longest common substring.

        Args:
            n: number of times to run the test
            genome_min/genome_max: the genome (sequence) size is
                randomly chosen between genome_min and genome_max
            num_probes: the number of probes generated from the random
                sequence
        """
        np.random.seed(1)
        for n in xrange(n):
            # Choose either lcf_thres=80 or lcf_thres=100
            lcf_thres = np.random.choice([80, 100])
            # Make a random sequence
            seq_length = np.random.randint(genome_min, genome_max)
            sequence = "".join(np.random.choice(['A', 'T', 'C', 'G'],
                                                size=seq_length,
                                                replace=True))
            desired_probe_cover_ranges = defaultdict(list)
            # Make num_probes random probes
            probes = []
            for m in xrange(num_probes):
                probe_length = 100
                subseq_start = np.random.randint(0, seq_length - probe_length)
                subseq_end = subseq_start + probe_length
                cover_length = np.random.randint(lcf_thres, 101)
                cover_start = subseq_start + \
                    np.random.randint(0, probe_length - cover_length + 1)
                cover_end = min(seq_length, cover_start + cover_length)
                probe_str_cover = sequence[cover_start:cover_end]
                # Add random bases before and after what the probe should
                # cover
                probe_str_start = "".join(
                    np.random.choice(['A', 'T', 'C', 'G'],
                                     size=cover_start - subseq_start,
                                     replace=True))
                probe_str_end = "".join(
                    np.random.choice(['A', 'T', 'C', 'G'],
                                     size=subseq_end - cover_end,
                                     replace=True))
                probe_str = probe_str_start + probe_str_cover + probe_str_end
                # Add 0, 1, 2, or 3 random mismatches
                for k in xrange(np.random.randint(0, 4)):
                    pos = np.random.randint(0, probe_length)
                    base_choices = [b for b in ['A', 'T', 'C', 'G']
                                    if b != probe_str[pos]]
                    probe_str = probe_str[:pos] + \
                        "".join(np.random.choice(base_choices, size=1)) + \
                        probe_str[(pos + 1):]
                p = probe.Probe.from_str(probe_str)
                desired_probe_cover_ranges[p].append((cover_start, cover_end))
                probes += [p]
            kmer_map = probe.construct_kmer_probe_map_to_find_probe_covers(
                probes, 3, lcf_thres)
            f = probe.probe_covers_sequence_by_longest_common_substring(
                3, lcf_thres)
            found = probe.find_probe_covers_in_sequence(
                sequence, kmer_map,
                cover_range_for_probe_in_subsequence_fn=f)
            # Check that this didn't find any extraneous probes and that
            # it found at least 95% of the original (it may miss some
            # due to false negatives in the approach)
            self.assertLessEqual(len(found), len(probes))
            self.assertGreaterEqual(len(found), 0.95 * len(probes))
            # Check that each desired probe was found correctly
            for p, cover_ranges in desired_probe_cover_ranges.iteritems():
                if p not in found:
                    continue
                found_cover_ranges = found[p]
                # This probe most likely was found once, but could have
                # been missed (due to false negatives in the approach) and
                # may have been found more than once due to chance (but
                # probably not too much more!)
                self.assertTrue(len(found_cover_ranges) in [1, 2])
                # The cover ranges should have been captured, and the ones
                # found may extend past what was desired by a small amount due
                # to allowing mismatches and chance
                # Because of mismatches possibly added to the end of the
                # desired cover range, what was recaptured may not always
                # encompass the entire cover range, so allow some small
                # tolerance
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
