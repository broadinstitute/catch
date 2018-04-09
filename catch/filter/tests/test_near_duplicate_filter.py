"""Tests for near_duplicate_filter module.
"""

import random
import unittest

from catch.filter import near_duplicate_filter as ndf
from catch import probe

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class TestNearDuplicateFilterWithHammingDistance(unittest.TestCase):
    """Tests output of near duplicate filter according to Hamming distance.
    """

    def setUp(self):
        # Set a random seed so hash functions are always the same
        random.seed(0)

    def test_all_similar_no_exact_dup(self):
        input = ['ATCGTCGCGG', 'ATCGTGGCGG', 'TTCGTCGCGG', 'ATCGGCGCGG']
        input_probes = [probe.Probe.from_str(s) for s in input]

        f = ndf.NearDuplicateFilterWithHammingDistance(2, 10, k=3)
        f.filter(input_probes)
        self.assertEqual(len(f.output_probes), 1)
        self.assertIn(f.output_probes[0], input_probes)

    def test_all_similar_with_exact_dup(self):
        input = ['ATCGTCGCGG', 'ATCGTCGCGG', 'ATCGTGGCGG', 'TTCGTCGCGG',
                 'ATCGGCGCGG']
        input_probes = [probe.Probe.from_str(s) for s in input]

        f = ndf.NearDuplicateFilterWithHammingDistance(2, 10, k=3)
        f.filter(input_probes)
        self.assertEqual(len(f.output_probes), 1)
        # The first probe in input_probes is the most common, so this
        # should be the one that is kept
        self.assertEqual(f.output_probes[0], input_probes[0])

    def test_all_similar_but_zero_dist_thres(self):
        input = ['ATCGTCGCGG', 'ATCGTGGCGG', 'TTCGTCGCGG', 'ATCGGCGCGG']
        input_probes = [probe.Probe.from_str(s) for s in input]

        f = ndf.NearDuplicateFilterWithHammingDistance(0, 10, k=3)
        f.filter(input_probes)
        self.assertCountEqual(f.output_probes, input_probes)

    def test_all_similar_but_one_too_far(self):
        input = ['ATCGTCGCGG', 'ATCGTGGCGG', 'TTCGTCGCGG', 'ATCGGCGCCT']
        input_probes = [probe.Probe.from_str(s) for s in input]

        f = ndf.NearDuplicateFilterWithHammingDistance(2, 10, k=3)
        f.filter(input_probes)
        self.assertEqual(len(f.output_probes), 2)
        # The last probe in input_probes is barely >2 mismatches
        # from the others
        self.assertIn(input_probes[-1], f.output_probes)

    def test_two_clusters(self):
        cluster1 = ['ATCGTCGCGG', 'ATCGTGGCGG', 'TTCGTCGCGG', 'ATCGGCGCGG']
        cluster2 = ['GGCTTACTGA', 'GGCTTACTGA', 'GGCTTTCTGA', 'GGCTTACTAT']
        input = cluster1 + cluster2
        random.shuffle(input)
        input_probes = [probe.Probe.from_str(s) for s in input]

        f = ndf.NearDuplicateFilterWithHammingDistance(2, 10, k=3)
        f.filter(input_probes)
        self.assertEqual(len(f.output_probes), 2)
        self.assertTrue((f.output_probes[0].seq_str in cluster1 and
                         f.output_probes[1].seq_str in cluster2) or
                        (f.output_probes[0].seq_str in cluster2 and
                         f.output_probes[1].seq_str in cluster1))

