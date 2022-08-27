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

        f = ndf.NearDuplicateFilterWithHammingDistance(2, 10)
        f.k = 3
        output_probes = f.filter(input_probes)
        self.assertEqual(len(output_probes), 1)
        self.assertIn(output_probes[0], input_probes)

    def test_all_similar_with_exact_dup(self):
        input = ['ATCGTCGCGG', 'ATCGTCGCGG', 'ATCGTGGCGG', 'TTCGTCGCGG',
                 'ATCGGCGCGG']
        input_probes = [probe.Probe.from_str(s) for s in input]

        f = ndf.NearDuplicateFilterWithHammingDistance(2, 10)
        f.k = 3
        output_probes = f.filter(input_probes)
        self.assertEqual(len(output_probes), 1)
        # The first probe in input_probes is the most common, so this
        # should be the one that is kept
        self.assertEqual(output_probes[0], input_probes[0])

    def test_all_similar_but_zero_dist_thres(self):
        input = ['ATCGTCGCGG', 'ATCGTGGCGG', 'TTCGTCGCGG', 'ATCGGCGCGG']
        input_probes = [probe.Probe.from_str(s) for s in input]

        f = ndf.NearDuplicateFilterWithHammingDistance(0, 10)
        f.k = 3
        output_probes = f.filter(input_probes)
        self.assertCountEqual(output_probes, input_probes)

    def test_all_similar_but_one_too_far(self):
        input = ['ATCGTCGCGG', 'ATCGTGGCGG', 'TTCGTCGCGG', 'ATCGGCGCCT']
        input_probes = [probe.Probe.from_str(s) for s in input]

        f = ndf.NearDuplicateFilterWithHammingDistance(2, 10)
        f.k = 3
        output_probes = f.filter(input_probes)
        self.assertEqual(len(output_probes), 2)
        # The last probe in input_probes is barely >2 mismatches
        # from the others
        self.assertIn(input_probes[-1], output_probes)

    def test_two_clusters(self):
        cluster1 = ['ATCGTCGCGG', 'ATCGTGGCGG', 'TTCGTCGCGG', 'ATCGGCGCGG']
        cluster2 = ['GGCTTACTGA', 'GGCTTACTGA', 'GGCTTTCTGA', 'GGCTTACTAT']
        input = cluster1 + cluster2
        random.shuffle(input)
        input_probes = [probe.Probe.from_str(s) for s in input]

        f = ndf.NearDuplicateFilterWithHammingDistance(2, 10)
        f.k = 3
        f.reporting_prob = 0.90
        output_probes = f.filter(input_probes)
        self.assertEqual(len(output_probes), 2)
        self.assertTrue((output_probes[0].seq_str in cluster1 and
                         output_probes[1].seq_str in cluster2) or
                        (output_probes[0].seq_str in cluster2 and
                         output_probes[1].seq_str in cluster1))

    def test_parallelization_over_groups(self):
        group1_cluster1 = ['ATCGTCGCGG',
                           'ATCGTGGCGG',
                           'TTCGTCGCGG',
                           'ATCGGCGCGG']
        group1_cluster2 = ['GGCTTACTGA',
                           'GGCTTACTGA',
                           'GGCTTTCTGA',
                           'GGCTTACTAT']
        group1 = group1_cluster1 + group1_cluster2
        random.shuffle(group1)
        group1_probes = [probe.Probe.from_str(s) for s in group1]

        group2_cluster1 = ['ATATATATAT',
                           'ATATCGATAT']
        group2_cluster2 = ['CGCGCGCGCG',
                           'CGCGCGATCG']
        group2 = group2_cluster1 + group2_cluster2
        random.shuffle(group2)
        group2_probes = [probe.Probe.from_str(s) for s in group2]

        input_probes_grouped = [group1_probes, group2_probes]
        group_clusters = {0: {'cluster1': group1_cluster1,
                              'cluster2': group1_cluster2},
                          1: {'cluster1': group2_cluster1,
                              'cluster2': group2_cluster2}}

        f = ndf.NearDuplicateFilterWithHammingDistance(5, 10)
        f.k = 3
        f.reporting_prob = 0.95
        for np in [None, 1, 2, 4, 8]:
            output_probes = f.filter(input_probes_grouped,
                    input_is_grouped=True,
                    num_processes=np)

            # Check every grouping
            for i in [0, 1]:
                # Every grouping had 2 clusters, so the output should
                # have 2 probes each from one of the clusters
                self.assertEqual(len(output_probes[i]), 2)
                group_cluster1 = group_clusters[i]['cluster1']
                group_cluster2 = group_clusters[i]['cluster2']
                self.assertTrue((output_probes[i][0].seq_str in group_cluster1 and
                                 output_probes[i][1].seq_str in group_cluster2) or
                                (output_probes[i][0].seq_str in group_cluster2 and
                                 output_probes[i][1].seq_str in group_cluster1))


class TestNearDuplicateFilterWithMinHash(unittest.TestCase):
    """Tests output of near duplicate filter using MinHash.
    """

    def setUp(self):
        # Set a random seed so hash functions are always the same
        random.seed(0)

    def test_all_similar_no_exact_dup(self):
        input = ['ATCGTCGCGG', 'ATCGTGGCGG', 'TTCGTCGCGG', 'ATCGGCGCGG']
        input_probes = [probe.Probe.from_str(s) for s in input]

        f = ndf.NearDuplicateFilterWithMinHash(0.8, 3)
        f.k = 3
        f.reporting_prob = 0.90
        output_probes = f.filter(input_probes)
        self.assertEqual(len(output_probes), 1)
        self.assertIn(output_probes[0], input_probes)

    def test_all_similar_with_exact_dup(self):
        input = ['ATCGTCGCGG', 'ATCGTCGCGG', 'ATCGTGGCGG', 'TTCGTCGCGG',
                 'ATCGGCGCGG']
        input_probes = [probe.Probe.from_str(s) for s in input]

        f = ndf.NearDuplicateFilterWithMinHash(0.8, 3)
        f.k = 3
        f.reporting_prob = 0.90
        output_probes = f.filter(input_probes)
        self.assertEqual(len(output_probes), 1)
        # The first probe in input_probes is the most common, so this
        # should be the one that is kept
        self.assertEqual(output_probes[0], input_probes[0])

    def test_all_similar_but_zero_dist_thres(self):
        input = ['ATCGTCGCGG', 'ATCGTGGCGG', 'TTCGTCGCGG', 'ATCGGCGCGG']
        input_probes = [probe.Probe.from_str(s) for s in input]

        f = ndf.NearDuplicateFilterWithMinHash(0, 3)
        f.k = 3
        f.reporting_prob = 0.90
        output_probes = f.filter(input_probes)
        self.assertCountEqual(output_probes, input_probes)

    def test_all_similar_but_one_too_far(self):
        input = ['ATCGTCGCGG', 'ATCGTGGCGG', 'TTCGTCGCGG', 'ATTGGGGCCA']
        input_probes = [probe.Probe.from_str(s) for s in input]

        f = ndf.NearDuplicateFilterWithMinHash(0.8, 3)
        f.k = 3
        f.reporting_prob = 0.90
        output_probes = f.filter(input_probes)
        self.assertEqual(len(output_probes), 2)
        # The last probe in input_probes is far from the others
        self.assertIn(input_probes[-1], output_probes)

    def test_two_clusters(self):
        cluster1 = ['ATCGTCGCGG', 'ATCGTGGCGG', 'TTCGTCGCGG', 'ATCGGCGCGG']
        cluster2 = ['GGCTTACTGA', 'GGCTTACTGA', 'GGCTTTCTGA', 'GGCTTACTAT']
        input = cluster1 + cluster2
        random.shuffle(input)
        input_probes = [probe.Probe.from_str(s) for s in input]

        f = ndf.NearDuplicateFilterWithMinHash(0.8, 3)
        f.k = 3
        f.reporting_prob = 0.90
        output_probes = f.filter(input_probes)
        self.assertEqual(len(output_probes), 2)
        self.assertTrue((output_probes[0].seq_str in cluster1 and
                         output_probes[1].seq_str in cluster2) or
                        (output_probes[0].seq_str in cluster2 and
                         output_probes[1].seq_str in cluster1))

    def test_parallelization_over_groups(self):
        group1_cluster1 = ['ATCGTCGCGG',
                           'ATCGTGGCGG',
                           'TTCGTCGCGG',
                           'ATCGGCGCGG']
        group1_cluster2 = ['GGCTTACTGA',
                           'GGCTTACTGA',
                           'GGCTTTCTGA',
                           'GGCTTACTAT']
        group1 = group1_cluster1 + group1_cluster2
        random.shuffle(group1)
        group1_probes = [probe.Probe.from_str(s) for s in group1]

        group2_cluster1 = ['ACCCGAG']
        group2_cluster2 = ['GGGGGGG',
                           'TGGGGGG']
        group2 = group2_cluster1 + group2_cluster2
        random.shuffle(group2)
        group2_probes = [probe.Probe.from_str(s) for s in group2]

        group3_cluster1 = ['ACCCGAA',
                           'ACCCGTA']
        group3_cluster2 = ['CGTGTGA',
                           'CGTGAGA',
                           'CGTGTGA',
                           'CGTGCGA',
                           'CGTGTCA',
                           'CGTGTAA',
                           'CGTGAGA',
                           'AGTGTGA']
        group3 = group3_cluster1 + group3_cluster2
        random.shuffle(group3)
        group3_probes = [probe.Probe.from_str(s) for s in group3]

        group4_cluster1 = ['ATATATATAT',
                           'ATATCGATAT']
        group4_cluster2 = ['CGCGCGCGCG',
                           'CGCGCGATCG']
        group4 = group4_cluster1 + group4_cluster2
        random.shuffle(group4)
        group4_probes = [probe.Probe.from_str(s) for s in group4]

        input_probes_grouped = [group1_probes, group2_probes,
                group3_probes, group4_probes]
        group_clusters = {0: {'cluster1': group1_cluster1,
                              'cluster2': group1_cluster2},
                          1: {'cluster1': group2_cluster1,
                              'cluster2': group2_cluster2},
                          2: {'cluster1': group3_cluster1,
                              'cluster2': group3_cluster2},
                          3: {'cluster1': group4_cluster1,
                              'cluster2': group4_cluster2}}

        f = ndf.NearDuplicateFilterWithMinHash(0.9, 3)
        f.k = 3
        f.reporting_prob = 0.95
        for np in [None, 1, 2, 4, 8]:
            output_probes = f.filter(input_probes_grouped,
                    input_is_grouped=True,
                    num_processes=np)

            # Check every grouping
            for i in [0, 1, 2, 3]:
                # Every grouping had 2 clusters, so the output should
                # have 2 probes each from one of the clusters
                self.assertEqual(len(output_probes[i]), 2)
                group_cluster1 = group_clusters[i]['cluster1']
                group_cluster2 = group_clusters[i]['cluster2']
                self.assertTrue((output_probes[i][0].seq_str in group_cluster1 and
                                 output_probes[i][1].seq_str in group_cluster2) or
                                (output_probes[i][0].seq_str in group_cluster2 and
                                 output_probes[i][1].seq_str in group_cluster1))

