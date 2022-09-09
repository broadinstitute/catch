"""Tests for probe_designer module.
"""

import logging
import unittest

from catch.filter import duplicate_filter
from catch.filter import probe_designer
from catch import genome
from catch import probe

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class TestProbeDesigner(unittest.TestCase):
    """Tests the probe designer output on contrived input.
    """

    def setUp(self):
        # Disable logging
        logging.disable(logging.WARNING)

    def test_one_filter1(self):
        """A basic test with a duplicate filter and one input sequence.
        Note that this test is dependent on the default values for
        generating candidate probes: probe length of 100 bp with a stride
        of 50 bp.
        """
        seqs = [
            [genome.Genome.from_one_seq('A' * 100 + 'B' * 100 + 'A' * 100)]]
        desired_candidate_probes = \
            ['A' * 100, 'A' * 50 + 'B' * 50, 'B' * 100, 'B' * 50 + 'A' * 50,
             'A' * 100]
        desired_candidate_probes = \
            [probe.Probe.from_str(s) for s in desired_candidate_probes]
        desired_final_probes = ['A' * 100, 'A' * 50 + 'B' * 50, 'B' * 100,
                                'B' * 50 + 'A' * 50]
        desired_final_probes = \
            [probe.Probe.from_str(s) for s in desired_final_probes]
        df = duplicate_filter.DuplicateFilter()
        pb = probe_designer.ProbeDesigner(seqs, [df], probe_length=100,
            probe_stride=50)
        pb.design()
        self.assertCountEqual(pb.candidate_probes, desired_candidate_probes)
        self.assertCountEqual(pb.final_probes, desired_final_probes)

    def test_one_filter2(self):
        """A basic test with a duplicate filter and one input sequence.
        Note that this test uses a probe length of 75 bp and a stride of
        25 bp.
        """
        seqs = [
            [genome.Genome.from_one_seq('A' * 100 + 'B' * 100 + 'A' * 100)]]
        desired_candidate_probes = \
            ['A' * 75, 'A' * 75, 'A' * 50 + 'B' * 25, 'A' * 25 + 'B' * 50,
             'B' * 75, 'B' * 75, 'B' * 50 + 'A' * 25, 'B' * 25 + 'A' * 50,
             'A' * 75, 'A' * 75]
        desired_candidate_probes = \
            [probe.Probe.from_str(s) for s in desired_candidate_probes]
        desired_final_probes = ['A' * 75, 'A' * 50 + 'B' * 25,
                                'A' * 25 + 'B' * 50, 'B' * 75,
                                'B' * 50 + 'A' * 25, 'B' * 25 + 'A' * 50]
        desired_final_probes = \
            [probe.Probe.from_str(s) for s in desired_final_probes]
        df = duplicate_filter.DuplicateFilter()
        pb = probe_designer.ProbeDesigner(seqs, [df], probe_length=75,
            probe_stride=25)
        pb.design()
        self.assertCountEqual(pb.candidate_probes, desired_candidate_probes)
        self.assertCountEqual(pb.final_probes, desired_final_probes)

    def test_two_groupings(self):
        """Tests two groupings of input sequences in which the first
        grouping has two sequences and the second grouping has one
        sequence.
        Note that this test is dependent on the default values for
        generating candidate probes: probe length of 100 bp with a stride
        of 50 bp.
        """
        seqs = [[genome.Genome.from_one_seq('A' * 200),
                 genome.Genome.from_one_seq('B' * 150)],
                [genome.Genome.from_one_seq('C' * 300)]]
        desired_candidate_probes = \
            ['A' * 100, 'A' * 100, 'A' * 100, 'B' * 100, 'B' * 100,
             'C' * 100, 'C' * 100, 'C' * 100, 'C' * 100, 'C' * 100]
        desired_candidate_probes = \
            [probe.Probe.from_str(s) for s in desired_candidate_probes]
        desired_final_probes = ['A' * 100, 'B' * 100, 'C' * 100]
        desired_final_probes = \
            [probe.Probe.from_str(s) for s in desired_final_probes]
        df = duplicate_filter.DuplicateFilter()
        pb = probe_designer.ProbeDesigner(seqs, [df], probe_length=100,
            probe_stride=50)
        pb.design()
        self.assertCountEqual(pb.candidate_probes, desired_candidate_probes)
        self.assertCountEqual(pb.final_probes, desired_final_probes)

    def test_with_small_sequences(self):
        """A test with a duplicate filter and input sequences that are smaller
        than the probe length.
        """
        seqs = [[genome.Genome.from_one_seq('ABCDEFGHIJKLMN'),
                 genome.Genome.from_one_seq('ABCDEFGHIXKLMN'),
                 genome.Genome.from_one_seq('XYZAB')]]
        desired_candidate_probes = \
            ['ABCDEF', 'DEFGHI', 'GHIJKL', 'IJKLMN',
             'ABCDEF', 'DEFGHI', 'GHIXKL', 'IXKLMN',
             'XYZAB']
        desired_candidate_probes = \
            [probe.Probe.from_str(s) for s in desired_candidate_probes]
        desired_final_probes = ['ABCDEF', 'DEFGHI', 'GHIJKL', 'IJKLMN',
                                'GHIXKL', 'IXKLMN', 'XYZAB']
        desired_final_probes = \
            [probe.Probe.from_str(s) for s in desired_final_probes]
        df = duplicate_filter.DuplicateFilter()
        pb = probe_designer.ProbeDesigner(seqs, [df], probe_length=6,
                                          probe_stride=3,
                                          allow_small_seqs=5)
        pb.design()
        self.assertCountEqual(pb.candidate_probes, desired_candidate_probes)
        self.assertCountEqual(pb.final_probes, desired_final_probes)

    def tearDown(self):
        # Re-enable logging
        logging.disable(logging.NOTSET)


class TestProbeDesignerWithClustering(unittest.TestCase):
    """Tests the probe designer output, including cluster sequences, on
    contrived input.
    """

    def setUp(self):
        # Disable logging
        logging.disable(logging.WARNING)

    def test_cluster_genomes(self):
        """Test the _cluster_genomes() function.
        """
        g1 = genome.Genome.from_one_seq('ATTA'*500)
        g2 = genome.Genome.from_one_seq('TAAT'*500)
        g3 = genome.Genome.from_one_seq('GC'*500)
        g4 = genome.Genome.from_one_seq('CG'*500)
        seqs = [[g1, g3], [g2, g4]]
        pb = probe_designer.ProbeDesigner(seqs, [], probe_length=100,
            probe_stride=50, cluster_threshold=0.1, cluster_merge_after=None,
            cluster_method='simple')
        clustered_genomes = pb._cluster_genomes()

        # g1 and g2 should be in a cluster, and g3 and g4 should be in a
        # cluster
        self.assertEqual(len(clustered_genomes), 2)
        cluster1 = clustered_genomes[0]
        cluster2 = clustered_genomes[1]
        if g1 in cluster1:
            self.assertEqual(set(cluster1), set([g1, g2]))
            self.assertEqual(set(cluster2), set([g3, g4]))
        else:
            self.assertEqual(set(cluster1), set([g3, g4]))
            self.assertEqual(set(cluster2), set([g1, g2]))

    def test_cluster_genomes_with_fragmenting(self):
        """Test the cluster_genomes() function, with fragmenting sequences.
        """
        g1 = genome.Genome.from_one_seq('ATTA'*500)
        g2 = genome.Genome.from_one_seq('TAAT'*500)
        g3 = genome.Genome.from_one_seq('GC'*500)
        g4 = genome.Genome.from_one_seq('CG'*500)
        seqs = [[g1, g3], [g2, g4]]
        pb = probe_designer.ProbeDesigner(seqs, [], probe_length=100,
            probe_stride=50, cluster_threshold=0.1, cluster_merge_after=None,
            cluster_method='simple', cluster_fragment_length=500)
        clustered_genomes = pb._cluster_genomes()

        g1_frag = genome.Genome.from_one_seq('ATTA'*125)
        g2_frag = genome.Genome.from_one_seq('TAAT'*125)
        g3_frag = genome.Genome.from_one_seq('GC'*250)
        g4_frag = genome.Genome.from_one_seq('CG'*250)

        # g1_frag and g2_frag should be in a cluster (4 of each),
        # and g3_frag and g4_frag should be in a cluster (2 of each)
        self.assertEqual(len(clustered_genomes), 2)
        cluster1 = clustered_genomes[0]
        cluster2 = clustered_genomes[1]
        if g1_frag in cluster1:
            self.assertCountEqual(cluster1, [g1_frag]*4 + [g2_frag]*4)
            self.assertCountEqual(cluster2, [g3_frag]*2 + [g4_frag]*2)
        else:
            self.assertCountEqual(cluster1, [g3_frag]*2 + [g4_frag]*2)
            self.assertCountEqual(cluster2, [g1_frag]*4 + [g2_frag]*4)

    def test_filter_with_clusters(self):
        """Tests running two back-to-back DuplicateFilters, with the first
        run on each cluster and the second run on the merge of the output
        probes from the clusters.
        """
        g1 = genome.Genome.from_one_seq('ATTA'*500)
        g2 = genome.Genome.from_one_seq('TAAT'*500)
        g3 = genome.Genome.from_one_seq('GC'*500)
        g4 = genome.Genome.from_one_seq('CG'*500)
        seqs = [[g1, g3], [g2, g4]]

        df1 = duplicate_filter.DuplicateFilter()
        df2 = duplicate_filter.DuplicateFilter()
        pb = probe_designer.ProbeDesigner(seqs, [df1, df2], probe_length=100,
            probe_stride=8, cluster_threshold=0.1, cluster_merge_after=df1,
            cluster_method='simple')
        pb.design()

        desired_candidate_probes = ['ATTA'*25, 'TAAT'*25, 'CG'*50, 'GC'*50]
        desired_candidate_probes = \
            [probe.Probe.from_str(s) for s in desired_candidate_probes]
        desired_final_probes = ['A' * 75, 'A' * 50 + 'B' * 25,
                                'A' * 25 + 'B' * 50, 'B' * 75,
                                'B' * 50 + 'A' * 25, 'B' * 25 + 'A' * 50]
        desired_final_probes = ['ATTA'*25, 'TAAT'*25, 'CG'*50, 'GC'*50]
        desired_final_probes = \
            [probe.Probe.from_str(s) for s in desired_final_probes]

        self.assertEqual(set(pb.candidate_probes), set(desired_candidate_probes))
        self.assertCountEqual(pb.final_probes, desired_final_probes)

    def tearDown(self):
        # Re-enable logging
        logging.disable(logging.NOTSET)
