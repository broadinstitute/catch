"""Tests for candidate_probes module.
"""

import unittest

from hybseldesign.datasets import ebola_zaire_with_2014
from hybseldesign.filter import candidate_probes
from hybseldesign.utils import seq_io

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class TestCandidateProbesOnContrivedInput(unittest.TestCase):
    """Tests explicitly the generated candidate probes from contrived input.
    """

    def test_no_n(self):
        p = candidate_probes.make_candidate_probes_from_sequence(
            'ATCGTCGCGGATCG',
            probe_length=6,
            probe_stride=3,
            min_n_string_length=2)
        p = [''.join(x.seq) for x in p]
        self.assertCountEqual(p, ['ATCGTC', 'GTCGCG', 'GCGGAT', 'GGATCG'])

    def test_one_n(self):
        p = candidate_probes.make_candidate_probes_from_sequence(
            'ATCGNCGCGGATCG',
            probe_length=6,
            probe_stride=3,
            min_n_string_length=2)
        p = [''.join(x.seq) for x in p]
        self.assertCountEqual(p, ['ATCGNC', 'GNCGCG', 'GCGGAT', 'GGATCG'])

    def test_two_n(self):
        p = candidate_probes.make_candidate_probes_from_sequence(
            'ATNGNCGCGGATCG',
            probe_length=6,
            probe_stride=3,
            min_n_string_length=2)
        p = [''.join(x.seq) for x in p]
        self.assertCountEqual(p, ['ATNGNC', 'GNCGCG', 'GCGGAT', 'GGATCG'])

    def test_n_string1(self):
        p = candidate_probes.make_candidate_probes_from_sequence(
            'ATCGNCGNNTCG',
            probe_length=6,
            probe_stride=3,
            min_n_string_length=2)
        p = [''.join(x.seq) for x in p]
        self.assertCountEqual(p, ['ATCGNC', 'TCGNCG'])

    def test_n_string2(self):
        p = candidate_probes.make_candidate_probes_from_sequence(
            'ATCGNCGNNTCGATAT',
            probe_length=6,
            probe_stride=3,
            min_n_string_length=2)
        p = [''.join(x.seq) for x in p]
        self.assertCountEqual(p, ['ATCGNC', 'TCGNCG', 'TCGATA', 'TCGATA',
                                  'CGATAT'])

    def test_multiple_seqs(self):
        p = candidate_probes.make_candidate_probes_from_sequences(
            ['ATCGNCGNNTCG', 'ATCGNCGNNTCGATAT'],
            probe_length=6,
            probe_stride=3,
            min_n_string_length=2)
        p = [''.join(x.seq) for x in p]
        self.assertCountEqual(
            p, ['ATCGNC', 'TCGNCG'] + ['ATCGNC', 'TCGNCG', 'TCGATA', 'TCGATA',
                                       'CGATAT'])

    def test_buggy(self):
        """Tests the version with bugs, intended to replicate the original
        candidate probe designer.
        """
        p = candidate_probes.make_candidate_probes_from_sequences(
            ['ATCGNCGNNTCG', 'ATCGNCGNNTCGATAT'],
            probe_length=6,
            probe_stride=3,
            min_n_string_length=2,
            insert_bugs=True,
            move_all_n_string_flanking_probes_to_end=True)
        p = [''.join(x.seq) for x in p]
        # Use assertEqual rather than assertCountEqual to check order
        self.assertEqual(
            p,
            ['ATCGNC'] + ['ATCGNC', 'TCGATA'] + ['TCGNCG', 'TCGNCG', 'TCGATA'])


class TestCandidateProbesOnEbolaZaire(unittest.TestCase):
    """Tests the candidate probes from the Ebola Zaire (w/ 2014) dataset.
    """

    def setUp(self):
        seqs = [gnm.seqs[0]
                for gnm in seq_io.read_dataset_genomes(ebola_zaire_with_2014)]
        self.probes_100 = candidate_probes.make_candidate_probes_from_sequences(
            seqs,
            probe_length=100,
            probe_stride=50,
            min_n_string_length=2)
        self.probes_75 = candidate_probes.make_candidate_probes_from_sequences(
            seqs,
            probe_length=75,
            probe_stride=25,
            min_n_string_length=2)

    def test_probe_length(self):
        """Test that all probes are the correct length.
        """
        for probe in self.probes_100:
            self.assertEqual(len(probe.seq), 100)
        for probe in self.probes_75:
            self.assertEqual(len(probe.seq), 75)

    def test_probe_count(self):
        """Test that probe counts are in roughly a correct ratio.

        Since probes_75 has a stride of 25 and probes_100 has a stride of
        50, there should be roughly twice as many probes in probes_75 as
        in probes_100.
        """
        ratio = float(len(self.probes_75)) / len(self.probes_100)
        self.assertTrue(1.95 < ratio < 2.05)

    def test_n_string(self):
        """Test that no probes have a string of two or more 'N's.
        """
        for probe in self.probes_100:
            self.assertNotIn('NN', ''.join(probe.seq))
        for probe in self.probes_75:
            self.assertNotIn('NN', ''.join(probe.seq))
