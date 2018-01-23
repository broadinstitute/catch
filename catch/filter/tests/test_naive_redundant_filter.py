"""Tests for naive_redundant_filter module.
"""

import logging
import unittest

import numpy as np

from catch.filter import naive_redundant_filter as nrf
from catch import probe

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class TestNaiveRedundantFilterShiftAndMismatchCount(unittest.TestCase):
    """Tests the shift/mismatch method for determining if probes are redundant.
    """

    def setUp(self):
        # Disable logging
        logging.disable(logging.INFO)

    def compare_input_with_desired_output(self, shift, mismatches, input,
                                          desired_output):
        input_probes = [probe.Probe.from_str(s) for s in input]
        desired_output_probes = [probe.Probe.from_str(s)
                                 for s in desired_output]
        f = nrf.NaiveRedundantFilter(
            nrf.redundant_shift_and_mismatch_count(shift, mismatches))
        f.filter(input_probes)
        self.assertCountEqual(f.input_probes, input_probes)
        self.assertCountEqual(f.output_probes, desired_output_probes)

    def test_no_shift_no_mismatch(self):
        # This should detect and remove duplicates
        input = ['ATCGTCGCGG', 'ATCGTAGCGG', 'ATCGTCACGG', 'ATCGTAGCGG',
                 'ATTGTCGCGG', 'ATCGTCGCGG']
        desired_output = ['ATCGTCGCGG', 'ATCGTAGCGG', 'ATCGTCACGG',
                          'ATTGTCGCGG']
        self.compare_input_with_desired_output(0, 0, input, desired_output)

    def test_no_shift_one_mismatch(self):
        input = ['ATCGTCGCGG', 'ATCGTAGCGG', 'ATCGTCACGG', 'ATCGTAGCGG',
                 'ATTGTCGCGA', 'ATCGTCGCGG']
        desired_output = ['ATCGTCGCGG', 'ATTGTCGCGA']
        self.compare_input_with_desired_output(0, 1, input, desired_output)

    def test_one_shift_one_mismatch(self):
        input = ['ATCGTCGCGG', 'TCGTAGCGGA', 'TATCGTCACG', 'AGTCGTAGCG',
                 'ATTGTCGCGA', 'ATCGTCGCGG']
        desired_output = ['ATCGTCGCGG', 'AGTCGTAGCG', 'ATTGTCGCGA']
        self.compare_input_with_desired_output(1, 1, input, desired_output)

    def tearDown(self):
        # Re-enable logging
        logging.disable(logging.NOTSET)


class TestRedundantByLongestCommonSubstring(unittest.TestCase):
    """Tests the redundant_longest_common_substring function.
    """

    def setUp(self):
        # Set a seed because probe.shares_some_kmers uses randomness
        np.random.seed(1)
        # Disable logging
        logging.disable(logging.INFO)

    def test_not_redundant_without_heuristic(self):
        fn = nrf.redundant_longest_common_substring(1, 5,
            prune_with_heuristic_and_anchor=False)
        a = probe.Probe.from_str('ATCGATCGATCG')
        b = probe.Probe.from_str('CAGGCCGGCTGA')
        self.assertFalse(fn(a, b))

    def test_is_redundant_without_heuristic(self):
        fn = nrf.redundant_longest_common_substring(2, 5,
            prune_with_heuristic_and_anchor=False)
        a = probe.Probe.from_str('ATCGATCGATCG')
        b = probe.Probe.from_str('CGATAGCTAGAT')
        self.assertTrue(fn(a, b))

    def test_not_redundant_with_heuristic(self):
        # Sequences need to be long because probe.shares_some_kmers
        # will use a default k-mer size of 20
        fn = nrf.redundant_longest_common_substring(2, 25,
            prune_with_heuristic_and_anchor=True)

        # Test when they do not share k-mers
        a = probe.Probe.from_str('ATCGATCGATCGAAAAAAAAAAAAAAATTTTT')
        b = probe.Probe.from_str('CCGGCCGGCCGGTTTTTTTTTTTTTTTAAAAA')
        self.assertFalse(fn(a, b))

        # Test when they do share k-mers but are still not
        # redundant
        a = probe.Probe.from_str('ATCGATCGATCGAAAAAAAAAAAAAAATTTTT')
        b = probe.Probe.from_str('CCGGCCGGCCGGAAAAAAAAAAAAAAATTTTT')
        self.assertFalse(fn(a, b))

        # Test when they do share k-mers but are still not
        # redundant
        fn = nrf.redundant_longest_common_substring(2, 60,
            prune_with_heuristic_and_anchor=True)
        a = probe.Probe.from_str(('A'*40 + 'DEF')*4)
        b = probe.Probe.from_str(('A'*40 + 'XYZ')*4)
        self.assertFalse(fn(a, b))

        # Test when they do share k-mers but are still not
        # redundant
        a = probe.Probe.from_str('ATCG' + 'AATT'*10 + 'CC' + 'GGCC'*4 + 'AT' + 'ATCG')
        b = probe.Probe.from_str('AATT'*10 + 'GG' + 'GGCC'*4 + 'CG' + 'CCCC' + 'CCGG')
        self.assertFalse(fn(a, b))

    def test_is_redundant_with_heuristic(self):
        # Sequences need to be long because probe.shares_some_kmers
        # will use a default k-mer size of 20
        fn = nrf.redundant_longest_common_substring(2, 25,
            prune_with_heuristic_and_anchor=True)

        # Test when they do share k-mers and are redundant
        a = probe.Probe.from_str('ATCGATCGATCG' + 'A'*50)
        b = probe.Probe.from_str('CGATAGCTAGAT' + 'A'*50)
        self.assertTrue(fn(a, b))

        # Test when they do share k-mers and are redundant
        fn = nrf.redundant_longest_common_substring(2, 60,
            prune_with_heuristic_and_anchor=True)
        a = probe.Probe.from_str('ATCG' + 'AATT'*10 + 'CC' + 'GGCC'*4 + 'AT' + 'ATCG')
        b = probe.Probe.from_str('AATT'*10 + 'GG' + 'GGCC'*4 + 'AT' + 'CCCC' + 'CCGG')
        self.assertTrue(fn(a, b))

    def tearDown(self):
        # Re-enable logging
        logging.disable(logging.NOTSET)
