"""Tests for polya_filter module.
"""

import random
import unittest

from catch.filter import polya_filter
from catch import probe

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class TestPolyAFilter(unittest.TestCase):
    """Tests the polyA filter output on contrived input.
    """

    def check_output(self, input, desired_output, length, mismatches):
        input_probes = [probe.Probe.from_str(s) for s in input]
        desired_output_probes = [probe.Probe.from_str(s)
                                 for s in desired_output]
        f = polya_filter.PolyAFilter(length, mismatches,
            min_exact_length_to_consider=2)
        f.filter(input_probes)
        self.assertEqual(f.input_probes, input_probes)
        self.assertEqual(f.output_probes, desired_output_probes)

    def test_no_polya(self):
        input = ['CCGGAAGGCC', 'GCGCGCGCGC']
        desired_output = ['CCGGAAGGCC', 'GCGCGCGCGC']
        self.check_output(input, desired_output, 4, 0)
        self.check_output(input, desired_output, 4, 1)

    def test_exact_polya(self):
        input = ['CCGGAAAAACC', 'CCCCCCCCCCC']
        desired_output = ['CCCCCCCCCCC']
        self.check_output(input, desired_output, 4, 0)
        self.check_output(input, desired_output, 4, 1)
        self.check_output(input, desired_output, 4, 2)
        self.check_output(input, desired_output, 5, 0)
        self.check_output(input, desired_output, 5, 1)
        self.check_output(input, desired_output, 5, 2)

    def test_mismatched_polya_0_mismatches(self):
        input = ['CCGGAAGAACC', 'CCCCCCCCCCC']
        desired_output = ['CCGGAAGAACC', 'CCCCCCCCCCC']
        self.check_output(input, desired_output, 4, 0)
        self.check_output(input, desired_output, 5, 0)

    def test_mismatched_polya_with_mismatches(self):
        input = ['CCGGAAGAACC', 'CCCCCCCCCCC']
        desired_output = ['CCCCCCCCCCC']
        self.check_output(input, desired_output, 4, 1)
        self.check_output(input, desired_output, 4, 2)
        self.check_output(input, desired_output, 5, 1)
        self.check_output(input, desired_output, 5, 2)

    def test_polyt(self):
        input = ['CCGGTTTTTCC', 'CCCCCCCCCCC']
        desired_output = ['CCCCCCCCCCC']
        self.check_output(input, desired_output, 4, 0)
        self.check_output(input, desired_output, 4, 1)
        self.check_output(input, desired_output, 4, 2)

    def test_different_length_probes_1(self):
        input = ['CCGGAAAAACC', 'CCCCC', 'CCGGTTTT', 'TTTG']
        desired_output = ['CCCCC', 'TTTG']
        self.check_output(input, desired_output, 4, 0)

    def test_different_length_probes_2(self):
        input = ['CCGGAAAAACC', 'CCCCC', 'CCGGTTTT', 'TTTG']
        desired_output = ['CCCCC']
        self.check_output(input, desired_output, 4, 1)
