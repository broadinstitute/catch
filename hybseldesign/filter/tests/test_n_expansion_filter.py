"""Tests for n_expansion_filter module.
"""

import unittest

from hybseldesign.filter import n_expansion_filter
from hybseldesign import probe

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class TestNExpansionFilter(unittest.TestCase):
    """Tests the N expansion filter output on contrived input.
    """

    def check_output(self, input, desired_output):
        input_probes = [probe.Probe.from_str(s) for s in input]
        desired_output_probes = [probe.Probe.from_str(s)
                                 for s in desired_output]
        f = n_expansion_filter.NExpansionFilter()
        f.filter(input_probes)
        self.assertEqual(f.input_probes, input_probes)
        self.assertEqual(f.output_probes, desired_output_probes)

    def test_no_N(self):
        input = ['ATCG', 'ACTG']
        desired_output = ['ATCG', 'ACTG']
        self.check_output(input, desired_output)

    def test_single_N_in_probe(self):
        input = ['ATCG', 'ANCG', 'ACTG']
        desired_output = ['ATCG', 'AACG', 'ATCG', 'ACCG', 'AGCG', 'ACTG']
        self.check_output(input, desired_output)

    def test_two_N_in_probe1(self):
        input = ['ATCG', 'ANNG', 'ACTG']
        desired_output = ['ATCG', 'AAAG', 'AATG', 'AACG', 'AAGG',
                          'ATAG', 'ATTG', 'ATCG', 'ATGG', 'ACAG',
                          'ACTG', 'ACCG', 'ACGG', 'AGAG', 'AGTG',
                          'AGCG', 'AGGG', 'ACTG']
        self.check_output(input, desired_output)

    def test_two_N_in_probe2(self):
        input = ['ATCG', 'ANCN', 'ACTG']
        desired_output = ['ATCG', 'AACA', 'AACT', 'AACC', 'AACG',
                          'ATCA', 'ATCT', 'ATCC', 'ATCG', 'ACCA',
                          'ACCT', 'ACCC', 'ACCG', 'AGCA', 'AGCT',
                          'AGCC', 'AGCG', 'ACTG']
        self.check_output(input, desired_output)

    def test_single_N_at_start(self):
        input = ['ATCG', 'NACG', 'ACTG']
        desired_output = ['ATCG', 'AACG', 'TACG', 'CACG', 'GACG', 'ACTG']
        self.check_output(input, desired_output)

    def test_N_in_two_probes(self):
        input = ['ATCG', 'ANCG', 'ACNG']
        desired_output = ['ATCG', 'AACG', 'ATCG', 'ACCG', 'AGCG',
                          'ACAG', 'ACTG', 'ACCG', 'ACGG']
        self.check_output(input, desired_output)
