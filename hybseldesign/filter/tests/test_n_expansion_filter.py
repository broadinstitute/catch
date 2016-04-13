"""Tests for n_expansion_filter module.
"""

import random
import unittest

from hybseldesign.filter import n_expansion_filter
from hybseldesign import probe

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class TestNExpansionFilter(unittest.TestCase):
    """Tests the N expansion filter output on contrived input.
    """

    def setUp(self):
        random.seed(1)

    def check_output(self, input, desired_output,
            limit_n_expansion_randomly=3):
        input_probes = [probe.Probe.from_str(s) for s in input]
        desired_output_probes = [probe.Probe.from_str(s)
                                 for s in desired_output]
        f = n_expansion_filter.NExpansionFilter(
                limit_n_expansion_randomly=limit_n_expansion_randomly)
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

    def test_limit_expansion_None(self):
        input = ['ATCG', 'ANCG', 'ACNG']
        desired_output = ['ATCG', 'AACG', 'ATCG', 'ACCG', 'AGCG',
                          'ACAG', 'ACTG', 'ACCG', 'ACGG']
        self.check_output(input, desired_output,
                          limit_n_expansion_randomly=None)

    def test_limit_expansion_0(self):
        input = ['ATCG', 'ANCG', 'ACNG']
        desired_output = ['ATCG', 'AACG', 'ACAG']
        self.check_output(input, desired_output,
                          limit_n_expansion_randomly=0)

    def test_limit_expansion_0_with_two_N(self):
        input = ['ATCG', 'ANNG', 'ACTG']
        desired_output = ['ATCG', 'AAAG', 'ACTG']
        self.check_output(input, desired_output,
                          limit_n_expansion_randomly=0)

    def test_limit_expansion_1_with_single_N(self):
        input = ['ATCG', 'ANCG', 'ACTG']
        desired_output = ['ATCG', 'AACG', 'ATCG', 'ACCG', 'AGCG', 'ACTG']
        self.check_output(input, desired_output,
                          limit_n_expansion_randomly=1)

    def test_limit_expansion_1_with_two_N(self):
        input = ['ATCG', 'ANNG', 'ACTG']
        desired_output = ['ATCG', 'AAAG', 'AATG', 'AACG', 'AAGG', 'ACTG']
        self.check_output(input, desired_output,
                          limit_n_expansion_randomly=1)

    def test_limit_expansion_1_with_all_N(self):
        input = ['ATCG', 'NNNN', 'ACTG']
        desired_output = ['ATCG', 'AAAG', 'TAAG', 'CAAG', 'GAAG', 'ACTG']
        self.check_output(input, desired_output,
                          limit_n_expansion_randomly=1)

    def test_limit_expansion_2_with_all_N(self):
        input = ['ATCG', 'NNNN', 'ACTG']
        desired_output = ['ATCG',
                          'AAAA', 'AAAT', 'AAAC', 'AAAG',
                          'TAAA', 'TAAT', 'TAAC', 'TAAG',
                          'CAAA', 'CAAT', 'CAAC', 'CAAG',
                          'GAAA', 'GAAT', 'GAAC', 'GAAG',
                          'ACTG']
        self.check_output(input, desired_output,
                          limit_n_expansion_randomly=2)

    def test_limit_expansion_2_with_two_N(self):
        input = ['ATCG', 'ANNG', 'ACTG']
        desired_output = ['ATCG', 'AAAG', 'AATG', 'AACG', 'AAGG',
                          'ATAG', 'ATTG', 'ATCG', 'ATGG', 'ACAG',
                          'ACTG', 'ACCG', 'ACGG', 'AGAG', 'AGTG',
                          'AGCG', 'AGGG', 'ACTG']
        self.check_output(input, desired_output,
                          limit_n_expansion_randomly=2)

    def test_limit_expansion_4_with_two_N(self):
        input = ['ATCG', 'ANNG', 'ACTG']
        desired_output = ['ATCG', 'AAAG', 'AATG', 'AACG', 'AAGG',
                          'ATAG', 'ATTG', 'ATCG', 'ATGG', 'ACAG',
                          'ACTG', 'ACCG', 'ACGG', 'AGAG', 'AGTG',
                          'AGCG', 'AGGG', 'ACTG']
        self.check_output(input, desired_output,
                          limit_n_expansion_randomly=4)

    def test_limit_expansion_1_with_N_in_two_probes(self):
        input = ['ATCG', 'ANCG', 'ACNG']
        desired_output = ['ATCG', 'AACG', 'ATCG', 'ACCG', 'AGCG',
                          'ACAG', 'ACTG', 'ACCG', 'ACGG']
        self.check_output(input, desired_output,
                          limit_n_expansion_randomly=1)
