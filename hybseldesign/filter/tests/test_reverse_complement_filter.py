"""Tests for reverse_complement_filter module.
"""

import unittest

from catch.filter import reverse_complement_filter
from catch import probe

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class TestReverseComplementFilter(unittest.TestCase):
    """Tests the reverse complement filter output on contrived input.
    """

    def test_basic(self):
        input = ['ATCGTCGCGG', 'TGACGACACC', 'AGCTGCTCTT', 'AGCNGCTCTT']
        desired_output = ['ATCGTCGCGG', 'TGACGACACC', 'AGCTGCTCTT',
                          'AGCNGCTCTT', 'CCGCGACGAT', 'GGTGTCGTCA',
                          'AAGAGCAGCT', 'AAGAGCNGCT']
        input_probes = [probe.Probe.from_str(s) for s in input]
        desired_output_probes = [probe.Probe.from_str(s)
                                 for s in desired_output]
        f = reverse_complement_filter.ReverseComplementFilter()
        f.filter(input_probes)
        self.assertCountEqual(f.input_probes, input_probes)
        self.assertCountEqual(f.output_probes, desired_output_probes)
