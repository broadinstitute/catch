"""Tests for reverse_complement_filter module.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

import unittest

from hybseldesign import probe
from hybseldesign.filter import reverse_complement_filter


"""Tests the reverse complement filter output on contrived input.
"""
class TestReverseComplementFilter(unittest.TestCase):

  def test_basic(self):
    input = ['ATCGTCGCGG',
             'TGACGACACC',
             'AGCTGCTCTT',
             'AGCNGCTCTT']
    desired_output = ['ATCGTCGCGG',
                      'TGACGACACC',
                      'AGCTGCTCTT',
                      'AGCNGCTCTT',
                      'CCGCGACGAT',
                      'GGTGTCGTCA',
                      'AAGAGCAGCT',
                      'AAGAGCNGCT']
    input_probes = [probe.Probe.from_str(s) for s in input]
    desired_output_probes = [probe.Probe.from_str(s) for \
                              s in desired_output]
    f = reverse_complement_filter.ReverseComplementFilter()
    f.filter(input_probes)
    self.assertItemsEqual(f.input_probes, input_probes)
    self.assertItemsEqual(f.output_probes, desired_output_probes)
