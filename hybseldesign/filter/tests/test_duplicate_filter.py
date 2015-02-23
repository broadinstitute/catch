"""Tests for duplicate_filter module.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

import unittest

from hybseldesign import probe
from hybseldesign.filter import duplicate_filter


"""Tests the duplicate filter output on contrived input.
"""
class TestDuplicateFilter(unittest.TestCase):

  def test_basic(self):
    input = ['ATCGTCGCGG',
             'ATCGTAGCGG',
             'ATCGTCACGG',
             'ATCGTAGCGG',
             'ATTGTCGCGG',
             'ATCGTCGCGG']
    desired_output = ['ATCGTCGCGG',
                      'ATCGTAGCGG',
                      'ATCGTCACGG',
                      'ATTGTCGCGG']
    input_probes = [probe.Probe.from_str(s) for s in input]
    desired_output_probes = [probe.Probe.from_str(s) for \
                              s in desired_output]
    f = duplicate_filter.DuplicateFilter()
    f.filter(input_probes)
    self.assertItemsEqual(f.input_probes, input_probes)
    self.assertItemsEqual(f.output_probes, desired_output_probes)
