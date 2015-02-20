"""Tests for naive_redundant_filter module.
"""

# Author: Hayden Metsky <hayden@mit.edu>

import unittest
import logging

from hybseldesign import probe
from hybseldesign.filter import naive_redundant_filter as nrf


"""Tests the shift and mismatch method for determining whether two
probes are redundant.
"""
class TestNaiveRedundantFilterShiftAndMismatchCount(unittest.TestCase):

  def setUp(self):
    # Disable logging
    logging.disable(logging.INFO)

  def compare_input_with_desired_output(self, shift, mismatches,
      input, desired_output):
    input_probes = [probe.Probe.from_str(s) for s in input]
    desired_output_probes = [probe.Probe.from_str(s) for \
                              s in desired_output]
    f = nrf.NaiveRedundantFilter(
         nrf.redundant_shift_and_mismatch_count(shift, mismatches))
    f.filter(input_probes)
    self.assertItemsEqual(f.input_probes, input_probes)
    self.assertItemsEqual(f.output_probes, desired_output_probes)

  def test_no_shift_no_mismatch(self):
    # This should detect and remove duplicates
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
    self.compare_input_with_desired_output(0, 0, input, desired_output)

  def test_no_shift_one_mismatch(self):
    input = ['ATCGTCGCGG',
             'ATCGTAGCGG',
             'ATCGTCACGG',
             'ATCGTAGCGG',
             'ATTGTCGCGA',
             'ATCGTCGCGG']
    desired_output = ['ATCGTCGCGG',
                      'ATTGTCGCGA']
    self.compare_input_with_desired_output(0, 1, input, desired_output)

  def test_one_shift_one_mismatch(self):
    input = ['ATCGTCGCGG',
             'TCGTAGCGGA',
             'TATCGTCACG',
             'AGTCGTAGCG',
             'ATTGTCGCGA',
             'ATCGTCGCGG']
    desired_output = ['ATCGTCGCGG',
                      'AGTCGTAGCG',
                      'ATTGTCGCGA']
    self.compare_input_with_desired_output(1, 1, input, desired_output)

  def tearDown(self):
    # Re-enable logging
    logging.disable(logging.NOTSET)

