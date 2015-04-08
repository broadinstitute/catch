"""Tests for dominating_set_filter module.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

import unittest
import logging

from hybseldesign import probe
from hybseldesign.filter import dominating_set_filter as dsf
from hybseldesign.filter import naive_redundant_filter as nrf


class TestDominatingSetFilter(unittest.TestCase):

  """Tests the dominating set filter output on contrived input.
  """

  def setUp(self):
    # Disable logging
    logging.disable(logging.INFO)

  def compare_input_with_desired_output(self, shift, mismatches,
      input, desired_output):
    input_probes = [probe.Probe.from_str(s) for s in input]
    desired_output_probes = [probe.Probe.from_str(s) for \
                              s in desired_output]
    f = dsf.DominatingSetFilter(
         nrf.redundant_shift_and_mismatch_count(shift, mismatches))
    f.filter(input_probes)
    self.assertItemsEqual(f.input_probes, input_probes)
    self.assertItemsEqual(f.output_probes, desired_output_probes)

  def test_one_shift_one_mismatch(self):
    input = ['ATCGTCGCGG',
             'TCGTAGCGGA',
             'CGTAGCGGAT']
    desired_output = ['TCGTAGCGGA']
    self.compare_input_with_desired_output(1, 1, input, desired_output)

  def tearDown(self):
    # Re-enable logging
    logging.disable(logging.NOTSET)

