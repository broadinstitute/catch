"""Tests for adapter_filter module.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

import unittest
import logging

from hybseldesign import probe
from hybseldesign.filter import adapter_filter as af
from hybseldesign.utils import interval


"""Tests the adapter filter output on contrived input.
"""
class TestAdapterFilter(unittest.TestCase):

  def setUp(self):
    # Disable logging
    logging.disable(logging.INFO)

  def get_filter_and_output(self, lcf_thres, mismatches,
      target_genomes, input, k, num_kmers_per_probe):
    input_probes = [probe.Probe.from_str(s) for s in input]
    f = af.AdapterFilter(mismatches=mismatches, lcf_thres=lcf_thres,
          kmer_size=k, num_kmers_per_probe=num_kmers_per_probe)
    f.target_genomes = target_genomes
    f.filter(input_probes)
    return (f, f.output_probes)

  def test_votes(self):
    target_genomes = [ 'ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEF',
                       'ZYXWVFGHIJWUTSOPQRSTFEDCBAZYXWVF' ]
    input = []
    # Create probes of length 6 bp with a stride of 3 bp
    for tg in target_genomes:
      input += [ tg[i:(i+6)] for i in xrange(0, len(tg)-6+1, 3) ]
    #f, output = self.get_filter_and_output(6, 0, target_genomes,
    #    input, 3, 10)

  def tearDown(self):
    # Re-enable logging
    logging.disable(logging.NOTSET)

