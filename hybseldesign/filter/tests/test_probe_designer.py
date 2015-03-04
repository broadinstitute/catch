"""Tests for probe_designer module.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

import unittest

from hybseldesign import probe
from hybseldesign.filter import probe_designer
from hybseldesign.filter import duplicate_filter


"""Tests the probe designer output on contrived input.
"""
class TestProbeDesigner(unittest.TestCase):

  """A basic test with a duplicate filter and one input sequence.
  Note that this test is dependent on the default values for
  generating candidate probes: probe length of 100 bp with a stride
  of 50 bp.
  """
  def test_one_filter(self):
    seqs = [ 'A'*100 + 'B'*100 + 'A'*100 ]
    desired_candidate_probes = \
        [ 'A'*100, 'A'*50+'B'*50, 'B'*100, 'B'*50+'A'*50,
          'A'*100 ]
    desired_candidate_probes = \
        [probe.Probe.from_str(s) for s in desired_candidate_probes]
    desired_final_probes = [ 'A'*100, 'A'*50+'B'*50, 'B'*100,
                             'B'*50+'A'*50 ]
    desired_final_probes = \
        [probe.Probe.from_str(s) for s in desired_final_probes]
    df = duplicate_filter.DuplicateFilter()
    pb = probe_designer.ProbeDesigner(seqs, [df])
    pb.design()
    self.assertEqual(pb.candidate_probes, desired_candidate_probes)
    self.assertEqual(pb.final_probes, desired_final_probes)
