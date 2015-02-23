"""Tests for candidate_probes module.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

import unittest

from hybseldesign.utils import seq_io                                  
from hybseldesign.datasets import ebola2014
from hybseldesign.filter import candidate_probes


"""Tests explicitly the generated candidate probes on contrived
inputs."""
class TestCandidateProbesOnContrivedInput(unittest.TestCase):

  def test_no_n(self):
    p = candidate_probes.make_candidate_probes_from_sequence(
              'ATCGTCGCGGATCG', probe_length=6, probe_stride=3,
              min_n_string_length=2)
    p = [x.seq.tostring() for x in p]
    self.assertItemsEqual(p,
        ['ATCGTC', 'GTCGCG', 'GCGGAT', 'GGATCG'])

  def test_one_n(self):
    p = candidate_probes.make_candidate_probes_from_sequence(
              'ATCGNCGCGGATCG', probe_length=6, probe_stride=3,
              min_n_string_length=2)
    p = [x.seq.tostring() for x in p]
    self.assertItemsEqual(p,                                          
        ['ATCGNC', 'GNCGCG', 'GCGGAT', 'GGATCG'])

  def test_two_n(self):
    p = candidate_probes.make_candidate_probes_from_sequence(         
              'ATNGNCGCGGATCG', probe_length=6, probe_stride=3,       
              min_n_string_length=2)
    p = [x.seq.tostring() for x in p]                                 
    self.assertItemsEqual(p,                                          
        ['ATNGNC', 'GNCGCG', 'GCGGAT', 'GGATCG']) 

  def test_n_string1(self):
    p = candidate_probes.make_candidate_probes_from_sequence(
              'ATCGNCGNNTCG', probe_length=6, probe_stride=3,
              min_n_string_length=2)
    p = [x.seq.tostring() for x in p]
    self.assertItemsEqual(p, ['ATCGNC', 'TCGNCG'])

  def test_n_string2(self):
    p = candidate_probes.make_candidate_probes_from_sequence(
              'ATCGNCGNNTCGATAT', probe_length=6, probe_stride=3,
              min_n_string_length=2)
    p = [x.seq.tostring() for x in p]
    self.assertItemsEqual(p,
        ['ATCGNC', 'TCGNCG', 'TCGATA', 'TCGATA', 'CGATAT'])
    
  def test_multiple_seqs(self):
    p = candidate_probes.make_candidate_probes_from_sequences(
              ['ATCGNCGNNTCG', 'ATCGNCGNNTCGATAT'],
              probe_length=6, probe_stride=3, min_n_string_length=2)
    p = [x.seq.tostring() for x in p]
    self.assertItemsEqual(p, ['ATCGNC', 'TCGNCG'] + \
        ['ATCGNC', 'TCGNCG', 'TCGATA', 'TCGATA', 'CGATAT'])


"""Tests the candidate probes from the Ebola 2014 dataset.
"""
class TestCandidateProbesOnEbola2014(unittest.TestCase):

  def setUp(self):
    seqs = seq_io.read_fasta(ebola2014.fasta_path()).values()
    self.probes = candidate_probes.\
        make_candidate_probes_from_sequences(seqs,
          probe_length=100, probe_stride=50,
          min_n_string_length=2)

  """Test that all probes are 100 bp.
  """
  def test_probe_length(self):
    for probe in self.probes:
      self.assertEqual(len(probe.seq), 100)

  """Test that no probes have a string of two or more 'N's.
  """
  def test_n_string(self):
    for probe in self.probes:
      self.assertNotIn('NN', probe.seq.tostring())
