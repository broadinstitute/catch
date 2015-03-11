"""Tests for seq_io module.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

import unittest
import logging

from hybseldesign.utils import seq_io
from hybseldesign.datasets import ebola2014


"""Tests reading the Ebola 2014 dataset (FASTA file).
"""
class TestEbola2014FASTARead(unittest.TestCase):

  def setUp(self):
    # Disable logging
    logging.disable(logging.INFO)

    self.seqs_map = seq_io.read_fasta(ebola2014.fasta_path())
    self.seqs = self.seqs_map.values()

  """Test that there are 99 sequences.
  """
  def test_num_seqs(self):
    self.assertEqual(len(self.seqs), 99)

  """Test that all sequences are of length 18-19 kbp.
  """
  def test_seq_length(self):
    for seq in self.seqs:
      self.assertGreater(len(seq), 18000)
      self.assertLess(len(seq), 19000)

  """Test that all sequences contain only the characters
  'A', 'T', 'C', 'G', and 'N'.
  """
  def test_seq_content(self):
    for seq in self.seqs:
      n = sum([seq.count(char) for char in ['A','T','C','G','N']])
      self.assertEqual(len(seq), n)

  """Test that the generator works correctly by comparing to
  the output from read_fasta.
  """
  def test_generator(self):
    generator_seqs = []
    for seq in seq_io.iterate_fasta(ebola2014.fasta_path()):
      generator_seqs += [seq]
    self.assertEqual(generator_seqs, self.seqs)

  def tearDown(self):
    # Re-enable logging
    logging.disable(logging.NOTSET)

