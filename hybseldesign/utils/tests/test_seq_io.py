"""Tests for seq_io module.
"""

import logging
import unittest

from hybseldesign.datasets import ebola2014
from hybseldesign import genome
from hybseldesign.utils import seq_io

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class TestEbola2014FASTARead(unittest.TestCase):

  """Tests reading the Ebola 2014 dataset (FASTA file).
  """

  def setUp(self):
    # Disable logging
    logging.disable(logging.INFO)

    self.seqs_map = seq_io.read_fasta(ebola2014.fasta_path)
    self.seqs = self.seqs_map.values()

  def test_num_seqs(self):
    """Test that there are 99 sequences.
    """
    self.assertEqual(len(self.seqs), 99)

  def test_seq_length(self):
    """Test that all sequences are of length 18-19 kbp.
    """
    for seq in self.seqs:
      self.assertGreater(len(seq), 18000)
      self.assertLess(len(seq), 19000)

  def test_seq_content(self):
    """Test that all sequences contain valid characters.

    They should only contain 'A', 'T', 'C', 'G', and 'N'.
    """
    for seq in self.seqs:
      n = sum([seq.count(char) for char in ['A','T','C','G','N']])
      self.assertEqual(len(seq), n)

  def test_generator(self):
    """Test that the generator works correctly.

    Tests this by comparing to the output from read_fasta.
    """
    generator_seqs = []
    for seq in seq_io.iterate_fasta(ebola2014.fasta_path):
      generator_seqs += [seq]
    self.assertEqual(generator_seqs, self.seqs)

  def tearDown(self):
    # Re-enable logging
    logging.disable(logging.NOTSET)


class TestDatasetGenomeRead(unittest.TestCase):
  
  def setUp(self):
    # Disable logging
    logging.disable(logging.INFO)

  def test_single_chr_dataset(self):
    """Tests that the genomes obtained from reading the ebola2014
    dataset are the same as those obtained from directly reading the
    FASTA.

    This is effectively executing most of the same code as
    seq_io.read_dataset_genomes() but does check that it correctly
    enters the condition of reading just one sequence per genome.
    """
    genomes = seq_io.read_dataset_genomes(ebola2014)
    desired_genomes = [ genome.Genome.from_one_seq(s) for s in \
                          seq_io.read_fasta(ebola2014.fasta_path).\
                          values() ]
    self.assertEqual(genomes, desired_genomes)

  def tearDown(self):
    # Re-enable logging
    logging.disable(logging.NOTSET)

