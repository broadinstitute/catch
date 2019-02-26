"""Tests for genome module.
"""

from collections import OrderedDict
import unittest

from catch import genome

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class TestGenome(unittest.TestCase):
    """Tests methods in the Genome class.
    """

    def test_type_error(self):
        with self.assertRaises(TypeError):
            genome.Genome.from_one_seq(['ATCG'])
        with self.assertRaises(TypeError):
            genome.Genome.from_chrs(
                OrderedDict([("chr1", 'ATCG'), ("chr2", ['A', 'T'])]))

    def test_divided_into_chrs(self):
        genome_one_seq = genome.Genome.from_one_seq('ATCGCGAT')
        genome_two_chrs = genome.Genome.from_chrs(
            OrderedDict([("chr1", 'TAT'), ("chr2", 'GCG')]))
        self.assertFalse(genome_one_seq.divided_into_chrs())
        self.assertTrue(genome_two_chrs.divided_into_chrs())

    def test_size(self):
        genome_one_seq = genome.Genome.from_one_seq('ATCGCNGAT')
        genome_two_chrs = genome.Genome.from_chrs(
            OrderedDict([("chr1", 'TAT'), ("chr2", 'GCNG')]))
        # count ambiguous bases ('N') too
        self.assertEqual(genome_one_seq.size(), 9)
        self.assertEqual(genome_two_chrs.size(), 7)
        # run again to ensure the cached value is correct
        self.assertEqual(genome_one_seq.size(), 9)
        self.assertEqual(genome_two_chrs.size(), 7)

    def test_size_unambig(self):
        genome_one_seq = genome.Genome.from_one_seq('ATCGCNGAT')
        genome_two_chrs = genome.Genome.from_chrs(
            OrderedDict([("chr1", 'TAT'), ("chr2", 'GCNG')]))
        # do not count ambiguous bases (i.e., not 'N')
        self.assertEqual(genome_one_seq.size(only_unambig=True), 8)
        self.assertEqual(genome_two_chrs.size(only_unambig=True), 6)
        # run again to ensure the cached value is correct
        self.assertEqual(genome_one_seq.size(only_unambig=True), 8)
        self.assertEqual(genome_two_chrs.size(only_unambig=True), 6)

    def test_hash(self):
        genome_one_seq = genome.Genome.from_one_seq('ATCGCNGAT')
        genome_two_chrs = genome.Genome.from_chrs(
            OrderedDict([("chr1", 'TAT'), ("chr2", 'GCNG')]))
        s = {genome_one_seq, genome_two_chrs}
        self.assertEqual(len(s), 2)
        self.assertIn(genome_one_seq, s)
        self.assertIn(genome_two_chrs, s)
        # run again to ensure the cached value is correct
        self.assertIn(genome_one_seq, s)
        self.assertIn(genome_two_chrs, s)

    def test_equals(self):
        genome_one_seq1 = genome.Genome.from_one_seq('ATCGCNGAT')
        genome_one_seq2 = genome.Genome.from_one_seq('ATCGCNGAA')
        genome_two_chrs1 = genome.Genome.from_chrs(
            OrderedDict([("chr1", 'TAT'), ("chr2", 'GCNG')]))
        genome_two_chrs2 = genome.Genome.from_chrs(
            OrderedDict([("chr1", 'TAT'), ("chr2", 'GCNA')]))
        self.assertEqual(genome_one_seq1, genome_one_seq1)
        self.assertNotEqual(genome_one_seq1, genome_one_seq2)
        self.assertNotEqual(genome_one_seq1, genome_two_chrs1)
        self.assertNotEqual(genome_one_seq1, genome_two_chrs2)
        self.assertEqual(genome_one_seq2, genome_one_seq2)
        self.assertNotEqual(genome_one_seq2, genome_two_chrs1)
        self.assertNotEqual(genome_one_seq2, genome_two_chrs2)
        self.assertEqual(genome_two_chrs1, genome_two_chrs1)
        self.assertNotEqual(genome_two_chrs1, genome_two_chrs2)

    def test_break_into_fragments_from_one_seq(self):
        genome_one = genome.Genome.from_one_seq('ATCGTTAA')
        broken = genome_one.break_into_fragments(3)
        expected_genome = genome.Genome.from_chrs(
                OrderedDict([('0', 'ATC'), ('1', 'GTT'), ('2', 'AA')]))
        self.assertEqual(broken, expected_genome)

    def test_break_into_fragments_from_one_seq_with_full_end(self):
        genome_one = genome.Genome.from_one_seq('ATCGTTAA')
        broken = genome_one.break_into_fragments(3, include_full_end=True)
        expected_genome = genome.Genome.from_chrs(
                OrderedDict([('0', 'ATC'), ('1', 'GTT'), ('2', 'TAA')]))
        self.assertEqual(broken, expected_genome)

    def test_break_into_fragments_from_chrs(self):
        genome_two_chrs = genome.Genome.from_chrs(
                OrderedDict([("chr1", 'ATCGTTAA'), ("chr2", 'AATTCCGGG')]))
        broken = genome_two_chrs.break_into_fragments(3)
        expected_genome = genome.Genome.from_chrs(
                OrderedDict([("chr1-0", 'ATC'), ("chr1-1", 'GTT'), ("chr1-2", 'AA'),
                             ("chr2-0", 'AAT'), ("chr2-1", 'TCC'), ("chr2-2", 'GGG')]))
        self.assertEqual(broken, expected_genome)
