"""Tests for seq_io module.
"""

from collections import OrderedDict
import logging
import pathlib
import tempfile
import unittest

from catch import genome
from catch.utils import seq_io

__author__ = 'Hayden Metsky <hayden@mit.edu>'


ZAIRE_EBOLAVIRUS_PATH = str(pathlib.Path(__file__).parent.joinpath(
        'data/zaire_ebolavirus.fasta.gz'))


class TestEbolaZaireFASTARead(unittest.TestCase):
    """Tests reading the Ebola Zaire (w/ 2014) dataset (FASTA file).
    """

    def setUp(self):
        # Disable logging
        logging.disable(logging.INFO)

        self.seqs_map = seq_io.read_fasta(ZAIRE_EBOLAVIRUS_PATH)
        self.seqs = list(self.seqs_map.values())

    def test_num_seqs(self):
        """Test that there are 1525 sequences.
        """
        self.assertEqual(len(self.seqs), 1525)

    def test_seq_length(self):
        """Test that all sequences are of length 17-20 kbp.
        """
        for seq in self.seqs:
            self.assertGreater(len(seq), 17000)
            self.assertLess(len(seq), 20000)

    def test_seq_content(self):
        """Test that all sequences contain valid characters.

        They should only contain 'A', 'T', 'C', 'G', and 'N'.
        """
        for seq in self.seqs:
            n = sum([seq.count(char) for char in ['A', 'T', 'C', 'G', 'N']])
            self.assertEqual(len(seq), n)

    def test_generator(self):
        """Test that the generator works correctly.

        Tests this by comparing to the output from read_fasta.
        """
        generator_seqs = []
        for seq in seq_io.iterate_fasta(ZAIRE_EBOLAVIRUS_PATH):
            generator_seqs += [seq]
        self.assertEqual(generator_seqs, self.seqs)

    def tearDown(self):
        # Re-enable logging
        logging.disable(logging.NOTSET)


class TestFastaEmptyLineRead(unittest.TestCase):
    """Tests reading a fasta file with empty lines.

    This has previously led to errors being thrown.
    """

    def setUp(self):
        # Disable logging
        logging.disable(logging.INFO)

        # Write the temporary fasta file
        self.fasta = tempfile.NamedTemporaryFile(mode='w')
        self.fasta.write(">genome_1\n")
        self.fasta.write("ATACG\n")
        self.fasta.write("TATGC\n")
        self.fasta.write(">genome_2\n")
        self.fasta.write("ATCG\n")
        self.fasta.write("TT\n")
        self.fasta.write("GG\n")
        self.fasta.write("\n")
        self.fasta.write(">genome_3\n")
        self.fasta.write("AAA\n")
        self.fasta.write("CCC\n")
        self.fasta.write("\n")
        self.fasta.write("\n")
        self.fasta.write(">genome_4\n")
        self.fasta.write("ATA\n")
        self.fasta.write("CGC\n")
        self.fasta.write("\n")
        self.fasta.write("\n")
        self.fasta.write("\n")
        self.fasta.write(">genome_5\n")
        self.fasta.write("AGGA\n")
        self.fasta.write("CAAT\n")
        self.fasta.write("\n")
        self.fasta.write("\n")
        self.fasta.seek(0)

        self.expected = OrderedDict()
        self.expected["genome_1"] = "ATACGTATGC"
        self.expected["genome_2"] = "ATCGTTGG"
        self.expected["genome_3"] = "AAACCC"
        self.expected["genome_4"] = "ATACGC"
        self.expected["genome_5"] = "AGGACAAT"

    def test_read(self):
        seqs = seq_io.read_fasta(self.fasta.name)
        self.assertEqual(seqs, self.expected)

    def test_iterate(self):
        seqs = list(seq_io.iterate_fasta(self.fasta.name))
        self.assertEqual(seqs, list(self.expected.values()))

    def tearDown(self):
        self.fasta.close()

        # Re-enable logging
        logging.disable(logging.NOTSET)
