"""Tests for seq_io module.
"""

from collections import OrderedDict
import logging
import tempfile
import unittest

from hybseldesign.datasets import ebola_zaire_with_2014
from hybseldesign.datasets import lassa
from hybseldesign import genome
from hybseldesign.utils import seq_io

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class TestEbolaZaireFASTARead(unittest.TestCase):
    """Tests reading the Ebola Zaire (w/ 2014) dataset (FASTA file).
    """

    def setUp(self):
        # Disable logging
        logging.disable(logging.INFO)

        assert len(ebola_zaire_with_2014.fasta_paths) == 1
        self.seqs_map = seq_io.read_fasta(ebola_zaire_with_2014.fasta_paths[0])
        self.seqs = list(self.seqs_map.values())

    def test_num_seqs(self):
        """Test that there are 883 sequences.
        """
        self.assertEqual(len(self.seqs), 883)

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
        assert len(ebola_zaire_with_2014.fasta_paths) == 1
        for seq in seq_io.iterate_fasta(ebola_zaire_with_2014.fasta_paths[0]):
            generator_seqs += [seq]
        self.assertEqual(generator_seqs, self.seqs)

    def tearDown(self):
        # Re-enable logging
        logging.disable(logging.NOTSET)


class TestDatasetGenomeRead(unittest.TestCase):
    """Tests reading a dataset.
    """

    def setUp(self):
        # Disable logging
        logging.disable(logging.INFO)

    def test_single_chr_dataset(self):
        """Tests that the genomes obtained from reading the
        ebola_zaire_with_2014 dataset are the same as those obtained
        from directly reading the FASTA.

        This is effectively executing most of the same code as
        seq_io.read_dataset_genomes() but does check that it correctly
        enters the condition of reading just one sequence per genome.
        """
        genomes = seq_io.read_dataset_genomes(ebola_zaire_with_2014)
        assert len(ebola_zaire_with_2014.fasta_paths) == 1
        desired_genomes = [
            genome.Genome.from_one_seq(s)
            for s in seq_io.read_fasta(ebola_zaire_with_2014.fasta_paths[0]).\
                values()
        ]
        self.assertEqual(genomes, desired_genomes)

    def test_multi_chr_dataset(self):
        """Tests that the lassa dataset can be read.

        This does not test that the genomes are read correctly -- just
        that they can be read without issues.
        """
        genomes = seq_io.read_dataset_genomes(lassa)

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
