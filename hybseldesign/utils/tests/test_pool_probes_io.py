"""Tests for pool_probes_io module.
"""

import logging
import tempfile
import unittest

from hybseldesign.utils import pool_probes_io

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class TestReadProbeCountTable(unittest.TestCase):
    """Tests reading a probe count table from a file.
    """

    def setUp(self):
        # Disable logging
        logging.disable(logging.INFO)

    def test_valid_file(self):
        # Write a temporary tsv file
        self.tsv = tempfile.NamedTemporaryFile(mode='w')
        self.tsv.write("dataset\tmismatches\tcover_extension\tnum_probes\n")
        self.tsv.write("ebola\t1\t10\t1000\n")
        self.tsv.write("ebola\t2\t30\t500\n")
        self.tsv.write("zika\t3\t20\t200\n")
        self.tsv.write("zika\t4\t40\t100\n")
        self.tsv.seek(0)

        param_names, probe_counts = pool_probes_io.read_table_of_probe_counts(
            self.tsv.name)
        self.assertEqual(param_names, ('mismatches', 'cover_extension'))
        self.assertEqual(probe_counts["ebola"][(2, 30)], 500)
        self.assertEqual(probe_counts["zika"][(3, 20)], 200)

        self.tsv.close()

    def test_invalid_file(self):
        # Write a temporary tsv file
        self.tsv = tempfile.NamedTemporaryFile(mode='w')
        self.tsv.write("dataset\tmismatches\tnum_probes\tcover_extension\n")
        self.tsv.write("ebola\t1\t1000\t10\n")
        self.tsv.seek(0)

        # The header must end with 'num_probes', so reading it should
        # raise an exception
        with self.assertRaises(Exception):
            probe_counts = pool_probes_io.read_table_of_probe_counts(self.tsv.name)

        self.tsv.close()

    def tearDown(self):
        # Re-enable logging
        logging.disable(logging.NOTSET)


class TestReadDatasetWeightsTable(unittest.TestCase):
    """Tests reading a dataset weights table from a file.
    """

    def setUp(self):
        # Disable logging
        logging.disable(logging.INFO)

    def test_valid_file(self):
        # Write a temporary tsv file
        self.tsv = tempfile.NamedTemporaryFile(mode='w')
        self.tsv.write("dataset\tweight\n")
        self.tsv.write("ebola\t1.5\n")
        self.tsv.write("zika\t1\n")
        self.tsv.seek(0)

        weights = pool_probes_io.read_table_of_dataset_weights(self.tsv.name,
            set(['ebola', 'zika']))
        self.assertEqual(weights, {'ebola': 1.5, 'zika': 1.0})

        self.tsv.close()

    def test_file_missing_dataset(self):
        # Write a temporary tsv file
        self.tsv = tempfile.NamedTemporaryFile(mode='w')
        self.tsv.write("dataset\tweight\n")
        self.tsv.write("ebola\t1.5\n")
        self.tsv.write("zika\t1\n")
        self.tsv.seek(0)

        # We will check for 'lassa', which is not in the table, so
        # reading it should raise an exception
        with self.assertRaises(Exception):
            weights = pool_probes_io.read_table_of_dataset_weights(self.tsv.name,
                set('ebola', 'zika', 'lassa'))

        self.tsv.close()

    def tearDown(self):
        # Re-enable logging
        logging.disable(logging.NOTSET)
