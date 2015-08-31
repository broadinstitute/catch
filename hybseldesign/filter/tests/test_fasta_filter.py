"""Tests for fasta_filter module.
"""

import logging
import tempfile
import unittest

from hybseldesign.filter import fasta_filter as ff
from hybseldesign import probe

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class TestFastaFilter(unittest.TestCase):
    """Tests the fasta filter output on contrived input.
    """

    def setUp(self):
        # Disable logging
        logging.disable(logging.WARNING)

    def test_basic(self):
        fasta_file = tempfile.NamedTemporaryFile(mode='w')
        fasta_file.write(">probe1\n")
        fasta_file.write("ATCGATCG\n")
        fasta_file.write(">probe2\n")
        fasta_file.write("GGGGG\n")
        fasta_file.seek(0)

        p1 = probe.Probe.from_str('CCCCC')
        p2 = probe.Probe.from_str('GGGGG')
        p3 = probe.Probe.from_str('ATCGTTTT')
        p4 = probe.Probe.from_str('ATCGATCG')
        p5 = probe.Probe.from_str('NNNNATCG')
        input_probes = [p1, p2, p3, p4, p5]

        fasta_filter = ff.FastaFilter(fasta_file.name)
        fasta_filter.filter(input_probes)
        self.assertEqual(fasta_filter.output_probes, [p4, p2])

        fasta_file.close()

    def test_skip_reverse_complements(self):
        fasta_file = tempfile.NamedTemporaryFile(mode='w')
        fasta_file.write(">probe1\n")
        fasta_file.write("ATCGATCG\n")
        fasta_file.write(">probe2 | reverse complement of probe1\n")
        fasta_file.write("CGATCGAT\n")
        fasta_file.write(">probe3\n")
        fasta_file.write("GGGGG\n")
        fasta_file.seek(0)

        p1 = probe.Probe.from_str('CCCCC')
        p2 = probe.Probe.from_str('GGGGG')
        p3 = probe.Probe.from_str('ATCGTTTT')
        p4 = probe.Probe.from_str('ATCGATCG')
        p5 = probe.Probe.from_str('CGATCGAT')
        p6 = probe.Probe.from_str('NNNNATCG')
        input_probes = [p1, p2, p3, p4, p5, p6]

        fasta_filter = ff.FastaFilter(fasta_file.name,
                                      skip_reverse_complements=True)
        fasta_filter.filter(input_probes)
        self.assertEqual(fasta_filter.output_probes, [p4, p2])

        fasta_file.close()

    def tearDown(self):
        # Re-enable logging 
        logging.disable(logging.NOTSET)
