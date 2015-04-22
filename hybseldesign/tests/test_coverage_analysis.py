"""Tests for coverage_analysis module.
"""

from collections import OrderedDict
import logging
import unittest

from hybseldesign import coverage_analysis as ca
from hybseldesign import genome
from hybseldesign import probe

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class TestAnalyzerWithTwoTargetGenomes(unittest.TestCase):
    """Tests methods in the Analyzer class.
    """

    def setUp(self):
        # Disable logging
        logging.disable(logging.INFO)

        # Create Analyzer instance with two target genomes
        genome_a = genome.Genome.from_one_seq('ATCCATCCATGGGTTTGAAGCG')
        genome_b = genome.Genome.from_chrs(OrderedDict([('chr1', 'CCCCCC'),
                                                        ('chr2', 'TGAAGCG')]))
        probes_str = ['ATCCAT', 'TTTGAA', 'GAAGCG', 'ATGGAT', 'AAACCC']
        probes = [probe.Probe.from_str(p) for p in probes_str]
        self.analyzer = ca.Analyzer(probes,
                                    [[genome_a], [genome_b]],
                                    target_genomes_names=["g_a", "g_b"],
                                    mismatches=0,
                                    lcf_thres=6,
                                    kmer_size=3,
                                    num_kmers_per_probe=10)
        self.analyzer.run()

    def test_probe_cover_ranges(self):
        """Test the probe cover ranges that are found.

        Check in given sequence and in its reverse complement.
        """
        self.assertEqual(len(self.analyzer.target_covers), 2)  # two groupings
        self.assertEqual(len(self.analyzer.target_covers[0]), 1)  # one genome
        self.assertEqual(len(self.analyzer.target_covers[1]), 1)  # one genome
        self.assertEqual(len(self.analyzer.target_covers[0][0]), 2)  # False/True
        self.assertEqual(len(self.analyzer.target_covers[1][0]), 2)  # False/True

        # genome_a
        self.assertItemsEqual(self.analyzer.target_covers[0][0][False],
                              [(0, 6), (4, 10), (13, 19), (16, 22)])
        self.assertItemsEqual(self.analyzer.target_covers[0][0][True],
                              [(6, 12), (12, 18), (16, 22)])

        # genome_b
        self.assertItemsEqual(self.analyzer.target_covers[1][0][False],
                              [(7, 13)])    # coverage in chr2
        self.assertItemsEqual(self.analyzer.target_covers[1][0][True],
                              [])   # no coverage in reverse complement

    def test_bp_covered(self):
        """Test the calculation of the number of bp covered.

        Check in given sequence and its reverse complement.
        """
        self.assertEqual(len(self.analyzer.bp_covered), 2)  # two groupings
        self.assertEqual(len(self.analyzer.bp_covered[0]), 1)  # one genome
        self.assertEqual(len(self.analyzer.bp_covered[1]), 1)  # one genome
        self.assertEqual(len(self.analyzer.bp_covered[0][0]), 2)  # False/True
        self.assertEqual(len(self.analyzer.bp_covered[1][0]), 2)  # False/True

        # genome_a
        # GGG is left out in the provided sequence (so 22-3=19 is covered)
        self.assertEqual(self.analyzer.bp_covered[0][0][False], 19)
        # CGCTTA is left out of the reverse complement (so 22-6=16 is covered)
        self.assertEqual(self.analyzer.bp_covered[0][0][True], 16)

        # genome_b
        # chr2 has coverage
        self.assertEqual(self.analyzer.bp_covered[1][0][False], 6)
        # no coverage in reverse complement
        self.assertEqual(self.analyzer.bp_covered[1][0][True], 0)

    def test_average_coverage(self):
        """Test the calculation of the average coverage.

        Check in the given sequence and its reverse complement.
        """
        self.assertEqual(len(self.analyzer.average_coverage), 2)  # two groupings
        self.assertEqual(len(self.analyzer.average_coverage[0]), 1)  # one genome
        self.assertEqual(len(self.analyzer.average_coverage[1]), 1)  # one genome
        self.assertEqual(len(self.analyzer.average_coverage[0][0]), 2)  # False/True
        self.assertEqual(len(self.analyzer.average_coverage[1][0]), 2)  # False/True

        # genome_a
        # in provided sequence, coverage across the bases is:
        # A T C C A T C C A T G G G T T T G A A G C G
        # 1 1 1 1 2 2 1 1 1 1 0 0 0 1 1 1 2 2 2 1 1 1
        # average = 24/22
        self.assertEqual(self.analyzer.average_coverage[0][0][False], 24./22)
        # in reverse complement, coverage across the bases is:
        # C G C T T C A A A C C C A T G G A T G G A T
        # 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 2 2 1 1 1 1
        # average = 18/22
        self.assertEqual(self.analyzer.average_coverage[0][0][True], 18./22)

        # genome_b
        # in provided sequence, coverage across the bases is:
        # C C C C C C ; T G A A G C G
        # 0 0 0 0 0 0 ; 0 1 1 1 1 1 1
        # average = 6/13
        self.assertEqual(self.analyzer.average_coverage[1][0][False], 6./13)
        # in reverse complement, coverage across the bases is:
        # G G G G G G ; C G C T T C A
        # 0 0 0 0 0 0 ; 0 0 0 0 0 0 0
        # average = 0/13
        self.assertEqual(self.analyzer.average_coverage[1][0][True], 0./13)

    def test_data_matrix(self):
        """Test the data matrix generated.
        """
        data = self.analyzer._make_data_matrix()
        expected = [["Genome",
                     "Num bases covered",
                     "Average coverage/depth"],
                    ["g_a, genome 0", "19 (86.36%)", "1.09"],
                    ["g_a, genome 0 (rc)", "16 (72.73%)", "0.82"],
                    ["g_b, genome 0", "6 (46.15%)", "0.46"],
                    ["g_b, genome 0 (rc)", "0 (<0.01%)", "<0.01"]]
        self.assertEqual(data, expected)

    def tearDown(self):
        # Re-enable logging
        logging.disable(logging.NOTSET)
