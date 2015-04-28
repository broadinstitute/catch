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
        genome_a = genome.Genome.from_one_seq('ATCCATCCATNGGGTTTGAAGCG')
        genome_b = genome.Genome.from_chrs(OrderedDict([('chr1', 'CCCCCC'),
                                                        ('chr2', 'NTGAAGCG')]))
        probes_str = ['ATCCAT', 'TTTGAA', 'GAAGCG', 'ATGGAT', 'AAACCC']
        probes = [probe.Probe.from_str(p) for p in probes_str]
        self.analyzer = ca.Analyzer(probes,
                                    [[genome_a], [genome_b]],
                                    target_genomes_names=["g_a", "g_b"],
                                    mismatches=0,
                                    lcf_thres=6,
                                    kmer_probe_map_k=3)
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
                              [(0, 6), (4, 10), (14, 20), (17, 23)])
        self.assertItemsEqual(self.analyzer.target_covers[0][0][True],
                              [(6, 12), (13, 19), (17, 23)])

        # genome_b
        self.assertItemsEqual(self.analyzer.target_covers[1][0][False],
                              [(8, 14)])    # coverage in chr2
        self.assertItemsEqual(self.analyzer.target_covers[1][0][True],
                              [])   # no coverage in reverse complement

    def test_bp_covered(self):
        """Test the calculation of the number of bp covered.

        Check in given sequence and its reverse complement.
        """
        # two groupings
        self.assertEqual(len(self.analyzer.bp_covered), 2)
        # one genome per grouping
        self.assertEqual(len(self.analyzer.bp_covered[0]), 1)
        self.assertEqual(len(self.analyzer.bp_covered[1]), 1)
        # False/True for each genome (reverse complement)
        self.assertEqual(len(self.analyzer.bp_covered[0][0]), 2)
        self.assertEqual(len(self.analyzer.bp_covered[1][0]), 2)

        # genome_a
        # NGGG is left out in the provided sequence (so 23-4=19 is covered)
        self.assertEqual(self.analyzer.bp_covered[0][0][False], 19)
        # CGCTTA and N is left out of the reverse complement
        # (so 23-7=16 is covered)
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
        # two groupings
        self.assertEqual(len(self.analyzer.average_coverage), 2)
        # one genome per grouping
        self.assertEqual(len(self.analyzer.average_coverage[0]), 1)
        self.assertEqual(len(self.analyzer.average_coverage[1]), 1)
        # False/True for each genome (reverse complement)
        self.assertEqual(len(self.analyzer.average_coverage[0][0]), 2)
        self.assertEqual(len(self.analyzer.average_coverage[1][0]), 2)

        # genome_a
        # in provided sequence, coverage across the bases is:
        # A T C C A T C C A T N G G G T T T G A A G C G
        # 1 1 1 1 2 2 1 1 1 1 0 0 0 0 1 1 1 2 2 2 1 1 1
        # average = 24/23 over all (24/22 over unambig)
        self.assertEqual(self.analyzer.average_coverage[0][0][False][0],
                         24. / 23)
        self.assertEqual(self.analyzer.average_coverage[0][0][False][1],
                         24. / 22)
        # in reverse complement, coverage across the bases is:
        # C G C T T C A A A C C C N A T G G A T G G A T
        # 0 0 0 0 0 0 1 1 1 1 1 1 0 1 1 1 1 2 2 1 1 1 1
        # average = 18/23 over all (18/22 over unambig)
        self.assertEqual(self.analyzer.average_coverage[0][0][True][0],
                         18. / 23)
        self.assertEqual(self.analyzer.average_coverage[0][0][True][1],
                         18. / 22)

        # genome_b
        # in provided sequence, coverage across the bases is:
        # C C C C C C ; N T G A A G C G
        # 0 0 0 0 0 0 ; 0 0 1 1 1 1 1 1
        # average = 6/14 over all (6/13 over unambig)
        self.assertEqual(self.analyzer.average_coverage[1][0][False][0],
                         6. / 14)
        self.assertEqual(self.analyzer.average_coverage[1][0][False][1],
                         6. / 13)
        # in reverse complement, coverage across the bases is:
        # G G G G G G ; C G C T T C A N
        # 0 0 0 0 0 0 ; 0 0 0 0 0 0 0 0
        # average = 0/13 over all (0/12 over unambig)
        self.assertEqual(self.analyzer.average_coverage[1][0][True][0],
                         0. / 13)
        self.assertEqual(self.analyzer.average_coverage[1][0][True][1],
                         0. / 12)

    def test_data_matrix(self):
        """Test the data matrix generated.
        """
        data = self.analyzer._make_data_matrix()
        expected = [["Genome",
                     "Num bases covered\n[over unambig]",
                     "Average coverage/depth\n[over unambig]"],
                    ["g_a, genome 0", "19 (82.61%) [86.36%]",
                     "1.04 [1.09]"],
                    ["g_a, genome 0 (rc)", "16 (69.57%) [72.73%]",
                     "0.78 [0.82]"],
                    ["g_b, genome 0", "6 (42.86%) [46.15%]",
                     "0.43 [0.46]"],
                    ["g_b, genome 0 (rc)", "0 (<0.01%) [<0.01%]",
                     "<0.01 [<0.01]"]]
        self.assertEqual(data, expected)

    def tearDown(self):
        # Re-enable logging
        logging.disable(logging.NOTSET)
