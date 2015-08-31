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
        self.analyzer.run(window_length=6, window_stride=3)

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
        self.assertCountEqual(self.analyzer.target_covers[0][0][False],
                              [(0, 6), (4, 10), (14, 20), (17, 23)])
        self.assertCountEqual(self.analyzer.target_covers[0][0][True],
                              [(6, 12), (13, 19), (17, 23)])

        # genome_b
        self.assertCountEqual(self.analyzer.target_covers[1][0][False],
                              [(8, 14)])    # coverage in chr2
        self.assertCountEqual(self.analyzer.target_covers[1][0][True],
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

    def test_sliding_coverage(self):
        """Check the calculation of a coverage across sliding windows.

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
        # in provided sequence, coverage across the bases is (underneath
        # is average across a 6 bp window, sliding by 3 bp):
        # A T C C A T C C A T N G G G T T T G A A G C G
        # 1 1 1 1 2 2 1 1 1 1 0 0 0 0 1 1 1 2 2 2 1 1 1
        # -----|-----|-----|-----|-----|-----|---|-----
        #     8/6   8/6   4/6   2/6   5/6   9/6 9/6
        #      3     6     9     12    15   18   20
        expected = {3: 8/6., 6: 8/6., 9: 4/6., 12: 2/6., 15: 5/6.,
                    18: 9/6., 20: 9/6.}
        self.assertEqual(self.analyzer.sliding_coverage[0][0][False],
                         expected)
        # in reverse complement, coverage across the bases is (underneath
        # is average across a 6 bp window, sliding by 3 bp):
        # C G C T T C A A A C C C N A T G G A T G G A T
        # 0 0 0 0 0 0 1 1 1 1 1 1 0 1 1 1 1 2 2 1 1 1 1
        # -----|-----|-----|-----|-----|-----|---|-----
        #      0    1/2    1    5/6    1    8/6 8/6
        #      3     6     9     12    15   18   20
        expected = {3: 0, 6: 1/2., 9: 1, 12: 5/6., 15: 1, 18: 8/6.,
                    20: 8/6.}
        self.assertEqual(self.analyzer.sliding_coverage[0][0][True],
                         expected)

        # genome_b
        # in provided sequence, coverage across the bases is (underneath
        # is average across a 6 bp window, sliding by 3 bp):
        # C C C C C C ; N T G A A G C G
        # 0 0 0 0 0 0 ; 0 0 1 1 1 1 1 1
        # -----|----- | -----|---|-----
        #      0     1/6    4/6  1
        #      3      6      9   11
        # (The window centered at 6 is meaningless because it is between
        # chromosomes; it's an artifact of the fact that the chromosomes
        # are concatenated for this analysis.)
        expected = {3: 0, 6: 1/6., 9: 4/6., 11: 1}
        self.assertEqual(self.analyzer.sliding_coverage[1][0][False],
                         expected)
        # in reverse complement, coverage across the bases is (underneath
        # is average across a 6 bp window, sliding by 3 bp):
        # G G G G G G ; C G C T T C A N
        # 0 0 0 0 0 0 ; 0 0 0 0 0 0 0 0
        # -----|----- | -----|---|-----
        #      0      0      0   0
        #      3      6      9   11
        # (The window centered at 6 is meaningless because it is between
        # chromosomes; it's an artifact of the fact that the chromosomes
        # are concatenated for this analysis.)
        expected = {3: 0, 6: 0, 9: 0, 11: 0}
        self.assertEqual(self.analyzer.sliding_coverage[1][0][True],
                         expected)

    def test_data_matrix_string(self):
        """Test the data matrix generated.
        """
        data = self.analyzer._make_data_matrix_string()
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


class TestAnalyzerCoversWithCoverExtension(unittest.TestCase):
    """Tests the probe covers found when extending a probe's coverage.
    """

    def setUp(self):
        # Disable logging
        logging.disable(logging.INFO)

        # Create Analyzer instance with two target genomes
        genome_a = genome.Genome.from_one_seq('ATCCATCCATNGGGTTTGAAGCG')
        genome_b = genome.Genome.from_chrs(OrderedDict([('chr1', 'CCCCCCA'),
                                                        ('chr2', 'ANTGAAGCG')]))
        probes_str = ['ATCCAT', 'TTTGAA', 'GAAGCG', 'ATGGAT', 
                      'CCCCCC', 'AAACCC']
        probes = [probe.Probe.from_str(p) for p in probes_str]
        self.analyzer = ca.Analyzer(probes,
                                    [[genome_a], [genome_b]],
                                    target_genomes_names=["g_a", "g_b"],
                                    mismatches=0,
                                    lcf_thres=6,
                                    cover_extension=2,
                                    kmer_probe_map_k=3)
        self.analyzer.run(window_length=6, window_stride=3)

    def test_probe_cover_ranges(self):
        """Test the probe cover ranges that are found.

        Check in given sequence and in its reverse complement.
        """
        # genome_a
        self.assertCountEqual(self.analyzer.target_covers[0][0][False],
                              [(0, 8), (2, 12), (12, 22), (15, 23)])
        self.assertCountEqual(self.analyzer.target_covers[0][0][True],
                              [(4, 14), (11, 21), (15, 23)])

        # genome_b
        self.assertCountEqual(self.analyzer.target_covers[1][0][False],
                              [(0, 7), (8, 16)])    # coverage in chr1 and chr2
        self.assertCountEqual(self.analyzer.target_covers[1][0][True],
                              [])   # no coverage in reverse complement
