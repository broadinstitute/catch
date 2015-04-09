"""Tests for adapter_filter module.
"""

import logging
import unittest

from hybseldesign.filter import adapter_filter as af
from hybseldesign.filter import candidate_probes as cp
from hybseldesign import genome
from hybseldesign import probe
from hybseldesign.utils import interval

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class TestAdapterFilter(unittest.TestCase):
    """Tests the adapter filter output on contrived input.
    """

    def setUp(self):
        # Disable logging
        logging.disable(logging.INFO)

    def get_filter_and_output(self, lcf_thres, mismatches, target_genomes,
                              input_probes, k, num_kmers_per_probe):
        f = af.AdapterFilter(mismatches=mismatches,
                             lcf_thres=lcf_thres,
                             kmer_size=k,
                             num_kmers_per_probe=num_kmers_per_probe)
        f.target_genomes = target_genomes
        f.filter(input_probes)
        return (f, f.output_probes)

    def assert_has_adapter(self, probe, adapter):
        """Assert that probe has a particular adapter.

        Args:
            probe: probe to check
            adapter: check if 'probe' has this adapter

        Returns:
            whether 'probe' starts and ends with either an 'A' or 'B'
            adapter, as specified by 'adapter'
        """
        if adapter == 'A':
            start = af.ADAPTER_A_5END
            end = af.ADAPTER_A_3END
        elif adapter == 'B':
            start = af.ADAPTER_B_5END
            end = af.ADAPTER_B_3END
        else:
            raise ValueError("Unknown adapter %s" % adapter)
        self.assertTrue(self.probe.seq_str.startswith(start))
        self.assertTrue(self.probe.seq_str.endswith(end))

    def make_probes_with_adapters(self, probe_str_a, probe_str_b):
        """Make probes with both 'A' and 'B' adapters.

        Args:
            probe_str_a: list of strings of the sequences of probes that
                should receive an 'A' adapter
            probe_str_b: list of strings of the sequences of probes that
                should receive a 'B' adapter

        Returns:
            list of probes (instances of probe.Probe) from the input
            with the corresponding adapters added
        """
        probes = []
        for p_str in probe_str_a:
            probes += [probe.Probe.from_str(p_str).with_prepended_str(
                af.ADAPTER_A_5END).with_appended_str(af.ADAPTER_A_3END)]
        for p_str in probe_str_b:
            probes += [probe.Probe.from_str(p_str).with_prepended_str(
                af.ADAPTER_B_5END).with_appended_str(af.ADAPTER_B_3END)]
        return probes

    def convert_target_genomes(self, target_genomes):
        """Convert genomes to instances of genome.Genome.

        Args:
            target_genomes: nested list of genomes, as strings, to be
                converted

        Returns:
            nested list of genomes, with the same structure as the input,
            in which each genome is an instance of genome.Genome instead
            of a string
        """
        r = []
        for genomes_from_group in target_genomes:
            rg = []
            for g in genomes_from_group:
                rg += [genome.Genome.from_one_seq(g)]
            r += [rg]
        return r

    def test_one_genome(self):
        target_genomes = [['ABCDEFGHIJKLMNOPQRSTUVWXYZ']]
        target_genomes = self.convert_target_genomes(target_genomes)
        # Create probes of length 6 bp with a stride of 3 bp
        input = []
        for genomes_from_group in target_genomes:
            for g in genomes_from_group:
                input += cp.make_candidate_probes_from_sequences(
                    g.seqs,
                    probe_length=6,
                    probe_stride=3)

        f, output = self.get_filter_and_output(6, 0, target_genomes, input, 3,
                                               10)
        desired_output = self.make_probes_with_adapters(['ABCDEF', 'GHIJKL',
                                                         'MNOPQR', 'STUVWX'],
                                                        ['DEFGHI', 'JKLMNO',
                                                         'PQRSTU', 'UVWXYZ'])
        self.assertItemsEqual(output, desired_output)

    def test_two_genomes(self):
        target_genomes = [['ABCDEFGHIJKLMNOPQRSTUVWXYZ'],
                          ['ZYXWVUTSRQPONMLKJIHGFEDCBA']]
        target_genomes = self.convert_target_genomes(target_genomes)
        # Create probes of length 6 bp with a stride of 3 bp
        input = []
        for genomes_from_group in target_genomes:
            for g in genomes_from_group:
                input += cp.make_candidate_probes_from_sequences(
                    g.seqs,
                    probe_length=6,
                    probe_stride=3)

        f, output = self.get_filter_and_output(6, 0, target_genomes, input, 3,
                                               10)
        desired_output = self.make_probes_with_adapters(
            ['ABCDEF', 'GHIJKL', 'MNOPQR', 'STUVWX', 'ZYXWVU', 'TSRQPO',
             'NMLKJI', 'HGFEDC'], ['DEFGHI', 'JKLMNO', 'PQRSTU', 'UVWXYZ',
                                   'WVUTSR', 'QPONML', 'KJIHGF', 'FEDCBA'])
        self.assertItemsEqual(output, desired_output)

    def test_almost_identical_probe(self):
        """Test four probes that align like:
           ------    ------
                ------
                ------
        where the bottom two are the same up to one mismatch. The top
        two probes should be assigned adapter 'A' and the bottom two should
        be assigned adapter 'B'.
        """
        target_genomes = [['ABCDEFGHIJKLMNOP', 'ABCDEFGHXJKLMNOP']]
        target_genomes = self.convert_target_genomes(target_genomes)
        input = ['ABCDEF', 'FGHIJK', 'FGHXJK', 'KLMNOP']
        input = [probe.Probe.from_str(s) for s in input]

        for allowed_mismatches in [0, 1]:
            f, output = self.get_filter_and_output(6, allowed_mismatches,
                                                   target_genomes, input, 3,
                                                   100)
            desired_output = self.make_probes_with_adapters(['ABCDEF',
                                                             'KLMNOP'],
                                                            ['FGHIJK',
                                                             'FGHXJK'])
            self.assertItemsEqual(output, desired_output)

            # Check votes too
            votes = f._make_votes_across_target_genomes(input)
            if allowed_mismatches == 0:
                # Each middle probe should align to one genome
                self.assertEqual(votes, [(2, 0), (0, 1), (0, 1), (2, 0)])
            if allowed_mismatches == 1:
                # Both middle probes should align to both genomes
                self.assertEqual(votes, [(2, 0), (0, 2), (0, 2), (2, 0)])

    def test_misaligned(self):
        """Test probes that align to two genomes, but in which the ones
        aligning to one genome are offset from the other.
        """
        target_genomes = [['ABCDEFGHIJKLMNOPQR', 'XYZABCDEFGHIJKLMNOPQR']]
        target_genomes = self.convert_target_genomes(target_genomes)
        input = ['XYZABC', 'ABCDEF', 'DEFGHI', 'GHIJKL', 'JKLMNO', 'MNOPQR']
        input = [probe.Probe.from_str(s) for s in input]

        f, output = self.get_filter_and_output(6, 0, target_genomes, input, 3,
                                               10)

        # Assume 'ABCDEF' gets 'A' adapter and 'XYZABC' gets 'B' adapter,
        # and so on. But flipping the 'A' and 'B' adapters would also
        # be OK.

        desired_output = self.make_probes_with_adapters(['ABCDEF', 'GHIJKL',
                                                         'MNOPQR'],
                                                        ['XYZABC', 'DEFGHI',
                                                         'JKLMNO'])
        self.assertItemsEqual(output, desired_output)

        # Check votes too
        votes = f._make_votes_across_target_genomes(input)
        self.assertEqual(votes, [(0, 1), (2, 0), (0, 2), (2, 0), (0, 2),
                                 (2, 0)])

    def test_three_genomes(self):
        """Test probes that align adjacent to each other in one genome,
        but overlapping in two others. One should be assigned adapter 'A'
        and the other should be assigned adapter 'B'.
        """
        target_genomes = [['ABCDEFGHEFKLMN', 'ABCDEFKLMN', 'ABCDEFKLMNO']]
        target_genomes = self.convert_target_genomes(target_genomes)
        input = ['ABCDEF', 'EFKLMN']
        input = [probe.Probe.from_str(s) for s in input]

        f, output = self.get_filter_and_output(6, 0, target_genomes, input, 3,
                                               10)

        desired_output = self.make_probes_with_adapters(['ABCDEF'], ['EFKLMN'])
        self.assertItemsEqual(output, desired_output)

        # Check votes too
        votes = f._make_votes_across_target_genomes(input)
        self.assertEqual(votes, [(3, 0), (1, 2)])

    def test_with_mismatches(self):
        target_genomes = [['ABCDEFGHIJKLMNO', 'ABCXEFGXIJKXMNO',
                           'ABCDEFGYYJKLMNO', 'ABCDEXGHIJKLXNO',
                           'ABCDEFGHIJKLMNX', 'AXCDEFGHIJKLMNO',
                           'ABCDEFGHIYYLMNO']]
        target_genomes = self.convert_target_genomes(target_genomes)
        input = ['ABCDEF', 'DEFGHI', 'GHIJKL', 'JKLMNO', 'DEFGYY', 'GYYJKL',
                 'IYYLMN']
        input = [probe.Probe.from_str(s) for s in input]

        f, output = self.get_filter_and_output(6, 1, target_genomes, input, 3,
                                               10)
        desired_output = self.make_probes_with_adapters(['ABCDEF', 'GHIJKL',
                                                         'GYYJKL', 'IYYLMN'],
                                                        ['DEFGHI', 'JKLMNO',
                                                         'DEFGYY'])
        self.assertItemsEqual(output, desired_output)

    def tearDown(self):
        # Re-enable logging
        logging.disable(logging.NOTSET)
