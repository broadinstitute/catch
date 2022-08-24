"""Tests for set_cover_filter module.
"""

from collections import OrderedDict
import logging
import os
import tempfile
import unittest

from catch.filter import set_cover_filter as scf
from catch import genome
from catch import probe
from catch.utils import interval

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class TestSetCoverFilter(unittest.TestCase):
    """Tests the set cover filter output on contrived input.
    """

    def setUp(self):
        # Disable logging
        logging.disable(logging.WARNING)

    def get_filter_and_output(self, lcf_thres, mismatches,
                              target_genomes_grouped,
                              input, coverage,
                              mismatches_tolerant=-1,
                              lcf_thres_tolerant=-1,
                              cover_extension=0,
                              identify=False,
                              avoided_genomes=[],
                              force_num_processes=None):
        input_probes_grouped = []
        for input_group in input:
            i = [probe.Probe.from_str(s) for s in input_group]
            # Remove duplicates
            input_probes_grouped += [list(OrderedDict.fromkeys(i))]
        f = scf.SetCoverFilter(
            mismatches=mismatches,
            lcf_thres=lcf_thres,
            coverage=coverage,
            cover_extension=cover_extension,
            mismatches_tolerant=mismatches_tolerant,
            lcf_thres_tolerant=lcf_thres_tolerant,
            identify=identify,
            avoided_genomes=avoided_genomes,
            kmer_probe_map_k=3)
        if force_num_processes is not None:
            f._force_num_processes = force_num_processes
        output_probes = f.filter(input_probes_grouped, target_genomes_grouped,
                input_is_grouped=True)
        # Flatten output_probes across groupings
        output_probes = list(set(p for i in output_probes for p in i))
        return (f, output_probes)

    def verify_target_genome_coverage(self, selected_probes, target_genomes,
                                      filter, desired_coverage,
                                      cover_extension=0):
        kmer_probe_map = probe.SharedKmerProbeMap.construct(
            probe.construct_kmer_probe_map_to_find_probe_covers(
                selected_probes, filter.mismatches, filter.lcf_thres,
                min_k=3, k=3)
        )
        probe.open_probe_finding_pool(kmer_probe_map,
                                      filter.cover_range_fn)
        for tg in [tg for genomes_from_group in target_genomes for tg in
                genomes_from_group]:
            num_bp_covered = 0
            for seq in tg.seqs:
                probe_cover_ranges = probe.find_probe_covers_in_sequence(seq)
                all_cover_ranges = []
                for cover_ranges in probe_cover_ranges.values():
                    for cv in cover_ranges:
                        start = max(0, cv[0] - cover_extension)
                        end = min(len(seq), cv[1] + cover_extension)
                        all_cover_ranges += [(start, end)]
                all_cover_ranges = interval.merge_overlapping(all_cover_ranges)
                for cover_range in all_cover_ranges:
                    num_bp_covered += cover_range[1] - cover_range[0]
            if desired_coverage <= 1.0:
                # check fraction covered
                desired_bp_covered = desired_coverage * tg.size()
                self.assertGreaterEqual(num_bp_covered, desired_bp_covered)
            else:
                # directly check num bp covered
                desired_coverage_adjusted = min(desired_coverage, tg.size())
                self.assertGreaterEqual(num_bp_covered,
                                        desired_coverage_adjusted)
        probe.close_probe_finding_pool()

    def run_full_coverage_check_for_target_genomes(self,
            target_genomes_grouped, check_must_have=True,
            force_num_processes=None):
        input = []
        for genomes_from_group in target_genomes_grouped:
            input_group = []
            for tg in genomes_from_group:
                for seq in tg.seqs:
                    input_group += [seq[i:(i + 6)] for i in range(len(seq) - 6 + 1)]
            input += [input_group]
        f, output = self.get_filter_and_output(
            6, 0, target_genomes_grouped, input, 1.0,
            force_num_processes=force_num_processes)
        # output must have probes in must_have_output
        if check_must_have:
            must_have_output = ['OPQRST', 'UVWXYZ', 'FEDCBA', 'ABCDEF', 'ZYXWVF']
            for o in must_have_output:
                self.assertTrue(probe.Probe.from_str(o) in output)
        # verify that each of the target genomes is fully covered
        self.verify_target_genome_coverage(output, target_genomes_grouped, f, 1.0)

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

    def test_full_coverage_one_group(self):
        target_genomes = [['ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEF',
                           'ZYXWVFGHIJWUTSOPQRSTFEDCBAZYXWVF']]
        target_genomes = self.convert_target_genomes(target_genomes)
        for np in [None, 1, 2, 4]:
            # Test with one grouping of two sequences
            self.run_full_coverage_check_for_target_genomes(target_genomes,
                    force_num_processes=np)

    def get_6bp_probes(self, target_genomes_grouped,
                       cover=1.0,
                       cover_extension=0,
                       identify=False,
                       mismatches_tolerant=0,
                       lcf_thres_tolerant=6,
                       avoided_genomes=[]):
        input = []
        for genomes_from_group in target_genomes_grouped:
            input_group = []
            for tg in genomes_from_group:
                for seq in tg.seqs:
                    input_group += [seq[i:(i + 6)] for i in range(len(seq) - 6 + 1)]
            input += [input_group]
        f, output = self.get_filter_and_output(
            6, 0, target_genomes_grouped, input, cover,
            mismatches_tolerant=mismatches_tolerant,
            lcf_thres_tolerant=lcf_thres_tolerant,
            cover_extension=cover_extension,
            identify=identify,
            avoided_genomes=avoided_genomes)
        # Test that the output is the same even when particular numbers of
        # processes are set
        for np in [1, 2, 4]:
            _, o = self.get_filter_and_output(
                6, 0, target_genomes_grouped, input, cover,
                mismatches_tolerant=mismatches_tolerant,
                lcf_thres_tolerant=lcf_thres_tolerant,
                cover_extension=cover_extension,
                identify=identify,
                avoided_genomes=avoided_genomes,
                force_num_processes=np)
            self.assertEqual(output, o)
        return f, output

    def test_same_output_with_duplicated_species(self):
        """Tests that, when the same species shows twice (with the same
        target genomes), the output is the same as if it shows once. In
        this unit test we are not removing duplicate probes before
        testing the set cover filter; therefore, we should remove duplicates
        before comparing (since the output when the species shows twice
        will have more probes than when it shows once).
        """
        target_genomes = [['ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEF',
                           'ZYXWVFGHIJWUTSOPQRSTFEDCBAZYXWVF']]
        target_genomes = self.convert_target_genomes(target_genomes)
        _, probes_once = self.get_6bp_probes(target_genomes)
        target_genomes = [['ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEF',
                           'ZYXWVFGHIJWUTSOPQRSTFEDCBAZYXWVF'],
                          ['ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEF',
                           'ZYXWVFGHIJWUTSOPQRSTFEDCBAZYXWVF']]
        target_genomes = self.convert_target_genomes(target_genomes)
        _, probes_twice = self.get_6bp_probes(target_genomes)
        # Compare set to ignore duplicates
        self.assertEqual(set(probes_once), set(probes_twice))

    def test_fractional_coverage(self):
        target_genomes = [['ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEF',
                           'ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEF']]
        target_genomes = self.convert_target_genomes(target_genomes)
        for cover_frac in [0.1, 0.5, 0.8, 1.0]:
            f, probes = self.get_6bp_probes(target_genomes, cover_frac)
            # Note that with cover_frac==0.5, it can be done with just 2
            # probes ('ABCDEF' and one other)
            min_num_probes = {0.1: 1, 0.5: 2, 0.8: 4, 1.0: 5}
            self.assertEqual(len(probes), min_num_probes[cover_frac])
            self.verify_target_genome_coverage(probes, target_genomes,
                                               f, cover_frac)

    def test_explicit_bp_coverage(self):
        target_genomes = [['ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEF',
                           'ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEF']]
        target_genomes = self.convert_target_genomes(target_genomes)
        for num_bp in [2, 5, 10, 15, 20, 1000]:
            f, probes = self.get_6bp_probes(target_genomes, num_bp)
            # Note that with num_bp==10, it can be done with just 1 probe
            # ('ABCDEF')
            min_num_probes = {2: 1, 5: 1, 10: 1, 15: 2, 20: 3, 1000: 5}
            self.assertEqual(len(probes), min_num_probes[num_bp])
            self.verify_target_genome_coverage(probes, target_genomes,
                                               f, num_bp)

    def test_varying_probe_length(self):
        target_genomes = [['ABCDEFGHIJKLM',
                           'ABCXE',
                           'CXEGH']]
        target_genomes = self.convert_target_genomes(target_genomes)
        candidate_probes = [['ABCDEF', 'DEFGHI', 'GHIJKLM', 'ABCXE', 'CXEGH']]
        f, probes = self.get_filter_and_output(
            5, 0, target_genomes, candidate_probes, 1.0)
        probes_str = [''.join(x.seq) for x in probes]
        self.assertCountEqual(probes_str,
                              ['ABCDEF', 'GHIJKLM', 'ABCXE', 'CXEGH'])
        self.verify_target_genome_coverage(probes, target_genomes, f, 1.0)

    def test_cover_extension1(self):
        target_genomes = [['ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEF',
                           'ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEF']]
        target_genomes = self.convert_target_genomes(target_genomes)
        f, probes = self.get_6bp_probes(target_genomes, cover_extension=2)
        # It is possible to achieve full coverage (counting the extensions)
        # with 3 probes
        self.assertEqual(len(probes), 3)
        self.verify_target_genome_coverage(probes, target_genomes,
                                           f, 1.0,
                                           cover_extension=2)

    def test_cover_extension2(self):
        target_genomes = [['ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEF',
                           'ZYXWVFGHIJWUTSOPQRSTFEDCBAZYXWVF']]
        target_genomes = self.convert_target_genomes(target_genomes)
        f, probes = self.get_6bp_probes(target_genomes, cover_extension=3)
        # It is possible to achieve full coverage (counting the extensions)
        # with 5 probes
        self.assertEqual(len(probes), 5)
        self.verify_target_genome_coverage(probes, target_genomes,
                                           f, 1.0,
                                           cover_extension=3)

    def test_cover_extension3(self):
        target_genomes = [['ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEF',
                           'ZYXWVFGHIJWUTSOPQRSTFED']]
        target_genomes = self.convert_target_genomes(target_genomes)
        f, probes = self.get_6bp_probes(target_genomes, cover_extension=3)
        # It is possible to achieve full coverage (counting the extensions)
        # with 4 probes
        self.assertEqual(len(probes), 4)
        self.verify_target_genome_coverage(probes, target_genomes,
                                           f, 1.0,
                                           cover_extension=3)

    def test_cover_extension_with_partial_coverage(self):
        target_genomes = [['ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEF',
                           'ZYXWVFGHIJWUTSOPQRSTFEDCBAZYXWVF']]
        target_genomes = self.convert_target_genomes(target_genomes)
        f, probes = self.get_6bp_probes(target_genomes,
                                        cover=0.5, cover_extension=3)
        # It is possible to achieve 50% coverage (counting the extensions)
        # with 3 probes
        self.assertEqual(len(probes), 3)
        self.verify_target_genome_coverage(probes, target_genomes,
                                           f, 0.5,
                                           cover_extension=3)

    def test_identify_one_group(self):
        target_genomes = [['ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEF',
                           'ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEF']]
        target_genomes = self.convert_target_genomes(target_genomes)
        f, probes = self.get_6bp_probes(target_genomes, cover=6, identify=True)
        self.assertEqual(probes, [probe.Probe.from_str('ABCDEF')])

    def test_identify_two_groups(self):
        target_genomes = [['ABCDEFXXIJKXMNOPQRXTUXWXYXABCDEF',
                           'ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEF'],
                          ['ATATATABCDEFATATATATATATATATATAT']]
        target_genomes = self.convert_target_genomes(target_genomes)
        f, probes = self.get_6bp_probes(target_genomes, cover=6, identify=True)
        self.assertEqual(set(probes), {probe.Probe.from_str('MNOPQR'),
                                           probe.Probe.from_str('ATATAT')})

    def test_identify_three_groups(self):
        target_genomes = [['ABCDEFQRSQRSHIJKLMQRSQRSQRSQRSQR',
                           'XYZXYZATATATAXYZXYZXYZEEEEEEXYZX'],
                          ['ATATATABXDXFATATATACGCGCGTATATAT',
                           'CGCGCGABCDEFATXTATATATATATATATAT'],
                          ['XYZXYZAAAAAAXYZXYZXYZXYZXYZXYZXY',
                           'QRSQRSQRSQRAAAAAAQRSQRSQRSQRSQRS']]
        target_genomes = self.convert_target_genomes(target_genomes)
        f, probes = self.get_6bp_probes(target_genomes, cover=6, identify=True)
        # CGCGCG for second group
        self.assertIn(probe.Probe.from_str('CGCGCG'), probes)
        # AAAAAA for third group
        self.assertIn(probe.Probe.from_str('AAAAAA'), probes)
        # 2 probes needed for first group (one per genome)
        self.assertEqual(len(probes), 4)

    def test_identify_three_groups_forced_pick(self):
        target_genomes = [['ABCDEFXYZXYZIJKLMN', 'XYZXYZBCDEFMNOPQ'],
                          ['ABCDEFMNOPQR'], ['ABCDEF']]
        target_genomes = self.convert_target_genomes(target_genomes)
        f, probes = self.get_6bp_probes(target_genomes, cover=6, identify=True)
        # ABCDEF is forced to be chosen due to the third group,
        # but XYZXYZ and MNOPQR will first be chosen for the first and
        # second groups
        self.assertEqual(set(probes), {probe.Probe.from_str('ABCDEF'),
                                           probe.Probe.from_str('XYZXYZ'),
                                           probe.Probe.from_str('MNOPQR')})

    def test_identify_three_groups_two_hit_species(self):
        target_genomes = [['ABCDEFXYZXYZ', 'MNOPQRXYZXYZ'], ['ABCDEFXYZXYZ'],
                          ['ABCDEFMNOPQR']]
        target_genomes = self.convert_target_genomes(target_genomes)
        f, probes = self.get_6bp_probes(target_genomes, cover=6, identify=True)
        # ABCDEF should not be chosen because it hits all three groups,
        # but some probe(s) that are selected will have to hit two groups
        self.assertNotIn(probe.Probe.from_str('ABCDEF'), probes)
        # MNOPQR should not be chosen because it is possible to cover
        # the second genome of the first group with NOPQRX (or some other
        # probes) and it is possible to cover the genome in the last
        # group with BCDEFM
        self.assertNotIn(probe.Probe.from_str('MNOPQR'), probes)

    def test_identify_two_groups_two_probes(self):
        target_genomes = [['ABCDEFXXIJKXMNOPQRXTUVWXYXABCDEF',
                           'TUVWXYGHIJKLMNOPQRSABCDEFAABCDEF'],
                          ['ATATATABCDEFATATATATATATATATATAT']]
        target_genomes = self.convert_target_genomes(target_genomes)
        f, probes = self.get_6bp_probes(target_genomes,
                                        cover=10,
                                        identify=True)
        self.assertEqual(set(probes), {probe.Probe.from_str('MNOPQR'),
                                           probe.Probe.from_str('TUVWXY'),
                                           probe.Probe.from_str('ATATAT')})

    def test_identify_two_groups_tolerant(self):
        target_genomes = [['ABCDEFXXIJKXMNOPQRXTATXAYABCDEFATAXATXYZX',
                           'ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEF'],
                          ['ATATATABCDEFATATATATATATATXYZXYZ']]
        target_genomes = self.convert_target_genomes(target_genomes)
        f, probes = self.get_6bp_probes(target_genomes,
                                        cover=6,
                                        mismatches_tolerant=1,
                                        lcf_thres_tolerant=5,
                                        identify=True)
        self.assertEqual(set(probes), {probe.Probe.from_str('MNOPQR'),
                                           probe.Probe.from_str('XYZXYZ')})

    def test_identify_two_groups_reverse_complement(self):
        target_genomes = [['ATCGGGXXIJKXMNOPQRXTUXWXYXATCGGG',
                           'ATCGGGGHIJKLMNOPQRSTUVWXYZATCGGG'],
                          ['ATATATCCCGATATATATATATATATATATAT']]
        target_genomes = self.convert_target_genomes(target_genomes)
        f, probes = self.get_6bp_probes(target_genomes, cover=6, identify=True)
        self.assertEqual(set(probes), {probe.Probe.from_str('MNOPQR'),
                                           probe.Probe.from_str('ATATAT')})

    def test_avoid_one_genome1(self):
        bl_file = tempfile.NamedTemporaryFile(mode='w')
        bl_file.write(">n/a\n")
        bl_file.write("AAAAAAAAAAAAAAAAAAAAA\n")
        bl_file.seek(0)

        target_genomes = [['ABCDEFXXIJKXMNOPQRXTUXWXYXABCDEF',
                           'ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEF']]
        target_genomes = self.convert_target_genomes(target_genomes)
        f, probes = self.get_6bp_probes(target_genomes,
                                        cover=6,
                                        identify=False,
                                        avoided_genomes=[bl_file.name])
        # No candidate probe is avoided
        self.assertEqual(set(probes), {probe.Probe.from_str('ABCDEF')})

        bl_file.close()

    def test_avoid_one_genome2(self):
        bl_file = tempfile.NamedTemporaryFile(mode='w')
        bl_file.write(">n/a\n")
        bl_file.write("AAAAAAAAATCGGGAAAAAAAA\n")
        bl_file.seek(0)

        target_genomes = [['ATCGGGXXIJKXMNOPQRXTUXWXYXATCGGG',
                           'ATCGGGGHIJKLMNOPQRSTUVWXYZATCGGG']]
        target_genomes = self.convert_target_genomes(target_genomes)
        f, probes = self.get_6bp_probes(target_genomes,
                                        cover=6,
                                        identify=False,
                                        avoided_genomes=[bl_file.name])
        # ABCDEF is avoided, so go for the second most common probe
        self.assertEqual(set(probes), {probe.Probe.from_str('MNOPQR')})

        bl_file.close()

    def test_avoid_one_genome_reverse_complement(self):
        bl_file = tempfile.NamedTemporaryFile(mode='w')
        bl_file.write(">n/a\n")
        bl_file.write("AAAAAAAACCCGATAAAAAA\n")
        bl_file.seek(0)

        target_genomes = [['ATCGGGXXIJKXMNOPQRXTUXWXYXATCGGG',
                           'ATCGGGGHIJKLMNOPQRSTUVWXYZAYCGGG']]
        target_genomes = self.convert_target_genomes(target_genomes)
        f, probes = self.get_6bp_probes(target_genomes,
                                        cover=6,
                                        identify=False,
                                        avoided_genomes=[bl_file.name])
        # ATCGGG is avoided, so go for the second most common probe
        self.assertEqual(set(probes), {probe.Probe.from_str('MNOPQR')})

        bl_file.close()

    def test_avoid_one_genome_tolerant(self):
        bl_file = tempfile.NamedTemporaryFile(mode='w')
        bl_file.write(">n/a\n")
        bl_file.write("AAAAAAAATCCGCAAAAAAAA\n")
        bl_file.seek(0)

        target_genomes = [['ATCGGGXXIJKXMNOPQRXTUXWXYXATCGGG',
                           'ATCGGGGHIJKLMNOPQRSTUVWXYZAYCGGG']]
        target_genomes = self.convert_target_genomes(target_genomes)
        f, probes = self.get_6bp_probes(target_genomes,
                                        cover=6,
                                        identify=False,
                                        mismatches_tolerant=1,
                                        lcf_thres_tolerant=5,
                                        avoided_genomes=[bl_file.name])
        # ATCGGG is avoided, so go for the second most common probe
        self.assertEqual(set(probes), {probe.Probe.from_str('MNOPQR')})

        bl_file.close()

    def test_avoid_two_genomes_one_file(self):
        bl_file = tempfile.NamedTemporaryFile(mode='w')
        bl_file.write(">n/a 1\n")
        bl_file.write("AAAAAAAAATCGGGAAAAAAAA\n")
        bl_file.write(">n/a 2\n")
        bl_file.write("AATCGGGAAAAAAAAGGGGGGAAAA\n")
        bl_file.seek(0)

        target_genomes = [['ATCGGGXXIJKXGGGGGGXTUXWXYXATCGGG',
                           'ATCGGGGHIJKLGGGGGGSTUVWXYZATCGGG']]
        target_genomes = self.convert_target_genomes(target_genomes)
        f, probes = self.get_6bp_probes(target_genomes,
                                        cover=6,
                                        identify=False,
                                        avoided_genomes=[bl_file.name])
        self.assertNotIn(probe.Probe.from_str('ATCGGG'), probes)
        self.assertNotIn(probe.Probe.from_str('GGGGGG'), probes)

        bl_file.close()

    def test_avoid_two_genomes_two_files(self):
        bl_file1 = tempfile.NamedTemporaryFile(mode='w')
        bl_file1.write(">n/a 1\n")
        bl_file1.write("AAAAAAAAATCGGGAAAAAAAA\n")
        bl_file1.seek(0)
        bl_file2 = tempfile.NamedTemporaryFile(mode='w')
        bl_file2.write(">n/a 1\n")
        bl_file2.write("AATCGGGAAAAAAAAGGGGGGAAAA\n")
        bl_file2.seek(0)

        target_genomes = [['ATCGGGXXIJKXGGGGGGXTUXWXYXATCGGG',
                           'ATCGGGGHIJKLGGGGGGSTUVWXYZATCGGG']]
        target_genomes = self.convert_target_genomes(target_genomes)
        f, probes = self.get_6bp_probes(
            target_genomes,
            cover=6,
            identify=False,
            avoided_genomes=[bl_file1.name, bl_file2.name])
        self.assertNotIn(probe.Probe.from_str('ATCGGG'), probes)
        self.assertNotIn(probe.Probe.from_str('GGGGGG'), probes)

        bl_file1.close()
        bl_file2.close()

    def test_avoid_one_genome_forced_pick(self):
        bl_file = tempfile.NamedTemporaryFile(mode='w')
        bl_file.write(">n/a\n")
        bl_file.write("AAAAAAAAAAATCGGGAAAAA\n")
        bl_file.seek(0)

        target_genomes = [['ABCDEFABCDEF'], ['ABCDEFATCGGG']]
        target_genomes = self.convert_target_genomes(target_genomes)
        f, probes = self.get_6bp_probes(target_genomes,
                                        cover=1.0,
                                        identify=False,
                                        avoided_genomes=[bl_file.name])
        # Should choose ABCDEF
        self.assertIn(probe.Probe.from_str('ABCDEF'), probes)
        # Forced to choose ATCGGG at end to ensure full coverage
        self.assertIn(probe.Probe.from_str('ATCGGG'), probes)
        # Before choosing ATCGGG, it should have made one other try
        # (e.g., FATCGG)
        self.assertEqual(len(probes), 3)

        bl_file.close()

    def test_identify_and_avoid(self):
        bl_file = tempfile.NamedTemporaryFile(mode='w')
        bl_file.write(">n/a\n")
        bl_file.write("AAAAAAAAAAATCGGGATCGGGAAAAA\n")
        bl_file.seek(0)

        target_genomes = [['ABCDEFGGGGGGCCCCCC'],
                          ['ABCDEFATCGGGATCGGGXXX',
                           'ATCGGGBCDEFGGGGGCCCCCATCGGGYYY']]
        target_genomes = self.convert_target_genomes(target_genomes)
        f, probes = self.get_6bp_probes(target_genomes,
                                        cover=12,
                                        identify=True,
                                        avoided_genomes=[bl_file.name])
        # Should pick GGGGGG and CCCCCC for the first genome
        # Note there are just 5 G's and 5 C's in the last genome
        self.assertIn(probe.Probe.from_str('GGGGGG'), probes)
        self.assertIn(probe.Probe.from_str('CCCCCC'), probes)
        # Should avoid ABCDEF because it hits two groups
        self.assertNotIn(probe.Probe.from_str('ABCDEF'), probes)
        # Should avoid ATCGGG because it's avoided
        self.assertNotIn(probe.Probe.from_str('ATCGGG'), probes)
        self.verify_target_genome_coverage(probes, target_genomes, f, 12)

        bl_file.close()

    def test_full_coverage_two_groups(self):
        target_genomes = [['ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEF'],
                          ['ZYXWVFGHIJWUTSOPQRSTFEDCBAZYXWVF']]
        target_genomes = self.convert_target_genomes(target_genomes)
        for np in [None, 1, 2, 4]:
            # Do not check for the 'must have' probes because the targets will
            #   be covered separately, and thus probes may not reflect the
            #   most conserved sites
            self.run_full_coverage_check_for_target_genomes(target_genomes,
                    check_must_have=False, force_num_processes=np)

    def test_custom_cover_range_fn(self):
        custom_path = os.path.join(os.path.dirname(__file__),
                                   'input/custom_cover_range_fn.py')
        custom_cover_range_fn = (custom_path, 'covers_abc')

        target_genomes = [['AAAAAAAAABCBBBBBBBBBB',
                           'AAAAAAAAABCBBBBBBBBBB']]
        target_genomes = self.convert_target_genomes(target_genomes)
        candidate_probes= [[probe.Probe.from_str(p) for p in
                            ['AAAAAA', 'AAABCB', 'BBBBBB', 'XXXXXX']]]

        # Only try to cover 3 bp; the only probe that will hybridize
        # by this toy cover_range_fn is the one with 'ABC'
        f = scf.SetCoverFilter(0, 0,
                               coverage=3,
                               custom_cover_range_fn=custom_cover_range_fn,
                               kmer_probe_map_k=3)
        output_probes = f.filter(candidate_probes, target_genomes,
                input_is_grouped=True)
        output_probes = list(set(p for i in output_probes for p in i)) # flatten
        self.assertEqual(set(output_probes), {probe.Probe.from_str('AAABCB')})

    def tearDown(self):
        # Re-enable logging
        logging.disable(logging.NOTSET)
