"""Tests for set_cover_filter module.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

import unittest
import logging

from hybseldesign import genome
from hybseldesign import probe
from hybseldesign.filter import set_cover_filter as scf
from hybseldesign.utils import interval


"""Tests the set cover filter output on contrived input.
"""
class TestSetCoverFilter(unittest.TestCase):

  def setUp(self):
    # Disable logging
    logging.disable(logging.INFO)

  def get_filter_and_output(self, lcf_thres, mismatches,
      target_genomes, input, coverage, k,
      num_kmers_per_probe):
    input_probes = [probe.Probe.from_str(s) for s in input]
    f = scf.SetCoverFilter(mismatches=mismatches, lcf_thres=lcf_thres,
          coverage=coverage, kmer_size=k,
          num_kmers_per_probe=num_kmers_per_probe)
    f.target_genomes = target_genomes
    f.filter(input_probes)
    return (f, f.output_probes)

  def verify_target_genome_coverage(self, selected_probes,
      target_genomes, k, num_kmers_per_probe, filter,
      desired_coverage):
    kmer_probe_map = probe.construct_kmer_probe_map(
                      selected_probes, k=k,
                      num_kmers_per_probe=10,
                      include_positions=True)
    for tg in [g for genomes_from_group in target_genomes \
          for g in genomes_from_group]:
      num_bp_covered = 0
      for seq in tg.seqs:
        probe_cover_ranges = probe.find_probe_covers_in_sequence(
            seq, kmer_probe_map, k=k,
            cover_range_for_probe_in_subsequence_fn=filter.cover_range_fn)
        all_cover_ranges = []
        for cover_ranges in probe_cover_ranges.values():
          for cv in cover_ranges:
            all_cover_ranges += [cv]
        all_cover_ranges = interval.merge_overlapping(all_cover_ranges)
        for cover_range in all_cover_ranges:
          num_bp_covered += cover_range[1] - cover_range[0]
      if desired_coverage <= 1.0:
        # check fraction covered
        desired_bp_covered = desired_coverage * tg.size()
        self.assertGreaterEqual(num_bp_covered, desired_bp_covered)
      else:
        # directly check num bp covered
        self.assertGreaterEqual(num_bp_covered, desired_coverage)

  def run_full_coverage_check_for_target_genomes(self, target_genomes):
    input = []
    for tg in [g for genomes_from_group in target_genomes \
          for g in genomes_from_group]:
      for seq in tg.seqs:
        input += [ seq[i:(i+6)] for i in xrange(len(seq)-6+1) ]
    must_have_output = [ 'OPQRST',
                         'UVWXYZ',
                         'FEDCBA',
                         'ABCDEF',
                         'ZYXWVF' ]
    f, output = self.get_filter_and_output(6, 0, target_genomes,
        input, 1.0, 3, 10)
    # output must have probes in must_have_output
    for o in must_have_output:
      self.assertTrue(probe.Probe.from_str(o) in output)
    # verify that each of the target genomes is fully covered
    self.verify_target_genome_coverage(output,
      target_genomes, 3, 10, f, 1.0)

  """Returns a nested list of target_genomes in which each
  genome is an instance of genome.Genome, given a nested list in
  which the genomes are strings.
  """
  def convert_target_genomes(self, target_genomes):
    r = []
    for genomes_from_group in target_genomes:
      rg = []
      for g in genomes_from_group:
        rg += [ genome.Genome.from_one_seq(g) ]
      r += [ rg ]
    return r

  def test_full_coverage(self):
    target_genomes = [ [ 'ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEF',
                         'ZYXWVFGHIJWUTSOPQRSTFEDCBAZYXWVF' ] ]
    target_genomes = self.convert_target_genomes(target_genomes)
    # Test with one grouping of two sequences
    self.run_full_coverage_check_for_target_genomes(target_genomes)
    # Test again with two groupings of one sequence each
    target_genomes = [ [ 'ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEF' ],
                       [ 'ZYXWVFGHIJWUTSOPQRSTFEDCBAZYXWVF' ] ]
    target_genomes = self.convert_target_genomes(target_genomes)
    self.run_full_coverage_check_for_target_genomes(target_genomes)

  """Tests that, when the same species shows twice (with the same
  target genomes), the output is the same as if it shows once. In
  this unit test we are not removing duplicate probes before
  testing the set cover filter; therefore, we should remove duplicates
  before comparing (since the output when the species shows twice
  will have more probes than when it shows once).
  """
  def test_same_output_with_duplicated_species(self):
    def get_probes(target_genomes):
      input = []
      for tg in [g for genomes_from_group in target_genomes \
            for g in genomes_from_group]:
        for seq in tg.seqs:
          input += [ seq[i:(i+6)] for i in xrange(len(seq)-6+1) ]
      f, output = self.get_filter_and_output(6, 0, target_genomes,
          input, 1.0, 3, 10)
      return output
    target_genomes = [ [ 'ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEF',
                         'ZYXWVFGHIJWUTSOPQRSTFEDCBAZYXWVF' ] ]
    target_genomes = self.convert_target_genomes(target_genomes)
    probes_once = get_probes(target_genomes)
    target_genomes = [ [ 'ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEF',
                         'ZYXWVFGHIJWUTSOPQRSTFEDCBAZYXWVF' ],
                       [ 'ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEF',
                         'ZYXWVFGHIJWUTSOPQRSTFEDCBAZYXWVF' ] ]
    target_genomes = self.convert_target_genomes(target_genomes)
    probes_twice = get_probes(target_genomes)
    # Compare set to ignore duplicates
    self.assertEqual(set(probes_once), set(probes_twice))

  def test_fractional_coverage(self):
    def get_probes(target_genomes, cover_frac):
      input = []
      for tg in [g for genomes_from_group in target_genomes \
            for g in genomes_from_group]:
        for seq in tg.seqs:
          input += [ seq[i:(i+6)] for i in xrange(len(seq)-6+1) ]
      f, output = self.get_filter_and_output(6, 0, target_genomes,
          input, cover_frac, 3, 10)
      return f, output
    target_genomes = [ [ 'ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEF',
                         'ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEF' ] ]
    target_genomes = self.convert_target_genomes(target_genomes)
    for cover_frac in [0.1, 0.5, 0.8, 1.0]:
      f, probes = get_probes(target_genomes, cover_frac)
      # note that with cover_frac==0.5, it can be done with just 2
      # probes ('ABCDEF' and one other)
      min_num_probes = { 0.1: 1, 0.5: 2, 0.8: 4, 1.0: 5 }
      self.assertEqual(len(probes), min_num_probes[cover_frac])
      self.verify_target_genome_coverage(probes,
        target_genomes, 3, 10, f, cover_frac)

  def test_explicit_bp_coverage(self):
    def get_probes(target_genomes, num_bp):
      input = []
      for tg in [g for genomes_from_group in target_genomes \
            for g in genomes_from_group]:
        for seq in tg.seqs:
          input += [ seq[i:(i+6)] for i in xrange(len(seq)-6+1) ]
      f, output = self.get_filter_and_output(6, 0, target_genomes,
          input, num_bp, 3, 10)
      return f, output
    target_genomes = [ [ 'ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEF',
                         'ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEF' ] ]
    target_genomes = self.convert_target_genomes(target_genomes)
    for num_bp in [2, 5, 10, 15, 20]:
      f, probes = get_probes(target_genomes, num_bp)
      # note that with num_bp==10, it can be done with just 1 probe
      # ('ABCDEF')
      min_num_probes = { 2: 1, 5: 1, 10: 1, 15: 2, 20: 3 }
      self.assertEqual(len(probes), min_num_probes[num_bp])
      self.verify_target_genome_coverage(probes,
        target_genomes, 3, 10, f, num_bp)

  def tearDown(self):
    # Re-enable logging
    logging.disable(logging.NOTSET)

