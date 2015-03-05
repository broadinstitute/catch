"""Selects a subset of candidate probes by treating the problem as an
instance of the set cover problem.

Each candidate probe is treated as a set whose elements consist of the
bases of the target genomes that the probe 'covers' (i.e., should
hybridize to). This uses a version of the problem in which there are
multiple universes -- each universe corresponds to one target genome,
and the solution obtained is guaranteed to cover at least some
specified fraction of each universe. That is, the selected probes
collectively cover a specified fraction of the bases in each target
genome.

There is a 'blacklist' of genomes, whose sequences we do not want to
cover (e.g., human DNA when we are interested in extracting viral
sequences). Probes that cover a portion of these genomes are
penalized in finding the set cover solution, where the penalty for
a probe is based on the amount of the 'blacklist' sequences the
probe covers. When the blacklist is provided, the object is to find
probes that minimize the cover of the blacklist genomes while
achieving sufficient coverage of the target genomes (this has an
implicit goal of finding a minimal number of probes). When a
blacklist is not provided, the objective is simply to minimize the
number of probes needed to achieve sufficient coverage of the target
genomes.

The preprocessing determines whether a (portion of a) probe covers a
sequence, and if so the portion of the sequence covered, by
considering the longest common substring with some number of
mismatches between a sequence and a probe.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

import logging
import re
from collections import defaultdict
from array import array

from hybseldesign import probe
from hybseldesign.filter.base_filter import BaseFilter
from hybseldesign.utils import set_cover
from hybseldesign.utils import seq_io

# Setup logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())
logger.setLevel(logging.INFO)


class SetCoverFilter(BaseFilter):

  """A (portion of a) probe is said to cover a portion of a sequence
  if the longest common substring with at most 'mismatches'
  mismatches between the probe and the sequence is at least
  'lcf_thres' bases long.

  'blacklisted_genomes' is a list of paths to FASTA files of genomes
  that should be blacklisted (i.e., probes are penalized by the
  amount they cover these genomes).

  'coverage_frac' is a float in [0,1] that determines the fraction
  of each of the target genomes that must be covered by the selected
  probes.
  """
  def __init__(self, mismatches=0, lcf_thres=100,
      blacklisted_genomes=[], coverage_frac=1.0,
      kmer_size=15, num_kmers_per_probe=10):
    self.cover_range_fn = \
        probe.probe_covers_sequence_by_longest_common_substring(
            mismatches=mismatches, lcf_thres=lcf_thres)
    self.blacklisted_genomes = blacklisted_genomes
    self.coverage_frac = coverage_frac
    self.kmer_size = kmer_size
    self.num_kmers_per_probe = num_kmers_per_probe

  """Returns a collection of sets, in which each set corresponds to
  a candidate probe and contains the bases of the target genomes
  covered by the candidate probe.

  Specifically, the returned value is a dict mapping set_ids (from 0
  through len(candidate_probes)-1) to dicts, where the dict for a
  particular set_id maps universe_ids (0 through
  len(target_genomes)-1) to sets. set_id corresponds to a candidate
  probe in candidate_probes and universe_id corresponds to a target
  genome in self.target_genomes. The i'th target genome in that
  list is given universe_id equal to i. In the returned value, sets,
  sets[set_id][universe_id] is a set of all the bases (as integers)
  covered by probe set_id in the target genome universe_id.

  The target genomes must be in the list self.target_genomes.

  The output is intended for input to set_cover.approx_multiuniverse
  as the 'sets' input.
  """
  def _make_sets(self, candidate_probes, kmer_probe_map):
    probe_id = {}
    sets = {}
    for id, p in enumerate(candidate_probes):
      probe_id[p] = id
      # Store values in an array of type 'I' to be space efficient
      sets[id] = defaultdict(lambda: array('I'))

    for i, sequence in enumerate(self.target_genomes):
      logger.info("  Computing coverage across genome %d of %d",
          i, len(self.target_genomes))
      probe_cover_ranges = probe.find_probe_covers_in_sequence(
          sequence, kmer_probe_map, k=self.kmer_size,
          cover_range_for_probe_in_subsequence_fn=self.cover_range_fn)
      # Add the bases of sequence that are covered by all the probes
      # into sets with universe_id equal to i
      for p, cover_ranges in probe_cover_ranges.iteritems():
        set_id = probe_id[p]
        for cover_range in cover_ranges:
          for bp in xrange(cover_range[0], cover_range[1]):
            sets[set_id][i].append(bp)

    return sets

  """Returns a collection of costs, in which each cost corresponds
  to a candidate probe.

  Specifically, the returned value is a dict mapping set_ids (from 0
  through len(candidate_probes)-1) to numbers, where each set_id
  corresponds to a candidate probe. The cost is a penalty for the
  probe and is equal to one plus the number of bases across the
  blacklisted genomes that the probe covers. (The one is added so
  that, if a probe does not cover any portion of any blacklisted
  genome, it is not automatically taken. It also ensures that, if
  there are no blacklisted genomes, the objective is the minimize
  the number of probes.)

  The output is intended for input to set_cover.approx_multiuniverse
  as the 'costs' input.
  """
  def _make_costs(self, candidate_probes, kmer_probe_map):
    probe_id = {}
    costs = {}
    for id, p in enumerate(candidate_probes):
      probe_id[p] = id
      costs[id] = 1

    for fasta_path in self.blacklisted_genomes:
      # Use a generator to read the FASTA to avoid loading too much
      # into memory (e.g., only store one chromosome of the human
      # genome at a time)
      for sequence in seq_io.iterate_fasta(fasta_path):
        logger.info(("  Computing coverage across a blacklisted "
                     "sequence"))
        # Blacklist both sequence and its reverse complement
        for reverse_complement in [False, True]:
          if reverse_complement:
            rc_map = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
            sequence = ''.join([rc_map.get(b,b) for b in sequence[::-1]])

          probe_cover_ranges = probe.find_probe_covers_in_sequence(
              sequence, kmer_probe_map, k=self.kmer_size,
              cover_range_for_probe_in_subsequence_fn=self.cover_range_fn)
          # Add the number of bases of sequence that are covered by
          # all the probes into costs
          for p, cover_ranges in probe_cover_ranges.iteritems():
            set_id = probe_id[p]
            for cover_range in cover_ranges:
              num_bp_covered = cover_range[1] - cover_range[0]
              costs[set_id] += num_bp_covered

    return costs

  """Returns a dict mapping each universe_id (representing a target
  genome) to self.coverage_frac, the fraction of the target genome
  that must be covered.

  The output is intended for input to set_cover.approx_multiuniverse
  as the 'universe_p' input.
  """
  def _make_universe_p(self):
    num_universes = len(self.target_genomes)
    return { i: self.coverage_frac for i in xrange(num_universes) }

  def _filter(self, input):
    # Ensure that the input is a list
    input = list(input)

    logger.info("Building map from k-mers to probes")
    kmer_probe_map = probe.construct_kmer_probe_map(input,
        k=self.kmer_size,
        num_kmers_per_probe=self.num_kmers_per_probe,
        include_positions=True)

    logger.info("Building set cover sets input")
    sets = self._make_sets(input, kmer_probe_map)
    logger.info("Building set cover costs input")
    costs = self._make_costs(input, kmer_probe_map)
    logger.info("Building set cover universe_p input")
    universe_p = self._make_universe_p()

    # Run the set cover approximation algorithm
    logger.info(("Approximating the solution to the set cover "
                 "instance"))
    set_ids_in_cover = set_cover.approx_multiuniverse(
                        sets, costs=costs, universe_p=universe_p,
                        use_arrays=True)

    # Save costs if there are blacklisted genomes
    if len(self.blacklisted_genomes) > 0:
      self.probe_costs = costs.values()

    return [input[id] for id in set_ids_in_cover]

