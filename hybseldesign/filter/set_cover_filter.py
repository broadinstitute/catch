"""Chooses candidate probes using a set cover approach.

In particular, reduces a set of candidate probes by treating the problem
as an instance of the set cover problem.

Each candidate probe is treated as a set whose elements consist of the
bases of the target genomes that the probe 'covers' (i.e., should
hybridize to). This uses a version of the problem in which there are
multiple universes -- each universe corresponds to one target genome,
and the solution obtained is guaranteed to cover at least some
specified fraction of each universe. That is, the selected probes
collectively cover a specified fraction of the bases in each target
genome.

There is an identification option that, when enabled, instructs the
filter to choose candidate probes that will make it easy -- after
hybrid selection and sequencing are performed -- to identify one
or more grouping(s) of the target genomes from a sample. (The input
target genomes are grouped. Each grouping could represent, for
example, a species; the candidate probes would then make it easy to
identify one or more species present in a sample.) Ideally, the
chosen candidate probes each "hit" just one grouping of the target
genomes, where a probe "hits" a grouping if it covers at least one
target genome in that grouping. Candidate probes that hit more than
one grouping (these would make identification more difficult) are
only chosen if absolutely necessary to achieve the desired coverage.
When identification is enabled, the specified coverage to achieve
should typically be small.

There is a 'blacklist' of genomes, whose sequences we do not want to
cover (e.g., human DNA when we are interested in extracting viral
sequences). Probes that cover a portion of these genomes are
penalized in finding the set cover solution. That is, this filter
tries to choose candidate probes that do NOT cover a portion of
the blacklisted genomes, and will only choose candidate probes that
DO cover a portion of these genomes if doing so is absolutely
necessary to achieve the desired coverage.

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

logger = logging.getLogger(__name__)


class SetCoverFilter(BaseFilter):

  """Filter that selects candidate probes using a set cover approach.
  """

  def __init__(self, mismatches=0, lcf_thres=100,
      mismatches_tolerant=None, lcf_thres_tolerant=None,
      identify=False, blacklisted_genomes=[], coverage=1.0,
      kmer_size=15, num_kmers_per_probe=10):
    """
    Args:
        mismatches/lcf_thres: consider a probe to hybridize to a sequence
            if a stretch of 'lcf_thres' or more bp aligns with
            'mismatches' or fewer mismatched bp; used to compute whether
            a probe "covers" a portion of a sequence
        mismatches_tolerant/lcf_thres_tolerant: more tolerant
            values corresponding to 'mismatches' and 'lcf_thres'. It should
            generally be true that 'mismatches_tolerant' > 'mismatches' and
            'lcf_thres_tolerant' < 'lcf_thres'. These values are used in
            determining the overlap that a candidate probe has with
            different groupings when the identification option is enabled.
            They are also used when determining the coverage of each
            candidate probe with the blacklisted genomes. They are meant
            to capture more potential hybridizations (i.e., be more
            sensitive). When not set, they are by default equal to
            'mismatches' and 'lcf_thres'.
        identity: when True, indicates that probes should be designed
            with the identification option enabled (default is False)
        blacklisted_genomes: list of paths to FASTA files of genomes
            that should be blacklisted (i.e., probes are penalized by the
            amount they cover these genomes).
        coverage: either a float in [0,1] or an int > 1. When it is a
            float in [0,1], it determines the fraction of each of the
            target genomes that must be covered by the selected probes.
            When it is an int > 1, it determines the number of bp of each
            of the target genomes that must be covered by the selected
            probes.
    """
    self.cover_range_fn = \
        probe.probe_covers_sequence_by_longest_common_substring(
            mismatches, lcf_thres)

    if not mismatches_tolerant:
      mismatches_tolerant = mismatches
    if not lcf_thres_tolerant:
      lcf_thres_tolerant = lcf_thres
    self.cover_range_tolerant_fn = \
        probe.probe_covers_sequence_by_longest_common_substring(
            mismatches_tolerant, lcf_thres_tolerant)

    # Warn if identification is enabled but the coverage is high
    if identify:
      if (coverage <= 1.0 and coverage >= 0.25) or \
         (coverage > 1 and coverage >= 5000):
        logger.warning(("Identification is enabled but the required "
                        "coverage is high; generally coverage should "
                        "be small when performing identification"))

    self.identify = identify
    self.blacklisted_genomes = blacklisted_genomes
    self.coverage = coverage
    self.kmer_size = kmer_size
    self.num_kmers_per_probe = num_kmers_per_probe

  def _make_sets(self, candidate_probes, kmer_probe_map):
    """Return a collection of sets to use in set cover.

    In the returned collection of sets, each set corresponds to a
    candidate probe and contains the bases of the target genomes
    covered by the candidate probe. The target genomes must be in
    grouped lists inside the list self.target_genomes.

    The output is intended for input to set_cover.approx_multiuniverse
    as the 'sets' input.

    Args:
        candidate_probes: list of candidate probes
        kmer_probe_map: dict mapping kmers to probes

    Returns:
        a dict mapping set_ids (from 0 through
        len(candidate_probes)-1) to dicts, where the dict for a
        particular set_id maps universe_ids to sets. set_id
        corresponds to a candidate probe in candidate_probes and
        universe_id is a tuple that corresponds to a target genome in
        a grouping from self.target_genomes. The j'th target genome
        from the i'th grouping in self.target_genomes is given
        universe_id equal to (i,j). That is, i ranges from 0 through
        len(self.target_genomes)-1 (i.e., the number of groupings) and
        j ranges from 0 through (n_i)-1 where n_i is the number of
        target genomes in the i'th group. In the returned value
        (sets), sets[set_id][universe_id] is a set of all the bases
        (as integers) covered by probe set_id in the target genome
        universe_id.
    """
    probe_id = {}
    sets = {}
    for id, p in enumerate(candidate_probes):
      probe_id[p] = id
      # Store values in an array of type 'I' to be space efficient
      sets[id] = defaultdict(lambda: array('I'))

    for i, genomes_from_group in enumerate(self.target_genomes):
      for j, gnm in enumerate(genomes_from_group):
        logger.info(("Computing coverage in grouping %d (of %d), "
                     "with target genome %d (of %d)"),
                    i, len(self.target_genomes),
                    j, len(genomes_from_group))
        universe_id = (i,j)
        length_so_far = 0
        for sequence in gnm.seqs:
          probe_cover_ranges = probe.find_probe_covers_in_sequence(
              sequence, kmer_probe_map, k=self.kmer_size,
              cover_range_for_probe_in_subsequence_fn=self.cover_range_fn)
          # Add the bases of sequence that are covered by all the
          # probes into sets with universe_id equal to (i,j)
          for p, cover_ranges in probe_cover_ranges.iteritems():
            set_id = probe_id[p]
            for cover_range in cover_ranges:
              for bp in xrange(cover_range[0], cover_range[1]):
                # bp gives a position in just this sequence
                # (chromosome), so adding the lengths of all the
                # sequences previously iterated (length_so_far) onto
                # bp gives a unique integer position in the genome gnm
                sets[set_id][universe_id].append(length_so_far + bp)
          length_so_far += len(sequence)

    # Convert each defaultdict to dict
    for set_id in sets.keys():
      sets[set_id] = dict(sets[set_id])

    return sets

  def _compute_tolerant_bp_covered_within_sequence(self,
      candidate_probes, kmer_probe_map, sequence, rc_too=True):
    """Compute number of bp captured in sequence by each input probe.

    Uses self.coverage_range_tolerant_fn for determining coverage (i.e.,
    the coverage is determined in a relatively tolerant way so that
    more potential hybridizations are included).

    Args:
        candidate_probes: list of candidate probes for which to
            determine coverage
        kmer_probe_map: dict mapping kmers to probes
        sequence: sequence as a string in which to determine the
            coverage of the probes
        rc_too: when True, the returned values also include bp that
            are captured in the reverse complement of sequence

    Returns:
        dict mapping each candidate probe to the number of bp it
        covers, for only the candidate probes that cover at least
        one bp; candidate probes that do not cover any bp are not
        included as keys in the returned dict
    """
    reverse_complement = [False]
    if rc_too:
      reverse_complement += [True]
    rc_map = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

    num_bp_covered = defaultdict(int)

    for rc in reverse_complement:
      if rc:
        sequence = ''.join([rc_map.get(b,b) for b in sequence[::-1]])
      probe_cover_ranges = probe.find_probe_covers_in_sequence(
          sequence, kmer_probe_map, k=self.kmer_size,
          cover_range_for_probe_in_subsequence_fn=\
            self.cover_range_tolerant_fn)

      all_cover_ranges = []
      for p, cover_ranges in probe_cover_ranges.iteritems():
        for cover_range in cover_ranges:
          num_bp_covered[p] += cover_range[1] - cover_range[0]

    return dict(num_bp_covered)

  def _count_num_groupings_hit(self, candidate_probes,
      kmer_probe_map):
    """Compute number of genome groupings hit by each candidate probe.

    A probe is said to "hit" a grouping of target genomes if it covers
    at least one bp of at least one target genome in the grouping.
    This decides whether a probe covers part of a target genome in
    a tolerant way (i.e., using self.cover_range_tolerant_fn) so that
    more potential hits are counted.

    Args:
        candidate_probes: list of candidate probes
        kmer_probe_map: dict mapping kmers to probes

    Returns:
        dict mapping each candidate probe to the number of target
        genome groupings it hits
    """
    num_groupings_hit = { p: 0 for p in candidate_probes }
    for i, genomes_from_group in enumerate(self.target_genomes):
      logger.info(("Computing coverage in grouping %d (of %d) to "
                   "count number of groupings hit"),
                  i, len(self.target_genomes))
      num_bp_covered_in_grouping = defaultdict(int)
      for j, gnm in enumerate(genomes_from_group):
        for sequence in gnm.seqs:
          # Count hits in both sequence and its reverse complement
          num_bp = self._compute_tolerant_bp_covered_within_sequence(
            candidate_probes, kmer_probe_map, sequence, rc_too=True)
          for p in num_bp.keys():
            num_bp_covered_in_grouping[p] += num_bp[p]
      # If a probe covers at least one bp in this grouping (i),
      # then it hits this grouping
      for p in num_bp_covered_in_grouping.keys():
        if num_bp_covered_in_grouping[p] >= 1:
          num_groupings_hit[p] += 1

    # Check that each candidate probe hits at least one grouping
    for p, hit in num_groupings_hit.iteritems():
      if hit == 0:
        # Something's strange! Every probe should hit one or more
        # target genome groupings because each candidate probe
        # is created from a target genome.
        logger.critical(("There is a probe that does not 'hit' "
                         "any target genome grouping, but every "
                         "candidate probe should hit at least one"))

    return num_groupings_hit

  def _count_blacklisted_bp_covered(self, candidate_probes,
      kmer_probe_map):
    """Compute number of blacklisted genome bp covered by each probe.

    This decides whether a candidate probe captures a portion of a
    blacklisted genome in a tolerant way (i.e., using
    self.cover_range_tolerant_fn) so that more potential hybridizations
    are counted. Also, the total number of bp includes bases of
    the reverse complement of each blacklisted genome that are covered
    by a probe, so that both a blacklisted genome and its reverse
    complement are blacklisted.

    Args:
        candidate_probes: list of candidate probes
        kmer_probe_map: dict mapping kmers to probes

    Returns:
        dict mapping each candidate probe to the total number of bp
        it covers in the blacklisted genomes and their reverse
        complements
    """
    total_num_bp = { p: 0 for p in candidate_probes }
    for fasta_path in self.blacklisted_genomes:
      # Use a generator to read the FASTA to avoid loading too much
      # into memory (e.g., only store one chromosome of the human
      # genome at a time)
      for sequence in seq_io.iterate_fasta(fasta_path):
        logger.info(("Computing coverage across a blacklisted "
                     "sequence"))
        # Blacklist both sequence and its reverse complement
        num_bp = self._compute_tolerant_bp_covered_within_sequence(
          candidate_probes, kmer_probe_map, sequence, rc_too=True)
        for p in num_bp.keys():
          total_num_bp[p] += num_bp[p]
    return total_num_bp

  def _make_ranks(self, candidate_probes, kmer_probe_map):
    """Return a rank for each candidate probe to use in set cover.

    The "rank" of a candidate probe is a level of penalty for that
    probe, where higher ranks are more penalized. A set cover is sought
    that uses as many candidate probes from rank i as possible before
    considering probes with rank i+1. There are two considerations in
    computing ranks:
      - When identification is turned on (i.e., self.identify is True),
        the number of species that a probe "hits". Fewer hit species
        yields a smaller rank.
      - The number of bases in blacklisted genomes that the probe
        covers. Fewer covered bases yields a smaller rank.
    A probe that covers any part of a blacklisted genome will always
    receive a higher rank than a probe that does not. (This is achieved
    by first computing ranks using tuples of the form (x,y) where x=0
    for any probe that does not cover a blacklisted genome and x=1
    for a probe that does; y determines relative rank among those probes
    with the same x value. The tuple ranks are then converted into
    integer ranks by sorting the tuples.) When identification is
    enabled, a probe that hits more than one grouping (e.g., species)
    will always receive a higher rank than a probe that only hits one
    grouping (and does not cover any blacklisted genomes). 

    When identification is not turned on, weighted set cover
    effectively does the following: (1) Covers as much of the target
    genomes as possible while minimizing the number of probes, without
    using any probe that covers any part of a blacklisted genome. (2)
    Covers whatever portions of the target genomes remain to be covered
    by using probes that cover parts of blacklisted genomes, while
    first seeking probes that cover less of the blacklisted genomes
    (i.e., even if probe B covers much more of the target genomes
    than probe A, A will be chosen before B if B covers a tiny bit
    more of the blacklisted genomes than A).
    When identification is turned on, weighted set cover: (1) Covers
    as much of the target genomes as possible while minimizing the
    number of probes, only using probes that hit one grouping. (2)
    Covers whatever portions of the target genomes remain to be covered
    while minimizing the number of probes, only using probes that hit
    two groupings, etc. (3) Considers probes that cover parts of
    blacklisted genomes, if there remains more of the target genomes to
    cover.

    The output is intended for input to set_cover.approx_multiuniverse
    as the 'ranks' input.

    Args:
        candidate_probes: list of candidate probes
        kmer_probe_map: dict mapping kmers to probes

    Returns:
        dict mapping set_ids (0 through len(candidate_probes)-1, each
        corresponding to a candidate probe) to a rank (integer) for
        that candidate probe
    """
    if self.identify:
      # Find the number of target genome groupings (e.g., species)
      # that each probe "hits". (A probe "hits" a grouping if it
      # covers a part of at least one target genome in that grouping.)
      # A probe that hits just one grouping is good for
      # identification and is therefore ranked relatively low (a
      # rank of 1); probes that hit more than one grouping are poor
      # for identification and their ranks are equal to the number
      # of groupings they hit.
      num_groupings_hit = self._count_num_groupings_hit(
          candidate_probes, kmer_probe_map)
      rank_val = { p: (0, hit) for p, hit \
                    in num_groupings_hit.iteritems() }
    else:
      # Start each probe with the same rank
      rank_val = { p: (0, 0) for p in candidate_probes }

    # Find probes that cover part of a blacklisted genome.
    # All of these get a higher rank than any probe that does not
    # cover any part of a blacklisted genome (since the first element
    # of the tuple put into rank_val is 1, but 0 was the first
    # element of the tuple above) and the rank among these is based
    # on the number of bp they cover.
    blacklisted_bp_covered = self._count_blacklisted_bp_covered(
        candidate_probes, kmer_probe_map)
    for p, bp in blacklisted_bp_covered.iteritems():
      if bp > 0:
        rank_val[p] = (1, bp)

    # Convert the ranks, specified as tuples, into ranks from 0
    # upward. The probe(s) with the smallest tuple rank get(s)
    # rank 0, the probe(s) with the next smallest tuple rank get(s)
    # rank 1, and so on..
    all_rank_tuples = sorted(set(rank_val.values()))
    tuple_rank_idx = {}
    for i in xrange(len(all_rank_tuples)):
      tuple_rank_idx[all_rank_tuples[i]] = i
    ranks = {}
    for set_id, p in enumerate(candidate_probes):
      ranks[set_id] = tuple_rank_idx[rank_val[p]]

    return ranks

  def _make_costs(self, candidate_probes):
    """Return a cost for each candidate probe to use in set cover.

    The ranks computed by self._make_ranks(..) effectively serve as
    costs for this application. (That is, a rank can be thought of as a
    sufficiently large cost.) This returns a cost of 1 for each
    candidate probe, so that the only factor distinguishing any two
    probes with the same rank is the number of bp of the target genomes
    they cover -- i.e., costs do not play a factor among probes with
    the same rank.

    The output is intended for input to set_cover.approx_multiuniverse
    as the 'costs' input.

    Args:
        candidate_probes: list of candidate probes

    Returns:
      dict mapping set_ids (0 through len(candidate_probes)-1, each
      corresponding to a candidate probe) to a cost (integer) for that
      candidate probe; currently the cost is always 1
    """
    return { set_id: 1 for set_id in xrange(len(candidate_probes)) }

  def _make_universe_p(self):
    """Return a required coverage for each universe to use in set cover.

    The output is intended for input to set_cover.approx_multiunvierse
    as the 'universe_p' input.

    Returns:
        dict mapping each universe_id (representing a target genome) to
        the fraction of the target genome that must be covered, as
        determined from self.coverage
    """
    universe_p = {}
    if self.coverage <= 1.0:
      # self.coverage should explicitly be the fraction of each
      # target genome to cover
      logger.info(("Building universe_p directly from desired "
                   "fractional coverage"))
      for i in xrange(len(self.target_genomes)):
        for j in xrange(len(self.target_genomes[i])):
          universe_p[(i,j)] = self.coverage
    else:
      # self.coverage should be an int representing the number of
      # bp of each target genome to cover; convert it into a
      # fraction using the size of each target genome
      logger.info(("Building universe_p from desired number of bp "
                   "to cover"))
      for i in xrange(len(self.target_genomes)):
        for j, gnm in enumerate(self.target_genomes[i]):
          universe_p[(i,j)] = float(self.coverage) / gnm.size()
    return universe_p

  def _filter(self, input):
    """Return a subset of the input probes.
    """
    # Ensure that the input is a list
    input = list(input)

    logger.info("Building map from k-mers to probes")
    kmer_probe_map = probe.construct_kmer_probe_map(input,
        k=self.kmer_size,
        num_kmers_per_probe=self.num_kmers_per_probe,
        include_positions=True)

    logger.info("Building set cover sets input")
    sets = self._make_sets(input, kmer_probe_map)
    logger.info("Building set cover ranks input")
    ranks = self._make_ranks(input, kmer_probe_map)
    logger.info("Building set cover costs input")
    costs = self._make_costs(input)
    logger.info("Building set cover universe_p input")
    universe_p = self._make_universe_p()

    # Run the set cover approximation algorithm
    logger.info(("Approximating the solution to the set cover "
                 "instance"))
    set_ids_in_cover = set_cover.approx_multiuniverse(
                        sets, costs=costs, universe_p=universe_p,
                        ranks=ranks, use_arrays=True)

    # Save ranks and costs
    self.probe_ranks = ranks
    self.probe_costs = costs

    # Warn when less-than-ideal probes are chosen (i.e., probes
    # whose ranks exceed 0)
    num_bad_probes = sum([True for set_id in set_ids_in_cover \
                            if ranks[set_id] > 0])
    if num_bad_probes > 0:
      logger.warning(("Forced to choose %d less-than-ideal probe%s "
                      "(i.e., probes that 'hit' more than one "
                      "grouping during identification or probes that "
                      "cover a blacklisted genome)"), num_bad_probes,
                      ('' if num_bad_probes == 1 else 's'))

    return [input[id] for id in set_ids_in_cover]

