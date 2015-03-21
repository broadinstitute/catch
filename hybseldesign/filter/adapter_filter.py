"""Filters probes by returning each input probe with an adapter
(short DNA sequence) prepended to its 5' end and an adapter appended
to its 3' end.

We need to place adapters (primers) on the end of each probe for use
in PCR. However, probes that "overlap" (roughly speaking) can cause
problems (explained later) in PCR if they are PCR'd together. Thus,
we wish to use k different adapters such that a set of probes that
do not "overlap" each other get a unique adapter and this set is
PCR'd together. Typically, in this use case k=2 is sufficient --
each adapter is either an 'A' adapter or a 'B' adapter -- and each
probe is given either an 'A' or 'B' adapter (on both ends); the
'A' adapter probes are PCR'd separately from the 'B' adapter probes.
If we could define "overlap" more precisely, we could construct a
graph in which each vertex corresponds to a probe and edges give
the probes that "overlap". We would then wish to color the vertices
with k colors so as to minimize the number of monochromatic edges
(i.e., minimize the number of edges whose two vertices have the
same color). The colors correspond to adapters, and this would
effectively provide us with sets of probes for each adapter that
minimize the overlap among probes that share an adapter. This is
the min k-partition problem (dual of the max k-cut problem).

However, the problem comes only from a more specific meaning of
"overlap". In particular, suppose a probe is 100 bp and an adapter
(on each end) is 20 bp, so that in total the synthesized probe is
140 bp. Problems arise after PCR when probes we desire have annealed
so as to deviate from 140 bp. Namely, this happens when probes
overlap in a way that allows them to "chain" together (in a sort of
ladder structure) that allows them to grow long. For example,
consider the following probes staggered across a target genome:
                   ----- ----- ----- -----
                      ----- ----- -----
These overlap in a way that (considering that complements arise
during PCR) would allow them to chain together. Thus, we would want
all probes on the top row to be given adapter 'A' and all on the
bottom row to be given adapter 'B'. This way, they will be PCR'd
separately. Now consider the following probes along a target genome:
                      -----   -----
                          -----
                          -----
This could arise if the two probes on the two bottom rows differ
in just one position (e.g., together they cover a variant) but
may still hybridize to the same sequence. While the probes could
hybridize to each other (specifically, the complement of one to
the other), this is not a problem. Rather, the problem is that, like
the example above, they could form a chain. So we would want to
assign the two probes on the top row each adapter 'A', and the two
bottom probes each adapter 'B'.

The problem itself is difficult: we want to find subsets of the
probes such that the probes in each subset overlap in a way that
could allow them to form a chain, and then assign each of the probes
in these subsets a different adapter. But there is a heuristic that
should be suitable here: we can treat the problem as an interval
scheduling problem.

For each target genome we scan to find where all of the probes would
hybridize. Based on the "intervals" where these hybridize, we select
the maximum number of non-overlapping intervals that span the target
genome (i.e., the maximum number of non-overlapping probes). We
find these by solving an instance of the interval scheduling problem
with a simple greedy algorithm. We assign all of these probes adapter
'A', and then assign all of the remaining probes adapter 'B'. (If we
have more than two adapters, and thus more than two PCR experiments,
we could solve the interval scheduling problem again on all of the
remaining probes, assign the selected ones adapter 'B' and the
further remaining ones adapter 'C', and so on..) This has some
downsides. For example, we are implicitly assuming that each probe
hybridizes just once in a target genome, while in reality it may
appear more than once. (We assign it adapter 'A' if it is ever
selected by the interval scheduling algoritm, even though it may
overlap with another adapter 'A' probe.) However, in general for
this use case it should be rare for a probe to appear a large number
of times in a single target genome. Another downside is that --
assuming each probe appears just once in a target genome -- there
will be no overlap (or potential for chain formation) among probes
with adapter 'A', whereas there may be overlap among probes with
adapter 'B'. However, while possible, it should be rare for probes
with adapter 'B' to overlap because, in this use case, the coverage
at any bp should generally stay below 2x.

We solve interval scheduling for each target genome. Each time, we
vote once (either 'A' or 'B') for each probe. (Each target genome
gives just one vote to each probe, or none if that probe does not
hybridize to the target genome.) But the distinction between 'A' and
'B' is arbitrary -- at a particular target genome, the two can be
exchanged with each other for all the probes. (That is, at a target
genome we could vote for 'A' instead of 'B' and vice-versa for each
probe.) We seek to maximize the sum, across all probes, of the
majority vote for the probe (to ensure a clear choice of adapter for
each probe). For a probe p, let V^A_p be the number of 'A' adapter
votes p receives, and let V^B_p be the number of 'B' votes. Then, we
seek to maximize, summed over all probes p, max(V^A_p, V^B_p). When
we add target genome n+1, we can either exchange all 'A' and 'B'
votes from that target genome or not. We exchange them if and only if
doing so yields a higher total sum than not doing so. We only need
to consider exchanging the votes for target genome n+1; assuming this
has already been done earlier for each of the n preceding target
genomes, when they were added, there is no need consider exchanging
votes from those. This process yields the maximum sum considering
exchanges of votes in each of the target genomes, and this can be
proven by induction on the number of target genomes.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

import logging

from hybseldesign import probe
from hybseldesign.filter.base_filter import BaseFilter
from hybseldesign.utils import interval

logger = logging.getLogger(__name__)

ADAPTER_A_5END = 'ATACGCCATGCTGGGTCTCC'
ADAPTER_A_3END = 'CGTACTTGGGAGTCGGCCAT'
ADAPTER_B_5END = 'AGGCCCTGGCTGCTGATATG'
ADAPTER_B_3END = 'GACCTTTTGGGACAGCGGTG'


class AdapterFilter(BaseFilter):

  def __init__(self, mismatches=0, lcf_thres=100,
      kmer_size=15, num_kmers_per_probe=10):
    self.cover_range_fn = \
        probe.probe_covers_sequence_by_longest_common_substring(
            mismatches=mismatches, lcf_thres=lcf_thres)
    self.kmer_size = kmer_size
    self.num_kmers_per_probe = num_kmers_per_probe

  """Compute votes for probes based on their hybridization within
  the sequence 'sequence' (e.g., one target genome).

  We use the greedy interval scheduling algorithm and assign 'A'
  votes to all probes selected by this algorithm. All other probes
  that hybridize to 'sequence' but are not selected receive a 'B'
  vote.

  Returns a list L, in which L[i] corresponds to the probe probes[i].
  L[i] is either (1,0) [vote for 'A'], (0,1) [vote for 'B'], or (0,0)
  [the probe does not hybridize in 'sequence'].
  """
  def _votes_in_sequence(self, probes, sequence, kmer_probe_map):
    probe_cover_ranges = probe.find_probe_covers_in_sequence(
        sequence, kmer_probe_map, k=self.kmer_size,
        cover_range_for_probe_in_subsequence_fn=self.cover_range_fn)
    aligned_probes = probe_cover_ranges.keys()
    # Make a list of all the intervals covered by all the probes,
    # along with a reference to the probe with the interval
    intervals = []
    for p, cover_ranges in probe_cover_ranges.iteritems():
      for cover_range in cover_ranges:
        intervals += [ (cover_range, p) ]

    # Perform interval scheduling to choose probes that should be
    # assigned the 'A' adapter
    chosen_probes = set(interval.schedule(intervals))

    votes = []
    for p in probes:
      if p in chosen_probes:
        # vote for 'A'
        vote = (1,0)
      else:
        if p in aligned_probes:
          # p should have been skipped by the interval scheduling
          # algorithm
          # vote for 'B'
          vote = (0,1)
        else:
          # p does not hybridize to sequence
          vote = (0,0)
      votes += [vote]
    return votes

  """Exchanges, for each probe in 'votes', an 'A' vote with a 'B'
  vote and vice-versa.
  """
  def _flip_AB_votes(self, votes):
    votes_flipped = []
    for vote in votes:
      assert len(vote) == 2
      votes_flipped += [ (vote[1], vote[0]) ]
    return votes_flipped

  """Sums, over each probe in 'votes', the count of the vote that
  received the most number of votes.
  """
  def _sum_plurality_vote_across_probes(self, votes):
    return sum(max(v) for v in votes)

  """Adds, for each probe, all of the votes in 'votes_a' to the
  votes in 'votes_b', and returns the summed votes.

  For example, let votes_a=[(1,2), (0,3)] and votes_b=[(0,1), (1,1)].
  (These represent two probes.) The summed votes, as returned, are
  [(1,3), (1,4)].
  """
  def _sum_votes_per_probe(self, votes_a, votes_b):
    assert len(votes_a) == len(votes_b)
    votes_sum = []
    for i in xrange(len(votes_a)):
      vote_a, vote_b = votes_a[i], votes_b[i]
      # They must represent the same number of adapters (e.g., for
      # 'A' and 'B' adapters), this length should be 2
      assert len(vote_a) == len(vote_b)
      votes_sum += [ tuple(x+y for x,y in zip(vote_a, vote_b)) ]
    return votes_sum

  """Returns, for each probe in 'probes', votes for the adapter to
  be assigned to the probe.

  Specifically, the returned value is a list L such that L[i] is a
  tuple (A,B) where A gives the number of 'A' adapter votes for
  the probe probes[i] and B gives the number of 'B' adapter votes.
  These are computed, cumulatively, across all target genomes in
  self.target_genomes.
  """
  def _make_votes_across_target_genomes(self, probes):
    logger.info("Building map from k-mers to probes")
    kmer_probe_map = probe.construct_kmer_probe_map(probes,
        k=self.kmer_size,
        num_kmers_per_probe=self.num_kmers_per_probe,
        include_positions=True)

    # Store adapter votes for each probe in a list where the element
    # at index i is a tuple (A,B) that corresponds to the probe
    # probes[i] where A gives the 'A' votes for the probe and B gives
    # the 'B' votes
    cumulative_votes = [(0,0) for _ in xrange(len(probes))]
    for i, sequence in enumerate(self.target_genomes):
      # Compute votes for the adapters for each probe in 'sequence',
      # and also exchange all 'A' votes with 'B' votes and vice-versa.
      # Determine whether or not the exchange matches better with
      # cumulative_votes so far, and update cumulative_votes
      # accordingly.
      logger.info("Computing adapter votes across genome %d of %d",
          i, len(self.target_genomes))
      votes = self._votes_in_sequence(probes, sequence,
                kmer_probe_map)
      votes_flipped = self._flip_AB_votes(votes)
      cumulative_votes_with_nonflipped = self._sum_votes_per_probe(
          cumulative_votes, votes)
      sum_nonflipped = self._sum_plurality_vote_across_probes(
          cumulative_votes_with_nonflipped)
      cumulative_votes_with_flipped = self._sum_votes_per_probe(
          cumulative_votes, votes_flipped)
      sum_flipped = self._sum_plurality_vote_across_probes(
          cumulative_votes_with_flipped)
      if sum_flipped > sum_nonflipped:
        # Add onto cumulative votes the votes in 'votes_flipped'
        # because these could be said to yield a more decisive
        # choice of adapter for each probe (i.e., the sum, across
        # all probes, of the most common vote of adapter for the
        # probe is higher) than the (unflipped) votes in 'votes'
        cumulative_votes = cumulative_votes_with_flipped
      else:
        cumulative_votes = cumulative_votes_with_nonflipped
    return cumulative_votes

  def _filter(self, input):
    # Ensure that the input is a list
    input = list(input)

    logger.info("Computing adapter votes across all target genomes")
    votes = self._make_votes_across_target_genomes(input)

    # Using votes, select the adapter for each probe (the one with
    # the most number of votes) and create a new probe that has
    # this adapter on both ends
    logger.info("Adding adapters to probes based on votes")
    input_with_adapters = []
    for i in xrange(len(input)):
      p = input[i]
      vote = votes[i]
      # Only work with 'A' and 'B' adapters
      assert len(vote) == 2
      if vote[0] > vote[1]:
        # Add an 'A' adapter
        new_p = p.with_prepended_str(ADAPTER_A_5END).\
                  with_appended_str(ADAPTER_A_3END)
      else:
        # Add a 'B' adapter
        new_p = p.with_prepended_str(ADAPTER_B_5END).\
                  with_appended_str(ADAPTER_B_3END)
      input_with_adapters += [new_p]
    return input_with_adapters

