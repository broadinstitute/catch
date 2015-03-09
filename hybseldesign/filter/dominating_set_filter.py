"""Reduces a set of candidate probes by treating the problem as an
instance of the dominating set problem.

Each probe is treated as a vertex in the graph. If two probes are
redundant (i.e., are deemed similar by a specified criterion) then
that redundancy is treated as an edge between the vertices
representing the two probes. In seeking to minimize the number of
probes while covering all probes (i.e., ensuring that every probe
is either chosen or is redundant to a chosen probe), we are solving
an instance of the dominating set problem.

We solve the problem by reducing it to an instance of set cover
and calling the utils.set_cover.approx function to compute an
approximate set cover. Because there is an L-reduction between
set cover and dominating set, the approximation guarantees from set
cover apply to the solution here.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

import logging
from collections import defaultdict

from hybseldesign.filter.base_filter import BaseFilter
from hybseldesign.utils import set_cover

logger = logging.getLogger(__name__)


class DominatingSetFilter(BaseFilter):

  """are_redundant_fn is a function that takes as input two probes
  and returns True iff the two are deemed redundant.
  """
  def __init__(self, are_redundant_fn=None):
    if are_redundant_fn == None:
      # Use the shift and mismatch count method by default, with
      # parameters ensuring that two probes are deemed redundant
      # iff they are identical (i.e., no shift with no mismatches
      # between the two)
      are_redundant_fn = redundant_shift_and_mismatch_count(
                          shift=0, mismatch_thres=0)
    self.are_redundant_fn = are_redundant_fn

  """Treat the problem as an instance of dominating set and reduce
  it to an instance of set cover.

  While it is helpful to think of treating the problem as an instance
  of dominating set, we instead directly build the desired instance
  of set cover. Also, rather than working directly with probes we
  instead construct the set cover with integers representing the
  probes (where each integer is the index of the probe in the input
  list).

  Add an edge between two vertices if they are redundant, where
  redundancy is determined by self.are_redundant_fn.
  """
  def _filter(self, input):
    # Ensure that the input is a list
    input = list(input)

    # In the instance of set cover, we have a set S for each probe P
    # in which S consists of P as well as all other probes redundant
    # to P. Construct these sets.
    sets = defaultdict(set)
    for i in xrange(len(input)):
      if i % 100 == 0:
        logger.info("Making set for candidate probe %d of %d",
          i, len(input))
      probe_a = input[i]
      # Put probe_a into its set
      sets[i].add(probe_a)
      # Find all other probes redundant to probe_a
      for j in xrange(i+1, len(input)):
        probe_b = input[j]
        if self.are_redundant_fn(probe_a, probe_b):
          # Put probe_b into probe_a's set (set[i]), and also put
          # probe_a into probe_b's set (set[j])
          sets[i].add(probe_b)
          sets[j].add(probe_a)

    # Run the set cover approximation algorithm
    set_ids_in_cover = set_cover.approx(sets)

    return [input[id] for id in set_ids_in_cover]

