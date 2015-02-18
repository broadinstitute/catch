"""Reduces a set of candidate probes in a naive way: by iterating
through the set and, for each candidate probe P, removing all others
that satisfy a specified similarity criterion with P (i.e., are
deemed redundant to P).

Note that this is just a heuristic with no guarantee to minimize the
number of probes. The output is dependent on the order in which we
iterate through the probes and it is not obvious whether this even
yields a reasonable approximation in an average case. This was
implemented primarily to replicate previous software and results on
designing probes for hybrid selection.
"""

# Author: Hayden Metsky <hayden@mit.edu>

from hybseldesign.filter.base_filter import BaseFilter


class NaiveRedundantFilter(BaseFilter):

  def _filter(self, input):
