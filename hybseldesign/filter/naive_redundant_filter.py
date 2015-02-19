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
    pass


"""Returns a function that determines whether two probes are
redundant based on the minimum number of mismatches as one probe
is shifted relative to the other.

One probe is shifted between -shift and +shift relative to the
other. If the smallest number of mismatches encountered is
is <= mismatch_thres, then are_redundant outputs True; otherwise
it outputs False.
"""
def are_redundant_based_on_shift_and_mismatch_count(shift=0,
    mismatch_thres=0):
  def are_redundant(probe_a, probe_b):
    mismatches = probe_a.min_mismatches_within_shift(probe_b, shift)
    return mismatches <= mismatch_thres
  return are_redundant


"""Returns a function that determines whether two probes are
redundant based on the length of the longest common substring
with less than or equal to 'mismatches' mismatches.

If the length of the longest common substring with at most
'mismatches' mismatches is >= lcf_thres, then are_redundant outputs
True; otherwise it outputs False.
"""
def are_redundant_based_on_longest_common_substring(mismatches=0,
    lcf_thres=100):
  def are_redundant(probe_a, probe_b):
    lcf_length = probe_a.longest_common_substring_length(probe_b,
                    mismatches)
    return lcf_length >= lcf_thres
  return are_redundant
  
