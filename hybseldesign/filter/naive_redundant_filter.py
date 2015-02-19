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

  def _filter(self, input):
    probes_to_delete = set()
    for i, probe_a in enumerate(input):
      if probe_a in probes_to_delete:
        continue
      for j in xrange(i+1, len(input)):
        probe_b = input[j]
        if probe_b in probes_to_delete:
          continue
        if self.are_redundant_fn(probe_a, probe_b):
          probes_to_delete.add(probe_b)

    # Return all probes except those in probes_to_delete
    return [p for p in input if p not in probes_to_delete]


"""Returns a function that determines whether two probes are
redundant based on the minimum number of mismatches as one probe
is shifted relative to the other.

One probe is shifted between -shift and +shift relative to the
other. If the smallest number of mismatches encountered is
is <= mismatch_thres, then are_redundant outputs True; otherwise
it outputs False.
"""
def redundant_shift_and_mismatch_count(shift=0,
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
def redundant_longest_common_substring(mismatches=0,
    lcf_thres=100):
  def are_redundant(probe_a, probe_b):
    lcf_length = probe_a.longest_common_substring_length(probe_b,
                    mismatches)
    return lcf_length >= lcf_thres
  return are_redundant
  
