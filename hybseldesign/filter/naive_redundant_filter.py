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

import logging

from hybseldesign.filter.base_filter import BaseFilter

# Setup logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())
logger.setLevel(logging.INFO)

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

  """For each probe P, removes all subsequent probes that are redundant
  to P, where redundancy is determined by self.are_redundant_fn.

  It is necessary to necessary to keep track of probes to delete
  by their index (in probe_indices_to_delete) rather than by
  storing the probe itself (e.g., in probes_to_delete). The reason
  is that if there are two probes that are identical they will have
  the same hash and be considered equal by __eq__; if only one is
  intended to be deleted (i.e., the latter one in the list of input),
  they will both be deleted accidentally.
  """
  def _filter(self, input):
    probe_indices_to_delete = set()
    for i in xrange(len(input)):
      if i % 100 == 0:
        logger.info("Processing candidate probe %d of %d",
          i, len(input))

      if i in probe_indices_to_delete:
        continue
      probe_a = input[i]
      for j in xrange(i+1, len(input)):
        if j in probe_indices_to_delete:
          continue
        probe_b = input[j]
        if self.are_redundant_fn(probe_a, probe_b):
          probe_indices_to_delete.add(j)

    # Return all probes except those whose indices are in
    # probe_indices_to_delete
    return [p for i, p in enumerate(input) \
              if i not in probe_indices_to_delete]


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
  
