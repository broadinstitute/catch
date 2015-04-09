"""Chooses candidate probes using a naive approach.

In particular, reduces a set of candidate probes in a naive way: by
iterating through the set and, for each candidate probe P, removing all
others that satisfy a specified similarity criterion with P (i.e., are
deemed redundant to P).

Note that this is just a heuristic with no guarantee to minimize the
number of probes. The output is dependent on the order in which we
iterate through the probes and it is not obvious whether this even
yields a reasonable approximation in an average case. This was
implemented primarily to replicate previous software and results on
designing probes for hybrid selection.
"""

import logging

from hybseldesign.filter.base_filter import BaseFilter

__author__ = 'Hayden Metsky <hayden@mit.edu>'

logger = logging.getLogger(__name__)


class NaiveRedundantFilter(BaseFilter):

    """Filter that selects candidate probes with a naive approach.
    """

    def __init__(self, are_redundant_fn=None):
        """
        Args:
            are_redundant_fn: function that takes as input two probes
                and returns True iff the two are deemed redundant
        """
        if are_redundant_fn is None:
            # Use the shift and mismatch count method by default, with
            # parameters ensuring that two probes are deemed redundant
            # iff they are identical (i.e., no shift with no mismatches
            # between the two)
            are_redundant_fn = redundant_shift_and_mismatch_count(
                shift=0, mismatch_thres=0)
        self.are_redundant_fn = are_redundant_fn

    def _filter(self, input):
        """Return a subset of the input probes.
        """
        # For each probe P, removes all subsequent probes that are redundant
        # to P, where redundancy is determined by self.are_redundant_fn.
        # It is necessary to keep track of probes to delete by their
        # index (in probe_indices_to_delete) rather than by storing the
        # probe itself (e.g., in probes_to_delete). The reason is that
        # if there are two probes that are identical they will have the
        # same hash and be considered equal by __eq__; if only one is
        # intended to be deleted (i.e., the latter one in the list of input),
        # they will both be deleted accidentally.
        probe_indices_to_delete = set()
        for i in xrange(len(input)):
            if i % 100 == 0:
                logger.info("Processing candidate probe %d of %d",
                            i, len(input))

            if i in probe_indices_to_delete:
                continue
            probe_a = input[i]
            for j in xrange(i + 1, len(input)):
                if j in probe_indices_to_delete:
                    continue
                probe_b = input[j]
                if self.are_redundant_fn(probe_a, probe_b):
                    probe_indices_to_delete.add(j)

        # Return all probes except those whose indices are in
        # probe_indices_to_delete
        return [p for i, p in enumerate(input)
                if i not in probe_indices_to_delete]


def redundant_shift_and_mismatch_count(
        shift=0,
        mismatch_thres=0,
        quick=True,
        quick_mismatch_cutoff=10):
    """Return a function for determining probe redundancy.

    The returned function determines whether two probes are redundant
    based on the minimum number of mismatches as one probe is shifted
    relative to the other.

    Args:
        shift: one probe is shifted between -shift and +shift relative
            to the other
        mismatch_thres: if the smallest number of mismatches encountered
            during shifting is <= this value, then the returned function
            returns True; otherwise it returns False
        quick: when this is True and mismatch_thres is low (less than
            quick_mismatch_cutoff), the returned function determines
            redundancy in a way that directly looks at the mismatch
            threshold in an attempt to reduce the number of comparisons
        quick_mismatch_cutoff: see quick

    Returns:
        function that returns True or False depending on whether two
        probes are redundant
    """
    # The 'quick' are_redundant function will become slower than
    # the shorter one below as mismatch_thres grows larger, so
    # arbitrarily set a cutoff (quick_mismatch_cutoff) and only
    # return the 'quick' function when mismatch_thres is below
    # this cutoff.
    if quick and mismatch_thres < quick_mismatch_cutoff:
        def are_redundant(probe_a, probe_b):
            probe_a_len = len(probe_a.seq)
            probe_b_len = len(probe_b.seq)
            for s in xrange(-shift, shift + 1):
                mismatches = 0
                if s < 0:
                    probe_a_idx = 0
                    probe_b_idx = -s
                else:
                    probe_a_idx = s
                    probe_b_idx = 0
                while probe_a_idx < probe_a_len and \
                        probe_b_idx < probe_b_len:
                    # Step through the probes, and stop comparing for this
                    # shift if there are too many mismatches
                    if probe_a.seq[probe_a_idx] != probe_b.seq[probe_b_idx]:
                        mismatches += 1
                    if mismatches > mismatch_thres:
                        break
                    probe_a_idx += 1
                    probe_b_idx += 1
                if mismatches <= mismatch_thres:
                    # Found a shift with a small enough number of mismatches
                    return True
            # All of the shifts have > mismatch_thres mismatches
            return False
    else:
        def are_redundant(probe_a, probe_b):
            mismatches = probe_a.min_mismatches_within_shift(probe_b, shift)
            return mismatches <= mismatch_thres
    return are_redundant


def redundant_longest_common_substring(
        mismatches=0,
        lcf_thres=100,
        prune_with_heuristic=True):
    """Return a function for determining probe redundancy.

    The returned function determines whether two probes are redundant
    based on the length of the longest common substring with less than
    or equal to 'mismatches' mismatches.

    Args:
        mismatches/lcf_thres: if the length of the longest common
            substring with at most 'mismatches' mismatches is >=
            'lcf_thres', then the returned function returns True;
            otherwise it returns False
        prune_with_heuristic: when True, a heuristic
            (probe.shares_some_kmers) is called to determine whether the
            two probes might be redundant. If it outputs that they are
            likely not, then the returned function outputs False and the
            longest common substring (a costly operation) is not computed.
            Both false positives and false negatives are possible, but
            should be rare.

    Returns:
        function that returns True or False depending on whether two
        probes are redundant
    """
    def are_redundant(probe_a, probe_b):
        if prune_with_heuristic:
            if not probe_a.shares_some_kmers(probe_b):
                # probe_a and probe_b are likely not redundant, so don't
                # bother computing their longest common substring
                return False
        lcf_length = probe_a.longest_common_substring_length(probe_b,
                                                             mismatches)
        return lcf_length >= lcf_thres
    return are_redundant
