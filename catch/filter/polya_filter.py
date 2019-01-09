"""Removes probes with long stretches of A or T bases.

This acts as a filter on the probes by returning a list of
probes, among the input probes, that do _not_ have a long
stretch of A or T, tolerating some given number of
mismatches. It preserves the order of the input.
"""

from collections import OrderedDict

from catch.filter.base_filter import BaseFilter
from catch import probe
from catch.utils import longest_common_substring as lcf

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class PolyAFilter(BaseFilter):
    """Filter that removes probes with poly(A) or poly(T).
    """

    def __init__(self, length, mismatches, min_exact_length_to_consider=6):
        """
        Args:
            length/mismatches: remove probes that contain at least
                LENGTH 'A' bases in a row, tolerating up to MISMATCHES
                mismatches (and likewise for 'T' bases)
            min_exact_length_to_consider: only look for a stretch of 'A'
                or 'T' (according to length/mismatches) in a probe if that
                probe contains an exact stretch of 'A' or 'T' that is
                at least this length long. This is only meant to improve
                runtime, because the call to
                probe.Probe.longest_common_substring_length() is slow;
                this can reduce the number of times that needs to be called,
                but also result in false negatives. To always look
                for a stretch according only to length/mismatches, set
                the value of this argument to 0.
        """
        self.length = length
        self.mismatches = mismatches
        self.min_exact_length_to_consider = min_exact_length_to_consider

    def _filter(self, input):
        """Return a subset of the input probes.
        """
        if len(input) == 0:
            return input

        exact_a_stretch = 'A'*self.min_exact_length_to_consider
        exact_t_stretch = 'T'*self.min_exact_length_to_consider

        probe_len = len(input[0])
        for p in input:
            probe_len = max(probe_len, len(p))
        a_stretch = probe.Probe.from_str('A'*probe_len)
        t_stretch = probe.Probe.from_str('T'*probe_len)

        out = []
        for p in input:
            keep = True
            if exact_a_stretch in p.seq_str or exact_t_stretch in p.seq_str:
                for stretch in [a_stretch, t_stretch]:
                    lcf_len = p.longest_common_substring_length(
                        stretch, self.mismatches)
                    if lcf_len >= self.length:
                        # The stretch exceeds the limit
                        keep = False
                        break
            if keep:
                out += [p]
        return out
