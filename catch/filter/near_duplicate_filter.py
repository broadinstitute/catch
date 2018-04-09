"""Removes near-duplicates from an input list of probes.

This acts as a filter on the probes by removing ones that are
near-duplicates of another using LSH. There might be near-duplicates
in the output that are not detected, but every near-duplicate removed
should indeed be a near-duplicate as defined by the given criteria.
"""

from collections import defaultdict
import math
import operator

from catch.filter.base_filter import BaseFilter
from catch.utils import lsh

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class NearDuplicateFilterWithHammingDistance(BaseFilter):
    """Filter that removes near-duplicates according to Hamming distance.

    This constructs a concatenation of k hash functions, and does
    this multiple times so as to achieve a desired probability of
    reporting any probe as a near-duplicate of a queried probe. k
    can be a constant, and the number of (concatenated) hash functions
    to use is calculated to achieve the desired reporting probability.

    This sorts input probes by their multiplicity; therefore, the
    duplicate filter should *not* be run before this.
    """

    def __init__(self, dist_thres, probe_length, k=20, reporting_prob=0.95):
        """
        Args:
            dist_thres: only call two probes near-duplicates if their
                Hamming distance is within this value; this should be
                equal to or commensurate with (but not greater than)
                the number of mismatches at/below which a probe is
                considered to hybridize to a target sequence so that
                candidate probes further apart than this value are not
                collapsed as near-duplicates
            probe_length: length of probes
            k: number of hash functions to draw from a family of
                hash functions for the Hamming distance for amplification;
                each hash function is then the concatenation
                (h_1, h_2, ..., h_k)
            reporting_prob: ensure that any probe within dist_thres of
                a queried probe is detected as such; this constructs
                multiple hash functions (each of which is a concatenation
                of k functions drawn from the family) to achieve this
                probability
        """
        self.lsh_family = lsh.HammingDistanceFamily(probe_length)
        self.dist_thres = dist_thres
        self.k = k
        self.reporting_prob = reporting_prob

    def _filter(self, input):
        """Filter with LSH using family that works with Hamming distance.

        Args:
            input: collection of probes to filter

        Returns:
            subset of input
        """
        def hamming_dist_fn(a, b):
            # a and b are Probe objects
            return a.mismatches(b)

        # Sort the probes by their mulitiplicity (descending)
        occurrences = defaultdict(int)
        for p in input:
            occurrences[p] += 1
        input_sorted = [p for p, count in
            sorted(occurrences.items(), key=operator.itemgetter(1),
                   reverse=True)]

        # Remove exact duplicates from the input
        input = list(set(input))

        # Construct a collection of hash tables for looking up
        # near neighbors of each probe
        nnl = lsh.NearNeighborLookup(self.lsh_family, self.k, self.dist_thres,
            hamming_dist_fn, self.reporting_prob)
        nnl.add(input)

        # Iterate through all probes in order; for each p, remove others
        # that are near-duplicates (neighbors) of p. Since we iterate
        # in sorted order by multiplicity, the ones that hit more targets
        # should appear earlier and will be included in the filtered output
        to_include = set()
        to_exclude = set()
        for p in input_sorted:
            # p should not have already been included because input_sorted
            # should not contain duplicates
            assert p not in to_include

            if p in to_exclude:
                # p is already being filtered out
                continue

            # Include p in the output and exclude all near-duplicates of it
            to_include.add(p)
            for near_dup in nnl.query(p):
                if near_dup not in to_include:
                    to_exclude.add(near_dup)

        # Check that every probe is either included or excluded and
        # that none are both included and excluded
        assert len(to_include | to_exclude) == len(input_sorted)
        assert len(to_include & to_exclude) == 0

        return list(to_include)

