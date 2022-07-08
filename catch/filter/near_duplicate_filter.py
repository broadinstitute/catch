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


class NearDuplicateFilter(BaseFilter):
    """Filter that removes near-duplicates using LSH.

    This constructs a concatenation of k hash functions, and does
    this multiple times so as to achieve a desired probability of
    reporting any probe as a near-duplicate of a queried probe. k
    can be a constant, and the number of (concatenated) hash functions
    to use is calculated to achieve the desired reporting probability.

    This sorts input probes by their multiplicity; therefore, the
    duplicate filter should *not* be run before this.
    """

    def __init__(self, k, reporting_prob=0.95):
        """
        Args:
            k: number of hash functions to draw from a family of
                hash functions for amplification; each hash function is then
                the concatenation (h_1, h_2, ..., h_k)
            reporting_prob: ensure that any probe within dist_thres of
                a queried probe is detected as such; this constructs
                multiple hash functions (each of which is a concatenation
                of k functions drawn from the family) to achieve this
                probability
        """
        self.k = k
        self.reporting_prob = reporting_prob

    def _filter(self, input):
        """Filter with an arbitrary LSH family.

        This performs near neighbor lookups using self.lsh_family. It only
        calls probes near-duplicates if their distance, according to
        self.dist_fn, is within self.dist_thres.

        Args:
            input: collection of probes to filter

        Returns:
            subset of input
        """
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
            self.dist_fn, self.reporting_prob)
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


class NearDuplicateFilterWithHammingDistance(NearDuplicateFilter):
    """Filter that removes near-duplicates according to Hamming distance.
    """

    def __init__(self, dist_thres, probe_length):
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
        """
        super().__init__(k=20)
        self.lsh_family = lsh.HammingDistanceFamily(probe_length)
        self.dist_thres = dist_thres

        def hamming_dist(a, b):
            # a and b are probe.Probe objects
            return a.mismatches(b)
        self.dist_fn = hamming_dist

    def _filter(self, input):
        """Filter with LSH using family that works with Hamming distance.

        Args:
            input: collection of probes to filter

        Returns:
            subset of input
        """
        return NearDuplicateFilter._filter(self, input)


class NearDuplicateFilterWithMinHash(NearDuplicateFilter):
    """Filter that removes near-duplicates using MinHash.
    """

    def __init__(self, dist_thres, kmer_size=10):
        """
        Args:
            dist_thres: only call two probes near-duplicates if their
                Jaccard distance (1 minus Jaccard similarity) is within
                this value; the Jaccard similarity is measured by treating
                each probe sequence as a set of k-mers and measuring
                the overlap of those k-mers
            kmer_size: the length of each k-mer to use with MinHash; note
                that this is *not* the same as self.k
        """
        super().__init__(k=3)
        self.lsh_family = lsh.MinHashFamily(kmer_size,
                use_fast_str_hash=True) # safe as long as not parallelized
        self.dist_thres = dist_thres

        def jaccard_dist(a, b):
            a_kmers = [a[i:(i + kmer_size)] for i in range(len(a) - kmer_size + 1)]
            b_kmers = [b[i:(i + kmer_size)] for i in range(len(b) - kmer_size + 1)]
            a_kmers = set(a_kmers)
            b_kmers = set(b_kmers)
            jaccard_sim = float(len(a_kmers & b_kmers)) / len(a_kmers | b_kmers)
            return 1.0 - jaccard_sim
        self.dist_fn = jaccard_dist

    def _filter(self, input):
        """Filter with LSH using MinHash family.

        Args:
            input: collection of probes to filter

        Returns:
            subset of input
        """
        return NearDuplicateFilter._filter(self, input)

