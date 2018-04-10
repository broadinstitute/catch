"""Classes and methods for applying locality-sensitive hashing.
"""

from collections import defaultdict
import logging
import math
import random

__author__ = 'Hayden Metsky <hayden@mit.edu>'

logger = logging.getLogger(__name__)


class HammingDistanceFamily:
    """An LSH family that works with Hamming distance by sampling bases."""

    def __init__(self, dim):
        self.dim = dim

    def make_h(self):
        """Construct a random hash function for this family.

        Returns:
            hash function
        """
        i = random.randint(0, self.dim - 1)
        def h(x):
            assert len(x) == self.dim
            return x[i]
        return h

    def P1(self, dist):
        """Calculate lower bound on probability of collision for nearby sequences.

        Args:
            dist: Hamming distance; suppose two sequences are within this
                distance of each other

        Returns:
            lower bound on probability that two sequences (e.g., probes) hash
            to the same value if they are within dist of each other
        """
        return 1.0 - float(dist) / float(self.dim)


class HashConcatenation:
    """Concatenated hash functions (AND constructions)."""

    def __init__(self, family, k):
        self.family = family
        self.k = k
        self.hs = [family.make_h() for _ in range(k)]

    def g(self, x):
        """Evaluate random hash functions and concatenate the result.

        Args:
            x: point (e.g., probe) on which to evaluate hash functions

        Returns:
            concatenation of the result of the self.k random hash functions
            evaluated at x
        """
        return tuple([h(x) for h in self.hs])


class NearNeighborLookup:
    """Support for approximate near neighbor lookups.

    This implements the R-near neighbor reporting problem described in
    Andoni and Indyk 2008.
    """

    def __init__(self, family, k, dist_thres, dist_fn, reporting_prob):
        """
        This selects a number of hash tables (defined as L in the above
        reference) according to the strategy it outlines: we want any
        neighbor (within dist_thres) of a query to be reported with
        probability at least reporting_prob; thus, the number of
        tables should be [log_{1 - (P1)^k} (1 - reporting_prob)]. In
        the above reference, delta is 1.0 - reporting_prob.

        Args:
            family: object giving family of hash functions
            k: number of hash functions from family to concatenate
            dist_thres: consider any two objects within this threshold
                of each other to be neighbors
            dist_fn: function f(a, b) that calculates the distance between
                a and b, to compare against dist_thres
            reporting_prob: report any neighbor of a query with
                probability at least equal to this
        """
        self.family = family
        self.k = k
        self.dist_thres = dist_thres
        self.dist_fn = dist_fn

        P1 = self.family.P1(dist_thres)
        if P1 == 1.0:
            # dist_thres might be 0, and any number of hash tables can
            # satisfy the reporting probability
            self.num_tables = 1
        else:
            self.num_tables = math.log(1.0 - reporting_prob, 1.0 - math.pow(P1, k))
            self.num_tables = int(math.ceil(self.num_tables))

        # Setup self.num_tables hash tables, each with a corresponding
        # function for hashing into it (the functions are concatenations
        # of k hash functions from the given family)
        self.hashtables = []
        self.hashtables_g = []
        for j in range(self.num_tables):
            g = HashConcatenation(self.family, self.k)
            self.hashtables += [defaultdict(list)]
            self.hashtables_g += [g]

    def add(self, pts):
        """Insert given points into each of the hash tables.

        Args:
            pts: collection of points (e.g., probes) to add to the hash
                tables
        """
        for j in range(self.num_tables):
            ht = self.hashtables[j]
            g = self.hashtables_g[j].g
            for p in pts:
                ht[g(p)].append(p)

    def query(self, q):
        """Find neighbors of a query point.

        Args:
            q: query point (e.g., probe)

        Returns:
            collection of stored points that are within self.dist_thres of
            q; all returned points are within this distance, but the
            returned points might not include all that are
        """
        neighbors = set()
        for j in range(self.num_tables):
            ht = self.hashtables[j]
            g = self.hashtables_g[j].g
            for p in ht[g(q)]:
                if self.dist_fn(q, p) <= self.dist_thres:
                    neighbors.add(p)
        return neighbors

