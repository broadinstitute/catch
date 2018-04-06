"""Classes and methods for applying locality-sensitive hashing.
"""

import logging
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

    def __init__(self, k, family):
        self.k = k
        self.family = family
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

