"""Classes and methods for applying locality-sensitive hashing.
"""

from collections import defaultdict
import hashlib
import heapq
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


class MinHashFamily:
    """An LSH family that works by taking the minimum permutation of
    k-mers in a string/sequence (MinHash).

    See (Broder et al. 1997) and (Andoni and Indyk 2008) for details.

    The signature produced here also has similarity to Mash (Ondov et al.
    2016).
    """

    def __init__(self, kmer_size, N=1, use_fast_str_hash=False):
        """
        Args:
            kmer_size: length of each k-mer to hash
            N: represent the signature of a sequence using hash values of
                the N k-mers in the sequence that have the smallest hash
                values
            use_fast_str_hash: if True, use a faster (~10x faster) inner
                hash function for strings (k-mers) that is not deterministic
                and varies across Python processes; this may not be
                appropriate if the family is shared across processes
        """
        self.kmer_size = kmer_size
        self.N = N
        self.use_fast_str_hash = use_fast_str_hash

    def make_h(self):
        """Construct a random hash function for this family.

        Here, we treat a sequence as being a set of k-mers. We calculate
        a hash value for each k-mer and the hash function on the sequence
        returns the N smallest of these in sorted order.

        Returns:
            hash function
        """
        # First construct a random hash function for a k-mer that
        # is a universal hash function (effectively a "permutation"
        # of the k-mer)
        # Constrain all values to be in [0, 2^31 - 1] to have a bound
        # on the output of the universal hash function; this upper bound
        # is nice because it is also a prime, so we can simply work
        # modulo (2^31 - 1)
        p = 2**31 - 1
        # Let the random hash function be:
        #   (a*x + b) mod p
        # for random integers a, b (a in [1, p] and b in [0, p])
        a = random.randint(1, p)
        b = random.randint(0, p)
        if self.use_fast_str_hash:
            def kmer_hash(x):
                # Hash a k-mer x with the random hash function
                # hash(..) uses a random seed in Python 3.3+, so its output
                #   varies across Python processes; thus, this may not be
                #   suitable if used across processes
                x_hash = abs(hash(x))
                return (a * x_hash + b) % p
        else:
            def kmer_hash(x):
                # Hash a k-mer x with the random hash function
                # hashlib.md5(..) gives a deterministic hash value of the k-mer
                #   but is ~10x slower than hash(..)
                x_hash = int(hashlib.md5(x.encode('utf-8')).hexdigest(), 16)
                return (a * x_hash + b) % p

        # TODO Consider only allowing fully unambiguous k-mers to be
        #   in the signature, as a check below in kmer_hashes()

        def h(s):
            # For a string/sequence s, have the MinHash function be the minimum
            # N hashes, in sorted order, over all the k-mers in it
            assert self.kmer_size <= len(s)
            if self.kmer_size >= len(s) / 2:
                logger.warning(("The k-mer size %d is large (> (1/2)x) "
                    "compared to the size of a sequence to hash (%d), which "
                    "might make it difficult for MinHash to find similar "
                    "sequence"), self.kmer_size, len(s))
            num_kmers = len(s) - self.kmer_size + 1
            if num_kmers < self.N:
                raise ValueError(("The number of k-mers (%d) in the given "
                    "sequence is too small to produce a signature of "
                    "size %d; try setting --small-seq-skip") %
                    (num_kmers, self.N))

            def kmer_hashes():
                for i in range(num_kmers):
                    kmer = s[i:(i + self.kmer_size)]
                    yield kmer_hash(kmer)
            if self.N == 1:
                # Speed the special case where self.N == 1
                return tuple([min(kmer_hashes())])
            else:
                return tuple(sorted(heapq.nsmallest(self.N, kmer_hashes())))
        return h

    def P1(self, dist):
        """Calculate lower bound on probability of collision for nearby sequences.

        Args:
            dist: Jaccard distance (1 minus Jaccard similarity); suppose
                two sequences are within this distance of each other. The
                Jaccard similarity can be thought of as the overlap in k-mers
                between the two sequences

        Returns:
            lower bound on probability that two sequences (e.g., probes) hash
            to the same value if they are within dist of each other
        """
        # With MinHash, the collision probability is the Jaccard similarity
        return 1.0 - dist

    def estimate_jaccard_dist(self, hA, hB):
        """Estimate Jaccard distance between two MinHash signatures.
        Args:
            hA: signature output by h(A) for a hash function h given by
                make_h() and a sequence A (specifically, an ordered list of
                hash values of k-mers in A)
            hB: signature output by h(B) for the same hash function used
                to generate hA, and a sequence B (specifically, an ordered
                list of hash values of k-mers in B)
        Returns:
            estimate of Jaccard distance between signatures A and B
        """
        # Let X = h(hA \union hB). This is equivalent to h(A \union B)
        # where A and B, here, represent the set of k-mers in their
        # respective sequences. X is a random sample of the hashes of
        # the k-mers in (A \union B), and (Y = X \intersect hA \intersect hB)
        # are the hashes in both X and (A \intersect B). So |Y|/|X| is
        # an estimate of the Jaccard similarity. Since X has self.N
        # elements (it is the output of h(.)), this is |Y|/self.N

        # Iterate over hA and hB simultaneously, to count the number of
        # hash values in common between them, among the first self.N of
        # them
        # This is possible because hA and hB are in sorted order
        hA_i, hB_i = 0, 0
        intersect_count = 0
        union_count = 0
        while hA_i < len(hA) and hB_i < len(hB):
            if union_count == self.N:
                # We found all the hash values in X, as defined above
                break
            elif hA[hA_i] < hB[hB_i]:
                hA_i += 1
                union_count += 1
            elif hA[hA_i] > hB[hB_i]:
                hB_i += 1
                union_count += 1
            else:
                # We found a new unique hash value in common:
                # hA[hA_i] == hB[hB_i]
                intersect_count += 1
                union_count += 1
                hA_i += 1
                hB_i += 1

        # It should be that union_count == self.N, except perhaps in
        # some edge cases (e.g., if a hash value is repeated in a signature)
        # And |Y| == intersect_count, where Y is defined above
        similarity = float(intersect_count) / union_count
        return 1.0 - similarity


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

