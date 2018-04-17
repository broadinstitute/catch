"""Tests for lsh module.
"""

import random
import unittest

from catch.utils import lsh

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class TestHammingDistanceFamily(unittest.TestCase):
    """Tests family of hash functions for Hamming distance.
    """

    def setUp(self):
        # Set a random seed so hash functions are always the same
        random.seed(0)

        self.family = lsh.HammingDistanceFamily(20)

    def test_identical(self):
        a = 'ATCGATATGGGCACTGCTAT'
        b = str(a)

        # Identical strings should hash to the same value
        h1 = self.family.make_h()
        self.assertEqual(h1(a), h1(b))
        h2 = self.family.make_h()
        self.assertEqual(h2(a), h2(b))

    def test_similar(self):
        a = 'ATCGATATGGGCACTGCTAT'
        b = 'ATCGACATGGGCACTGGTAT'

        # a and b should probably collide
        collision_count = 0
        for i in range(10):
            h = self.family.make_h()
            if h(a) == h(b):
                collision_count += 1
        self.assertGreater(collision_count, 8)

    def test_not_similar(self):
        a = 'ATCGATATGGGCACTGCTAT'
        b = 'AGTTGTCACCCTTGACGATA'

        # a and b should rarely collide
        collision_count = 0
        for i in range(10):
            h = self.family.make_h()
            if h(a) == h(b):
                collision_count += 1
        self.assertLess(collision_count, 2)

    def test_collision_prob(self):
        # Collision probability for 2 mismatches should be
        # 1 - 2/20
        self.assertEqual(self.family.P1(2), 0.9)


class TestMinHashFamily(unittest.TestCase):
    """Tests family of hash functions for MinHash.
    """

    def setUp(self):
        # Set a random sseed so hash functions are always the same
        random.seed(0)

        self.family = lsh.MinHashFamily(3)

    def test_identical(self):
        a = 'ATCGATATGGGCACTGCTAT'
        b = str(a)

        # Identical strings should hash to the same value
        h1 = self.family.make_h()
        self.assertEqual(h1(a), h1(b))
        h2 = self.family.make_h()
        self.assertEqual(h2(a), h2(b))

    def test_similar(self):
        a = 'ATCGATATGGGCACTGCTATGTAGCGC'
        b = 'ATCGACATGGGCACTGGTATGTAGCGC'

        # a and b should probably collide; the Jaccard similarity
        # of a and b is ~67% (with 3-mers being the elements that
        # make up each sequence) so they should collide with that
        # probability (check that it is >60%)
        collision_count = 0
        for i in range(100):
            h = self.family.make_h()
            if h(a) == h(b):
                collision_count += 1
        self.assertGreater(collision_count, 60)

    def test_not_similar(self):
        a = 'ATCGATATGGGCACTGCTAT'
        b = 'AGTTGTCACCCTTGACGATA'

        # a and b should rarely collide
        collision_count = 0
        for i in range(100):
            h = self.family.make_h()
            if h(a) == h(b):
                collision_count += 1
        self.assertLess(collision_count, 30)

    def test_collision_prob(self):
        # Collision probability for two sequences with a Jaccard
        # distance of 0.2 should be 0.8
        self.assertEqual(self.family.P1(0.2), 0.8)


class TestHammingHashConcatenation(unittest.TestCase):
    """Tests concatenations of hash functions with Hamming distance.
    """

    def setUp(self):
        # Set a random seed so hash functions are always the same
        random.seed(0)

        self.family = lsh.HammingDistanceFamily(20)
        self.G = lsh.HashConcatenation(self.family, 100)

    def test_identical(self):
        # Identical a and b should collide even with large k
        a = 'ATCGATATGGGCACTGCTAT'
        b = str(a)
        self.assertEqual(self.G.g(a), self.G.g(b))

    def test_similar(self):
        # Similar (but not identical) a and b should rarely
        # collide when k is large
        a = 'ATCGATATGGGCACTGCTAT'
        b = 'ATCGACATGGGCACTGGTAT'

        collision_count = 0
        for i in range(10):
            if self.G.g(a) == self.G.g(b):
                collision_count += 1
        self.assertLess(collision_count, 2)

    def test_not_similar(self):
        a = 'ATCGATATGGGCACTGCTAT'
        b = 'AGTTGTCACCCTTGACGATA'

        # a and b should rarely collide
        collision_count = 0
        for i in range(10):
            if self.G.g(a) == self.G.g(b):
                collision_count += 1
        self.assertLess(collision_count, 2)


class TestHammingNearNeighborLookup(unittest.TestCase):
    """Tests approximate near neighbor lookups with Hamming distance."""

    def setUp(self):
        # Set a random seed so hash functions are always the same
        random.seed(0)

        self.family = lsh.HammingDistanceFamily(20)
        self.dist_thres = 5
        def f(a, b):
            assert len(a) == len(b)
            return sum(1 for i in range(len(a)) if a[i] != b[i])
        self.dist_fn = f

    def test_varied_k(self):
        a = 'ATCGATATGGGCACTGCTAT'
        b = str(a)  # identical to a
        c = 'ATCGACATGGGCACTGGTAT'  # similar to a
        d = 'AGTTGTCACCCTTGACGATA'  # not similar to a
        e = 'AGTTGTCACCCTTGACGATA'  # similar to d

        for k in [2, 5, 10]:
            nnl = lsh.NearNeighborLookup(self.family, k, self.dist_thres,
                self.dist_fn, 0.95)
            nnl.add([a, b, c, d])

            # b and c are within self.dist_thres of a, so only these
            # should be returned (along with a); note that since
            # a==b, {a,b,c}=={a,c}=={b,c} and nnl.query(a) returns
            # a set, which will be {a,c} or {b,c}
            self.assertCountEqual(nnl.query(a), {a, b, c})

            # Although e was not added, a query for it should return d
            self.assertCountEqual(nnl.query(e), {d})


class TestMinHashNearNeighborLookup(unittest.TestCase):
    """Tests approximate near neighbor lookups with MinHash."""

    def setUp(self):
        # Set a random seed so hash functions are always the same
        random.seed(0)

        kmer_size = 3
        self.family = lsh.MinHashFamily(kmer_size)
        self.dist_thres = 0.5
        def f(a, b):
            a_kmers = [a[i:(i + kmer_size)] for i in range(len(a) - kmer_size + 1)]
            b_kmers = [b[i:(i + kmer_size)] for i in range(len(b) - kmer_size + 1)]
            a_kmers = set(a_kmers)
            b_kmers = set(b_kmers)
            jaccard_sim = float(len(a_kmers & b_kmers)) / len(a_kmers | b_kmers)
            return 1.0 - jaccard_sim
        self.dist_fn = f

    def test_varied_k(self):
        a = 'ATCGATATGGGCACTGCTAT'
        b = str(a)  # identical to a
        c = 'ATCGACATGGGCACTGGTAT'  # similar to a
        d = 'AGTTGTCACCCTTGACGATA'  # not similar to a
        e = 'AGTTGTCACCCTTGACGATA'  # similar to d

        for k in [2, 5, 10]:
            nnl = lsh.NearNeighborLookup(self.family, k, self.dist_thres,
                self.dist_fn, 0.95)
            nnl.add([a, b, c, d])

            # b and c are within self.dist_thres of a, so only these
            # should be returned (along with a); note that since
            # a==b, {a,b,c}=={a,c}=={b,c} and nnl.query(a) returns
            # a set, which will be {a,c} or {b,c}
            self.assertCountEqual(nnl.query(a), {a, b, c})

            # Although e was not added, a query for it should return d
            self.assertCountEqual(nnl.query(e), {d})

