"""Tests for cluster module.
"""

import logging
import unittest

import numpy as np
import scipy

from catch.utils import cluster

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class TestClusterFromMatrix(unittest.TestCase):
    """Test basic clustering functions."""

    def setUp(self):
        # Disable logging
        logging.disable(logging.WARNING)

    def test_create_condensed_dist_matrix(self):
        # Have 3 elements: 0 and 1 are similar, 0 and 2 are
        # very dissimilar, and 1 and 2 are similar
        n = 3
        dist_matrix_2d = np.array(
            [[0,   1, 100],
             [1,   0, 2  ],
             [100, 2, 0  ]])
        def dist_fn(i, j):
            return dist_matrix_2d[i][j]

        condensed = cluster.create_condensed_dist_matrix(n, dist_fn)

        # Use scipy to create a condensed matrix, which is possible
        # since we already have the 2d matrix available
        scipy_condensed = scipy.spatial.distance.squareform(dist_matrix_2d)

        self.assertTrue(np.array_equal(condensed, scipy_condensed))

    def test_cluster_from_dist_matrix_1(self):
        # Have 3 elements: 0 and 1 are similar, 0 and 2 are
        # very dissimilar, and 1 and 2 are similar
        n = 3
        dists = {(0, 1): 1, (0, 2): 100, (1, 2): 2}
        def dist_fn(i, j):
            return dists[(i, j)]
        dist_matrix = cluster.create_condensed_dist_matrix(n, dist_fn)

        clusters = cluster.cluster_from_dist_matrix(dist_matrix, 10)
        self.assertEqual(clusters, [[0, 1], [2]])

    def test_cluster_from_dist_matrix_2(self):
        # Have 3 elements: 0 and 1 are very dissimilar similar, 0 and 2 are
        # very dissimilar, and 1 and 2 are similar
        n = 3
        dists = {(0, 1): 100, (0, 2): 100, (1, 2): 1}
        def dist_fn(i, j):
            return dists[(i, j)]
        dist_matrix = cluster.create_condensed_dist_matrix(n, dist_fn)

        clusters = cluster.cluster_from_dist_matrix(dist_matrix, 10)
        self.assertEqual(clusters, [[1, 2], [0]])

    def test_cluster_from_dist_matrix_3(self):
        # Have 3 elements: 0 and 1 are very dissimilar similar, 0 and 2 are
        # very dissimilar, and 1 and 2 are similar
        # Have 3 elements: all very dissimilar
        n = 3
        dists = {(0, 1): 20, (0, 2): 30, (1, 2): 40}
        def dist_fn(i, j):
            return dists[(i, j)]
        dist_matrix = cluster.create_condensed_dist_matrix(n, dist_fn)

        clusters = cluster.cluster_from_dist_matrix(dist_matrix, 10)
        self.assertEqual(sorted(clusters), [[0], [1], [2]])

    def tearDown(self):
        # Re-enable logging
        logging.disable(logging.NOTSET)


class TestClusterWithMinHashSignatures(unittest.TestCase):
    """Test cluster_with_minhash_signatures() function."""

    def setUp(self):
        # Disable logging
        logging.disable(logging.WARNING)

    def test_simple(self):
        seqs = {'a': 'AT'*500,
                'b': 'CG'*500,
                'c': 'AT'*500,
                'd': 'TA'*500,
                'e': 'TT'*500,
                'f': 'CG'*500}
        # The expected clusters should be: [a, c, d], [b, f], [e]

        clusters = cluster.cluster_with_minhash_signatures(seqs)
        self.assertEqual(len(clusters), 3)
        self.assertEqual(sorted(clusters[0]), ['a', 'c', 'd'])
        self.assertEqual(sorted(clusters[1]), ['b', 'f'])
        self.assertEqual(sorted(clusters[2]), ['e'])

    def tearDown(self):
        # Re-enable logging
        logging.disable(logging.NOTSET)

