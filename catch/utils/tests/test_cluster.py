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

    def create_condensed_dist_matrix(self, n, dist_fn):
        # Run cluster.create_condensed_dist_matrix() for multiple values
        # of num_processes; make sure the returned results are all equal
        result = None
        for num_processes in [1, 2, 4, 8]:
            condensed = cluster.create_condensed_dist_matrix(n, dist_fn,
                    num_processes=num_processes)
            if result is None:
                result = condensed
            else:
                self.assertTrue(np.array_equal(result, condensed))
        return result

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

        condensed = self.create_condensed_dist_matrix(n, dist_fn)

        # Use scipy to create a condensed matrix, which is possible
        # since we already have the 2d matrix available
        scipy_condensed = scipy.spatial.distance.squareform(dist_matrix_2d)

        self.assertTrue(np.array_equal(condensed, scipy_condensed))

    def test_cluster_hierarchically_from_dist_matrix_1(self):
        # Have 3 elements: 0 and 1 are similar, 0 and 2 are
        # very dissimilar, and 1 and 2 are similar
        n = 3
        dists = {(0, 1): 1, (0, 2): 100, (1, 2): 2}
        def dist_fn(i, j):
            return dists[(i, j)]
        dist_matrix = self.create_condensed_dist_matrix(n, dist_fn)

        clusters = cluster.cluster_hierarchically_from_dist_matrix(dist_matrix, 10)
        self.assertEqual(clusters, [[0, 1], [2]])

    def test_cluster_hierarchically_from_dist_matrix_2(self):
        # Have 3 elements: 0 and 1 are very dissimilar similar, 0 and 2 are
        # very dissimilar, and 1 and 2 are similar
        n = 3
        dists = {(0, 1): 100, (0, 2): 100, (1, 2): 1}
        def dist_fn(i, j):
            return dists[(i, j)]
        dist_matrix = self.create_condensed_dist_matrix(n, dist_fn)

        clusters = cluster.cluster_hierarchically_from_dist_matrix(dist_matrix, 10)
        self.assertEqual(clusters, [[1, 2], [0]])

    def test_cluster_hierarchically_from_dist_matrix_3(self):
        # Have 3 elements: 0 and 1 are very dissimilar similar, 0 and 2 are
        # very dissimilar, and 1 and 2 are similar
        # Have 3 elements: all very dissimilar
        n = 3
        dists = {(0, 1): 20, (0, 2): 30, (1, 2): 40}
        def dist_fn(i, j):
            return dists[(i, j)]
        dist_matrix = self.create_condensed_dist_matrix(n, dist_fn)

        clusters = cluster.cluster_hierarchically_from_dist_matrix(dist_matrix, 10)
        self.assertEqual(sorted(clusters), [[0], [1], [2]])

    def tearDown(self):
        # Re-enable logging
        logging.disable(logging.NOTSET)


class TestSimpleClusterWithThreshold(unittest.TestCase):
    """Test find_connected_components() function."""

    def setUp(self):
        # Disable logging
        logging.disable(logging.WARNING)

        # Make a distance function
        def make_dist_fn(seqs):
            def dist_fn(i, j):
                # Distance is the number of characters separating
                #   letter at index i from letter at index j
                return abs(ord(seqs[i]) - ord(seqs[j]))
            return dist_fn
        self.make_dist_fn = make_dist_fn


    def test_simple(self):
        # If any letter is next to another in the alphabet,
        #   consider them adjacent in a graph
        seqs = ['a', 'b', 'x', 'm', 'c', 'o', 'n', 'z', 'w', 'y', 'v', 'd', 'k']

        # clusters = {v,w,x,y,z}, {a,b,c,d}, {m,n,o}, {k}
        expected = [[2,7,8,9,10], [0,1,4,11], [3,5,6], [12]]

        dist = self.make_dist_fn(seqs)

        # Cluster with threshold of 1
        ccs = cluster.find_connected_components(len(seqs), dist, 1)
        self.assertEqual(ccs, expected)

    def test_one_connected_component(self):
        # If any letter is next to another in the alphabet,
        #   consider them adjacent in a graph
        seqs = ['a', 'c', 'b', 'c']

        # clusters = {a,b,c,c}
        expected = [[0,1,2,3]]

        dist = self.make_dist_fn(seqs)

        # Cluster with threshold of 1
        ccs = cluster.find_connected_components(len(seqs), dist, 1)
        self.assertEqual(ccs, expected)

    def test_all_different_connected_components(self):
        # If any letter is next to another in the alphabet,
        #   consider them adjacent in a graph
        seqs = ['a', 'z', 'm']

        # clusters = {a}, {z}, {m}
        expected = [[0], [1], [2]]

        dist = self.make_dist_fn(seqs)

        # Cluster with threshold of 1
        ccs = cluster.find_connected_components(len(seqs), dist, 1)

        # Use assertCountEqual because order may vary
        self.assertCountEqual(ccs, expected)

    def test_empty_input(self):
        # If any letter is next to another in the alphabet,
        #   consider them adjacent in a graph
        seqs = []

        expected = []

        dist = self.make_dist_fn(seqs)

        # Cluster with threshold of 1
        ccs = cluster.find_connected_components(len(seqs), dist, 1)

        self.assertEqual(ccs, expected)

    def tearDown(self):
        # Re-enable logging
        logging.disable(logging.NOTSET)


class TestClusterWithMinHashSignatures(unittest.TestCase):
    """Test cluster_with_minhash_signatures() function."""

    def setUp(self):
        # Disable logging
        logging.disable(logging.WARNING)

        self.seqs = {'a': 'AT'*500,
                     'b': 'CG'*500,
                     'c': 'AT'*500,
                     'd': 'TA'*500,
                     'e': 'TT'*500,
                     'f': 'CG'*500,
                     'g': 'GC'*500,
                     'h': 'AT'*505}
        # The expected clusters should be: [a, c, d, h], [b, f, g], [e]

    def test_wth_simple_cluster_method(self):
        clusters = cluster.cluster_with_minhash_signatures(self.seqs,
                cluster_method='simple')
        self.assertEqual(len(clusters), 3)
        self.assertEqual(sorted(clusters[0]), ['a', 'c', 'd', 'h'])
        self.assertEqual(sorted(clusters[1]), ['b', 'f', 'g'])
        self.assertEqual(sorted(clusters[2]), ['e'])

    def test_wth_hierarchical_cluster_method(self):
        clusters = cluster.cluster_with_minhash_signatures(self.seqs,
                cluster_method='hierarchical')
        self.assertEqual(len(clusters), 3)
        self.assertEqual(sorted(clusters[0]), ['a', 'c', 'd', 'h'])
        self.assertEqual(sorted(clusters[1]), ['b', 'f', 'g'])
        self.assertEqual(sorted(clusters[2]), ['e'])

    def tearDown(self):
        # Re-enable logging
        logging.disable(logging.NOTSET)

