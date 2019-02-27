"""Functions for clustering sequences before input.

This includes computing a distance matrix using MinHash, and
clustering that matrix.
"""

from collections import defaultdict
import logging
import operator

import numpy as np
from scipy.cluster import hierarchy

from catch.utils import lsh

__author__ = 'Hayden Metsky <hayden@mit.edu>'

logger = logging.getLogger(__name__)


def make_signatures_with_minhash(family, seqs):
    """Construct a signature using MinHash for each sequence.
    Args:
        family: lsh.MinHashFamily object
        seqs: dict mapping sequence header to sequences
    Returns:
        dict mapping sequence header to signature
    """
    # Construct a single hash function; use the same for all sequences
    h = family.make_h()

    signatures = {}
    for name, seq in seqs.items():
        signatures[name] = h(seq)
    return signatures


def create_condensed_dist_matrix(n, dist_fn):
    """Construct a 1d condensed distance matrix for scipy.
    Args:
        n: number of elements whose pairwise distances to store in the
            matrix
        dist_fn: function such that dist_fn(i, j) gives the distance
            between i and j, for all i<j<n
    Returns:
        condensed 1d distance matrix for input to scipy functions
    """
    def idx(i, j):
        # Compute index in 1d vector for pair (i, j)
        return int((-1 * i*i)/2 + i*n - 3*i/2 + j - 1)

    dist_matrix_len = int(n*(n-1)/2)
    dist_matrix = np.zeros(dist_matrix_len)

    num_pairs = n*(n-1) / 2
    pair_counter = 0
    logger.debug(("Condensed distance matrix has %d entries"), num_pairs)
    for j in range(n):
        for i in range(j):
            if (pair_counter + 1) % 100000 == 0:
                logger.debug(("Computing condensed distance matrix entry "
                    "%d of %d"), pair_counter + 1, num_pairs)
            pair_counter += 1
            dist_matrix[idx(i, j)] = dist_fn(i, j)

    return dist_matrix


def cluster_from_dist_matrix(dist_matrix, threshold):
    """Use scipy to cluster a distance matrix.
    Args:
        dist_matrix: distance matrix, represented in scipy's 1d condensed form
        threshold: maximum inter-cluster distance to merge clusters (higher
            results in fewer clusters)
    Returns:
        list c such that c[i] is a collection of all the observations
        (whose pairwise distances are indexed in dist) in the i'th
        cluster, in sorted order by descending cluster size
    """
    linkage = hierarchy.linkage(dist_matrix, method='average')
    clusters = hierarchy.fcluster(linkage, threshold, criterion='distance')

    # clusters are numbered starting at 1, but base the count on
    # first_clust_num just in case this changes
    first_clust_num = min(clusters)
    num_clusters = max(clusters) + 1 - first_clust_num
    elements_in_cluster = defaultdict(list)
    for i, clust_num in enumerate(clusters):
        elements_in_cluster[clust_num].append(i)
    cluster_sizes = {c: len(elements_in_cluster[c])
                     for c in range(first_clust_num,
                                    num_clusters + first_clust_num)}

    elements_in_cluster_sorted = []
    for clust_num, _ in sorted(cluster_sizes.items(),
            key=operator.itemgetter(1), reverse=True):
        elements_in_cluster_sorted += [elements_in_cluster[clust_num]]
    return elements_in_cluster_sorted


def cluster_with_minhash_signatures(seqs, k=12, N=100, threshold=0.1):
    """Cluster sequences based on their MinHash signatures.
    Args:
        seqs: dict mapping sequence header to sequences
        k: k-mer size to use for k-mer hashes (smaller is likely more
            sensitive for divergent genomes, but may lead to false positives
            in determining which genomes are close)
        N: number of hash values to use in a signature (higher is slower for
            clustering, but likely more sensitive for divergent genomes)
        threshold: maximum inter-cluster distance to merge clusters, in
            average nucleotide dissimilarity (1-ANI, where ANI is
            average nucleotide identity); higher results in fewer
            clusters
    Returns:
        list c such that c[i] gives a collection of sequence headers
        in the same cluster, and the clusters in c are sorted
        in descending order of size
    """
    num_seqs = len(seqs)

    logger.info(("Producing signatures of %d sequences"), num_seqs)
    family = lsh.MinHashFamily(k, N=N)
    signatures_map = make_signatures_with_minhash(family, seqs)

    # Map each sequence header to an index (0-based), and get
    # the signature for the corresponding index
    seq_headers = []
    signatures = []
    for name, seq in seqs.items():
        seq_headers += [name]
        signatures += [signatures_map[name]]

    # Eq. 4 of the Mash paper (Ondov et al. 2016) shows that the
    # Mash distance, which is shown to be closely related to 1-ANI, is:
    #  D = (-1/k) * ln(2*j/(1+j))
    # where j is a Jaccard similarity. Solving for j:
    #  j = 1/(2*exp(k*D) - 1)
    # So, for a desired distance D in terms of 1-ANI, the corresponding
    # Jaccard distance is:
    #  1.0 - 1/(2*exp(k*D) - 1)
    # We can use this to calculate a clustering threshold in terms of
    # Jaccard distance
    jaccard_dist_threshold = 1.0 - 1.0/(2.0*np.exp(k*threshold) - 1)

    def jaccard_dist(i, j):
        # Return estimated Jaccard dist between signatures at
        # index i and index j
        return family.estimate_jaccard_dist(
            signatures[i], signatures[j])

    logger.info(("Creating condensed distance matrix of %d sequences"), num_seqs)
    dist_matrix = create_condensed_dist_matrix(num_seqs, jaccard_dist)
    logger.info(("Clustering %d sequences at Jaccard distance threshold of %f"),
            num_seqs, jaccard_dist_threshold)
    clusters = cluster_from_dist_matrix(dist_matrix,
        jaccard_dist_threshold)

    seqs_in_cluster = []
    for cluster_idxs in clusters:
        seqs_in_cluster += [[seq_headers[i] for i in cluster_idxs]]
    return seqs_in_cluster

