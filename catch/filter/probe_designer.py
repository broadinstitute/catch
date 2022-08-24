"""Designs probes with a filtering approach.
"""

import itertools
import logging

from catch.filter import candidate_probes
from catch import genome
from catch.utils import cluster

__author__ = 'Hayden Metsky <hayden@mit.edu>'

logger = logging.getLogger(__name__)


class ProbeDesigner:
    """Generates a set of candidate probes and filters them.

    There is an (ordered) list of one or more provided filters.
    """

    def __init__(self, genomes, filters, probe_length,
            probe_stride, allow_small_seqs=None, seq_length_to_skip=None,
            cluster_threshold=None, cluster_merge_after=None,
            cluster_method=None, cluster_fragment_length=None):
        """
        Args:
            genomes: list [g_1, g_2, g_m] of m groupings of genomes, where
                each g_i is a list of genome.Genomes belonging to group i.
                For example, a group may be a species and each g_i would be
                a list of the target genomes of species i.
            filters: an (ordered) list of filters, each of which should be
                an instance of a subclass of BaseFilter
            probe_length: generate candidate probes with this number of bp
            probe_stride: generate probes from each sequence separated by
                this number of bp
            allow_small_seqs: if set, allow sequences that are smaller than the
                probe length by creating candidate probes equal to the sequence;
                the value gives the minimum allowed probe (sequence) length
            seq_length_to_skip: if set, skip sequences whose length is <=
                the given value (i.e., do not design candidate probes for
                them)
            cluster_threshold: if set, cluster genomes and run filters (except
                the final ones) separately on each cluster, and merge probes;
                when cluster_method is 'simple', this value is the maximum
                distance at which to consider two sequences as adjacent when
                determining connected components; when cluster_method is
                'hierarchical', this value is the inter-cluster distance to
                merge clusters; for both, expressed in average nucleotide
                dissimilarity (1-ANI, where ANI is average nucleotide
                identity); higher results in fewer clusters
            cluster_merge_after: a filter in filters such that output probes
                from each cluster are merged just after running this filter,
                and all subsequent filters are run on the merged list; must be
                set if cluster_threshold is set
            cluster_method: method for performing cluster -- 'simple' for
                determining clusters based on connected components decided
                based on a distance threshold (relative less resource
                intensive); 'hierarchical' for agglomerative hierarchical
                clustering (more resource intensive); must be set if
                cluster_threshold is set
            cluster_fragment_length: if set, break genomes into fragments of
                this length and cluster these fragments rather than the whole
                sequences
        """
        self.genomes = genomes
        self.filters = filters
        self.probe_length = probe_length
        self.probe_stride = probe_stride
        self.allow_small_seqs = allow_small_seqs
        self.seq_length_to_skip = seq_length_to_skip
        self.cluster_threshold = cluster_threshold
        self.cluster_merge_after = cluster_merge_after
        self.cluster_method = cluster_method
        self.cluster_fragment_length = cluster_fragment_length

    def _cluster_genomes(self):
        """Cluster genomes by nucleotide similarity using MinHash signatures.

        This collapses all sequences, across both groups and genomes, into one
        list of sequences in one group. Therefore, genomes will not be grouped
        as specified in the input and sequences will not be grouped by genome.

        if self.cluster_fragment_length is specified, this additionally breaks
        genomes into fragments of that length, and clusters those fragments
        instead of the whole sequences.

        Returns:
            list of Genome objects, where each consists of a single sequence
        """
        if len(self.genomes) > 1:
            logger.warning(("There are >1 groups of genomes in the input, but "
                "clustering these will result in just one group; differential "
                "identification or other tasks that rely on group separation "
                "will no longer work"))

        # Create a dict of sequences across all groupings and genomes, for
        # input to clustering
        seqs = {}
        seq_idx = 0
        for genomes_from_group in self.genomes:
            for g in genomes_from_group:
                if self.cluster_fragment_length is not None:
                    # Break g into fragments
                    g_fragments = g.break_into_fragments(
                            self.cluster_fragment_length,
                            include_full_end=True)
                    g_seqs = g_fragments.seqs
                else:
                    # Do not break into fragments; use whole sequences
                    g_seqs = g.seqs
                for s in g_seqs:
                    if (self.seq_length_to_skip is not None and
                            len(s) <= self.seq_length_to_skip):
                        # Skip this sequence, which is too short
                        continue
                    seqs[seq_idx] = s
                    seq_idx += 1

        logger.info(("Clustering %d sequences using MinHash signatures, at an "
            "average nucleotide dissimilarity threshold of %f"),
                seq_idx, self.cluster_threshold)
        clusters = cluster.cluster_with_minhash_signatures(seqs,
                threshold=self.cluster_threshold,
                cluster_method=self.cluster_method)

        logger.info(("Found %d clusters with sizes: %s"), len(clusters),
                [len(clust) for clust in clusters])

        clustered_genomes = []
        for clust in clusters:
            genomes_in_clust = []
            for seq_idx in clust:
                seq = seqs[seq_idx]
                genomes_in_clust += [genome.Genome.from_one_seq(seq)]
            clustered_genomes += [genomes_in_clust]
        return clustered_genomes

    def _pass_through_filters(self, probes, genomes, filters):
        """Pass candidate probes through a list of filters.

        Args:
            probes: list [p_1, p_2, ..., p_m] of m groupings or clusters of
                genomes, where each p_i is a collection of probes for group i
            genomes: list [g_1, g_2, ..., g_m] of m groupings or clusters of
                genomes, where each g_i is a list of genome.Genomes belonging
                to group i. For example, a group may be a species and each g_i
                would be a list of the target genomes of species i.
            filters: an (ordered) list of filters, each of which should be
                an instance of a subclass of BaseFilter

        Returns:
            length-m list of subsets of probes, after passing through the
            filters
        """
        assert len(probes) == len(genomes)
        for f in filters:
            logger.info("Starting filter %s", f.__class__.__name__)
            probes = f.filter(probes, genomes, input_is_grouped=True)
        return probes

    def _pass_through_filters_ungrouped(self, probes, genomes, filters):
        """Pass candidate probes through a list of filters, where probes are
           not grouped.

        Args:
            probes: list of probes
            genomes: list [g_1, g_2, ..., g_m] of m groupings or clusters of
                genomes, where each g_i is a list of genome.Genomes belonging
                to group i. For example, a group may be a species and each g_i
                would be a list of the target genomes of species i.
            filters: an (ordered) list of filters, each of which should be
                an instance of a subclass of BaseFilter

        Returns:
            list of subsets of probes, after passing through the filters
        """
        for f in filters:
            logger.info("Starting filter %s", f.__class__.__name__)
            probes = f.filter(probes, genomes, input_is_grouped=False)
        return probes

    def _design_for_genomes(self, genomes, filters):
        """Design probes on a subset of genomes using a subset of filters.

        Generates the set of candidate probes and runs these through the
        filters.

        Args:
            genomes: list [g_1, g_2, ..., g_m] of m groupings or clusters of
                genomes, where each g_i is a list of genome.Genomes belonging
                to group i. For example, a group may be a species and each g_i
                would be a list of the target genomes of species i.
            filters: an (ordered) list of filters, each of which should be
                an instance of a subclass of BaseFilter

        Returns:
            tuple (candidate_probes, output_probes) where candidate_probes is a
            list of candidate probes and output_probes is the list of candidate
            probes after passing through the provided filters (each is grouped)
        """
        logger.info("Building candidate probes from target sequences")
        candidates = []
        for genomes_from_group in genomes:
            candidates_for_group = []
            for g in genomes_from_group:
                candidates_for_group += candidate_probes.\
                    make_candidate_probes_from_sequences(
                        g.seqs, probe_length=self.probe_length,
                        probe_stride=self.probe_stride,
                        allow_small_seqs=self.allow_small_seqs,
                        seq_length_to_skip=self.seq_length_to_skip)
            if len(candidates_for_group) == 0:
                # There are no candidate probes, possibly because all input
                # sequences in genomes were skipped
                logger.warning(("There are no candidate probes for a grouping of "
                    "genomes; it is possible that --small-seq-skip or "
                    "--small-seq-min are incompatible with the input sequence "
                    "lengths, especially if --cluster-and-design-separately is "
                    "set small."))
            candidates += [candidates_for_group]

        probes = self._pass_through_filters(candidates, genomes, filters)
        return (candidates, probes)

    def design(self):
        """Design probes using the provided filters.

        If desired, this first clusters sequences and then designs separately
        for each cluster. This stores the resulting probes in
        self.final_probes and all the candidate probes in
        self.candidate_probes.
        """
        if self.cluster_threshold is None:
            # Do not cluster sequences; just design on self.genomes as given
            candidates, probes = self._design_for_genomes(self.genomes,
                    self.filters)

            # Flatten each list; they are currently grouped
            self.candidate_probes = list(itertools.chain(*candidates))
            self.final_probes = list(set(itertools.chain(*probes)))
            return

        # Find filters before and after merging
        assert self.cluster_merge_after is not None
        assert self.cluster_merge_after in self.filters
        filter_merge_idx = self.filters.index(self.cluster_merge_after) + 1
        filters_before_merge = self.filters[:filter_merge_idx]
        filters_after_merge = self.filters[filter_merge_idx:]

        # Cluster genomes and design probes on each cluster
        clustered_genomes = self._cluster_genomes()
        candidates_by_cluster, probes_by_cluster = \
                self._design_for_genomes(clustered_genomes,
                        filters_before_merge)

        # Flatten the candidate probes, which are grouped
        self.candidate_probes = list(itertools.chain(*candidates_by_cluster))

        # Flatten the probes, which are grouped
        probes = list(set(itertools.chain(*probes_by_cluster)))

        # Run the remaining filters (filters_after_merge) on all the probes
        # after the merge
        probes = self._pass_through_filters_ungrouped(probes, clustered_genomes,
                filters_after_merge)

        self.final_probes = probes
