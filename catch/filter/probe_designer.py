"""Designs probes with a filtering approach.
"""

import logging

from catch.filter import candidate_probes

__author__ = 'Hayden Metsky <hayden@mit.edu>'

logger = logging.getLogger(__name__)


class ProbeDesigner:
    """Generates a set of candidate probes and filters them.

    There is an (ordered) list of one or more provided filters.
    """

    def __init__(self, genomes, filters, probe_length,
            probe_stride, allow_small_seqs=None, seq_length_to_skip=None):
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
        """
        self.genomes = genomes
        self.filters = filters
        self.probe_length = probe_length
        self.probe_stride = probe_stride
        self.allow_small_seqs = allow_small_seqs
        self.seq_length_to_skip = seq_length_to_skip

    def design(self):
        """Design probes using the provided filters.

        Generates the set of candidate probes and runs these through the
        filters. Stores the candidate probes in self.candidate_probes and
        the probes processed by the filters in self.final_probes.
        """

        logger.info("Building candidate probes from target sequences")
        self.candidate_probes = []
        for genomes_from_group in self.genomes:
            for g in genomes_from_group:
                self.candidate_probes += candidate_probes.\
                    make_candidate_probes_from_sequences(
                        g.seqs, probe_length=self.probe_length,
                        probe_stride=self.probe_stride,
                        allow_small_seqs=self.allow_small_seqs,
                        seq_length_to_skip=self.seq_length_to_skip)

        probes = self.candidate_probes
        for f in self.filters:
            logger.info("Starting filter %s", f.__class__.__name__)
            f.target_genomes = self.genomes
            probes = f.filter(probes)
        self.final_probes = probes
