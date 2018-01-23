"""Chooses only probes present in a given FASTA file.

This acts as a filter on the probes by returning a list of
probes from the given input probes. This keeps all probes
that are equal to a sequence in a FASTA file, thereby removing
all probes that are not a sequence in the file. It preserves
the order of the sequences in the FASTA file; i.e., the input
probes are re-arranged to match their order in the file.
"""

from collections import OrderedDict

from catch.filter.base_filter import BaseFilter
from catch.utils import seq_io

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class FastaFilter(BaseFilter):
    """Filter that selects only probes equal to a sequence in a FASTA file.
    """

    def __init__(self, fasta_path, skip_reverse_complements=False):
        """
        Args:
            fasta_path: path to FASTA file whose sequences should be used
                in this filter
            skip_reverse_complements: does not read any sequences in the
                FASTA file that contain "reverse complement" in their
                header
        """
        self.fasta_path = fasta_path
        self.skip_reverse_complements = skip_reverse_complements

    def _filter(self, input):
        """Return a subset of the input probes.
        """
        # Read the FASTA file
        fasta = seq_io.read_fasta(self.fasta_path)

        # Construct a set of the sequences from the file
        seqs_to_keep = {}
        for i, (header, seq) in enumerate(fasta.items()):
            if self.skip_reverse_complements:
                if "reverse complement" not in header:
                    seqs_to_keep[seq] = i
            else:
                seqs_to_keep[seq] = i

        # Construct a list of tuples of the form
        # (probe's position in fasta file, probe)
        filtered = []
        for probe in input:
            if probe.seq_str in seqs_to_keep:
                order_in_file = seqs_to_keep[probe.seq_str]
                filtered += [(order_in_file, probe)]

        # Sort the filtered probes by their order
        filtered.sort()

        # Remove ordering information from the filtered list
        for i, (order_in_file, probe) in enumerate(filtered):
            filtered[i] = probe

        return filtered
