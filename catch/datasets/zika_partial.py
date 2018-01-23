"""Dataset with 'Zika virus' partial sequences.

A dataset with 221 'Zika virus' partial sequences. These consist mostly
of UTRs and partial coding sequences (i.e., particular proteins). This
dataset is distinct from the 'zika' dataset, which only contains full
genomes (or complete coding sequences). The sequences in the 'zika'
dataset are based on complete genomes in NCBI's viral accession list,
whereas the sequences in this dataset can be found with a search on
GenBank for:
    ("Zika virus"[porgn:__txid64320]) NOT complete

These partial sequences are kept separate from full genomes because
this software package is meant to use full genomes as input. Designing
probes with many partial sequences of the same region (e.g., the
same protein or the 3' UTR) may lead to more probes than are
necessary. The reason is that many of the partial sequences may be
slightly offset from each other (i.e., unlike in complete genomes,
their endpoints are not in alignment) and separate probes will be
needed to capture each of the ends of each of these partial sequences.
In actual hybrid selection, these probes are extraneous because it
is performed on complete genomes and so probes that might only
hybridize to a piece of the end of a partial sequence will actually
capture the desired fragment.
"""

import sys

from catch.datasets import GenomesDatasetSingleChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'


ds = GenomesDatasetSingleChrom(__name__, __file__, __spec__)
ds.add_fasta_path("data/zika_partial.fasta", relative=True)
sys.modules[__name__] = ds
