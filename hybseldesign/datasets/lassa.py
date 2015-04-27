"""Dataset with Lassa virus sequences.

A dataset with 107 Lassa virus samples. Each sequence corresponds
to either a small (S) or large (L) segment of the genome; for a
particular sample, the sequences for its two segments are adjacent
in the FASTA file.
"""

import sys

from hybseldesign.datasets import GenomesDatasetSingleChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'

ds = GenomesDatasetSingleChrom(__name__, __file__)
ds.set_fasta_path("data/lassa.fasta", relative=True)
sys.modules[__name__] = ds
