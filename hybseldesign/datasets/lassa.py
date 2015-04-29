"""Dataset with Lassa virus sequences.

A dataset with 214 Lassa virus segments. Each segment is either
a small (S) or large (L) segment of a genome. Most are paired
with another, but some are unmatched.
"""

import sys

from hybseldesign.datasets import GenomesDatasetSingleChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'

ds = GenomesDatasetSingleChrom(__name__, __file__)
ds.set_fasta_path("data/lassa.fasta", relative=True)
sys.modules[__name__] = ds
