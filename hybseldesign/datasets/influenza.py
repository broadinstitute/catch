"""Dataset with Influenza sequences.

A dataset with 8985 Influenza segments (both Influenza A and
Influenza B).
"""

import sys

from hybseldesign.datasets import GenomesDatasetSingleChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'

ds = GenomesDatasetSingleChrom(__name__, __file__)
ds.set_fasta_path("data/influenza.fasta", relative=True)
sys.modules[__name__] = ds
