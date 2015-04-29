"""Dataset with HIV sequences.

A dataset with 1779 HIV genomes.
"""

import sys

from hybseldesign.datasets import GenomesDatasetSingleChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'

ds = GenomesDatasetSingleChrom(__name__, __file__)
ds.set_fasta_path("data/hiv.fasta", relative=True)
sys.modules[__name__] = ds
