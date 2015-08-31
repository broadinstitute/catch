"""Dataset with HIV-1 sequences.

A dataset with 1749 HIV-1 genomes.
"""

import sys

from hybseldesign.datasets import GenomesDatasetSingleChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'

ds = GenomesDatasetSingleChrom(__name__, __file__, __spec__)
ds.set_fasta_path("data/hiv1.fasta", relative=True)
sys.modules[__name__] = ds
