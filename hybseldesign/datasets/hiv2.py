"""Dataset with HIV-2 sequences.

A dataset with 30 HIV-2 genomes.
"""

import sys

from hybseldesign.datasets import GenomesDatasetSingleChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'

ds = GenomesDatasetSingleChrom(__name__, __file__, __spec__)
ds.set_fasta_path("data/hiv2.fasta", relative=True)
sys.modules[__name__] = ds
