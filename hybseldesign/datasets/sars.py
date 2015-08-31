"""Dataset with SARS sequences.

A dataset with 143 SARS genomes.
"""

import sys

from hybseldesign.datasets import GenomesDatasetSingleChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'

ds = GenomesDatasetSingleChrom(__name__, __file__, __spec__)
ds.set_fasta_path("data/sars.fasta", relative=True)
sys.modules[__name__] = ds
