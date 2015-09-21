"""Dataset with GB virus C sequences.

A dataset with 17 GB virus C (GBV-C, formerly known as hepatitis
G virus) genomes.
"""

import sys

from hybseldesign.datasets import GenomesDatasetSingleChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'

ds = GenomesDatasetSingleChrom(__name__, __file__, __spec__)
ds.add_fasta_path("data/gbv_c.fasta", relative=True)
sys.modules[__name__] = ds
