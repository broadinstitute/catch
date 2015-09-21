"""Dataset with Hepatitis-A sequences.

A dataset with 34 Hepatitis-A genomes.
"""

import sys

from hybseldesign.datasets import GenomesDatasetSingleChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'

ds = GenomesDatasetSingleChrom(__name__, __file__, __spec__)
ds.add_fasta_path("data/hepatitis_a.fasta", relative=True)
sys.modules[__name__] = ds
