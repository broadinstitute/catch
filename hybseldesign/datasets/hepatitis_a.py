"""Dataset with Hepatitis-A sequences.

A dataset with 34 Hepatitis-A samples.
"""

import sys

from hybseldesign.datasets import GenomesDatasetSingleChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'

ds = GenomesDatasetSingleChrom(__name__, __file__)
ds.set_fasta_path("data/hepatitis_a.fasta", relative=True)
sys.modules[__name__] = ds
