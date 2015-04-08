"""Dataset with Marburg sequences.

A dataset with 78 Marburg samples.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

import sys
from hybseldesign.datasets import GenomesDatasetSingleChrom


ds = GenomesDatasetSingleChrom(__name__, __file__)
ds.set_fasta_path("data/marburg.fasta", relative=True)
sys.modules[__name__] = ds

