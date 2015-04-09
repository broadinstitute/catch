"""Dataset with Ebola sequences from within Zaire.

A dataset with 140 Ebola samples taken from Zaire.
"""

import sys

from hybseldesign.datasets import GenomesDatasetSingleChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'

ds = GenomesDatasetSingleChrom(__name__, __file__)
ds.set_fasta_path("data/ebola_zaire.fasta", relative=True)
sys.modules[__name__] = ds
