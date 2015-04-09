"""Dataset with Ebola sequences from outside Zaire.

A dataset with 53 Ebola samples taken from places that do not
include Zaire. These samples do not include those in the ebola2014
dataset.
"""

import sys

from hybseldesign.datasets import GenomesDatasetSingleChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'

ds = GenomesDatasetSingleChrom(__name__, __file__)
ds.set_fasta_path("data/ebola_nonzaire.fasta", relative=True)
sys.modules[__name__] = ds
