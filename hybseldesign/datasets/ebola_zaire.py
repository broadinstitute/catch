"""Dataset with Zaire Ebola sequences.

A dataset with 140 Ebola Zaire genomes. These do not include
samples from the ebola2014 dataset (i.e., from the 2014
West Africa Ebola outbreak).
"""

import sys

from hybseldesign.datasets import GenomesDatasetSingleChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'

ds = GenomesDatasetSingleChrom(__name__, __file__)
ds.set_fasta_path("data/ebola_zaire.fasta", relative=True)
sys.modules[__name__] = ds
