"""Dataset with Ebola sequences from 2014 outbreak.

A dataset with 99 genomes of Ebola taken from the 2014 outbreak.
"""

import sys

from hybseldesign.datasets import GenomesDatasetSingleChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'

ds = GenomesDatasetSingleChrom(__name__, __file__)
ds.set_fasta_path("data/ebola2014.fasta", relative=True)
sys.modules[__name__] = ds
