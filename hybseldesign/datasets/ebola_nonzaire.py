"""Dataset with Ebola sequences from outside Zaire.

A dataset with 53 Ebola genomes taken from places that do not
include Zaire. These are Ebola viruses that are not Zaire
ebolavirus.
"""

import sys

from hybseldesign.datasets import GenomesDatasetSingleChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'

ds = GenomesDatasetSingleChrom(__name__, __file__)
ds.set_fasta_path("data/ebola_nonzaire.fasta", relative=True)
sys.modules[__name__] = ds
