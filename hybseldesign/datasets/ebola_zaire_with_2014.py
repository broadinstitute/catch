"""Dataset with Zaire Ebola sequences, including from the 2014 outbreak.

A dataset with 239 Ebola Zaire genomes. These consist of samples
from the ebola_zaire dataset and from the ebola2014 dataset
(i.e., unlike the ebola_zaire dataset, these include samples
from the 2014 West Africa Ebola outbreak).
"""

import sys

from hybseldesign.datasets import GenomesDatasetSingleChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'

ds = GenomesDatasetSingleChrom(__name__, __file__, __spec__)
ds.add_fasta_path("data/ebola_zaire.fasta", relative=True)
ds.add_fasta_path("data/ebola2014.fasta", relative=True)
sys.modules[__name__] = ds
