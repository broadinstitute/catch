"""Dataset with Yellow fever virus sequences.

A dataset with 55 Yellow fever virus genomes.
"""

import sys

from hybseldesign.datasets import GenomesDatasetSingleChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'

ds = GenomesDatasetSingleChrom(__name__, __file__)
ds.set_fasta_path("data/yellow_fever.fasta", relative=True)
sys.modules[__name__] = ds
