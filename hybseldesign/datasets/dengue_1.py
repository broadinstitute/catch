"""Dataset with 'Dengue virus 1' sequences.

A dataset with 1599 'Dengue virus 1' genomes.

Note that the sequences in this dataset are a subset of those in the
'dengue' dataset.

THIS PYTHON FILE WAS GENERATED BY A COMPUTER PROGRAM! DO NOT EDIT!
"""

import sys

from hybseldesign.datasets import GenomesDatasetSingleChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'


ds = GenomesDatasetSingleChrom(__name__, __file__, __spec__)
ds.add_fasta_path("data/dengue_1.fasta", relative=True)
sys.modules[__name__] = ds
