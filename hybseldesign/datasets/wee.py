"""Dataset with 'Western equine encephalitis virus' sequences.

A dataset with 25 'Western equine encephalitis virus' genomes.

THIS PYTHON FILE WAS GENERATED BY A COMPUTER PROGRAM! DO NOT EDIT!
"""

import sys

from hybseldesign.datasets import GenomesDatasetSingleChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'


ds = GenomesDatasetSingleChrom(__name__, __file__, __spec__)
ds.add_fasta_path("data/wee.fasta", relative=True)
sys.modules[__name__] = ds
