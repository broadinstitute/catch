"""Dataset with 'Avian metapneumovirus' sequences.

A dataset with 11 'Avian metapneumovirus' genomes.

THIS PYTHON FILE WAS GENERATED BY A COMPUTER PROGRAM! DO NOT EDIT!
"""

import sys

from hybseldesign.datasets import GenomesDatasetSingleChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'


ds = GenomesDatasetSingleChrom(__name__, __file__, __spec__)
ds.add_fasta_path("data/ampv.fasta", relative=True)
sys.modules[__name__] = ds
