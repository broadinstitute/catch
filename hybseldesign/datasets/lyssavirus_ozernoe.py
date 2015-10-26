"""Dataset with 'Lyssavirus Ozernoe' sequences.

A dataset with 1 'Lyssavirus Ozernoe' genomes.

THIS PYTHON FILE WAS GENERATED BY A COMPUTER PROGRAM! DO NOT EDIT!
"""

import sys

from hybseldesign.datasets import GenomesDatasetSingleChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'


ds = GenomesDatasetSingleChrom(__name__, __file__, __spec__)
ds.add_fasta_path("data/lyssavirus_ozernoe.fasta", relative=True)
sys.modules[__name__] = ds
