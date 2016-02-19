"""Dataset with 'Rotavirus G' sequences.

A dataset with 19 'Rotavirus G' genomes.

THIS PYTHON FILE WAS GENERATED BY A COMPUTER PROGRAM! DO NOT EDIT!
"""

import sys

from hybseldesign.datasets import GenomesDatasetSingleChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'


ds = GenomesDatasetSingleChrom(__name__, __file__, __spec__)
ds.add_fasta_path("data/rotavirus_g.fasta", relative=True)
sys.modules[__name__] = ds