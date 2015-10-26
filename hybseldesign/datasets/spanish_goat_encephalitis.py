"""Dataset with 'Spanish goat encephalitis virus' sequences.

A dataset with 1 'Spanish goat encephalitis virus' genomes.

THIS PYTHON FILE WAS GENERATED BY A COMPUTER PROGRAM! DO NOT EDIT!
"""

import sys

from hybseldesign.datasets import GenomesDatasetSingleChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'


ds = GenomesDatasetSingleChrom(__name__, __file__, __spec__)
ds.add_fasta_path("data/spanish_goat_encephalitis.fasta", relative=True)
sys.modules[__name__] = ds
