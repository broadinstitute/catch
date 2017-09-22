"""Dataset with 'Yellow fever virus' sequences.

A dataset with 42 'Yellow fever virus' genomes.

THIS PYTHON FILE WAS GENERATED BY A COMPUTER PROGRAM! DO NOT EDIT!
"""

import sys

from hybseldesign.datasets import GenomesDatasetSingleChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'


ds = GenomesDatasetSingleChrom(__name__, __file__, __spec__)
ds.add_fasta_path("data/yellow_fever.fasta", relative=True)
sys.modules[__name__] = ds
