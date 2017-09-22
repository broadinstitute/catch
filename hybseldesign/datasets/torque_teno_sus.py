"""Dataset with 'Torque teno sus virus( [a-z0-9]+)?' sequences.

A dataset with 72 'Torque teno sus virus( [a-z0-9]+)?' genomes.

THIS PYTHON FILE WAS GENERATED BY A COMPUTER PROGRAM! DO NOT EDIT!
"""

import sys

from hybseldesign.datasets import GenomesDatasetSingleChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'


ds = GenomesDatasetSingleChrom(__name__, __file__, __spec__)
ds.add_fasta_path("data/torque_teno_sus.fasta", relative=True)
sys.modules[__name__] = ds
