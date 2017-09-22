"""Dataset with 'Rodent Torque teno virus( [0-9]+)?' sequences.

A dataset with 2 'Rodent Torque teno virus( [0-9]+)?' genomes.

THIS PYTHON FILE WAS GENERATED BY A COMPUTER PROGRAM! DO NOT EDIT!
"""

import sys

from hybseldesign.datasets import GenomesDatasetSingleChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'


ds = GenomesDatasetSingleChrom(__name__, __file__, __spec__)
ds.add_fasta_path("data/rodent_torque_teno.fasta", relative=True)
sys.modules[__name__] = ds
