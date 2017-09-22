"""Dataset with 'Human mastadenovirus B' sequences.

A dataset with 111 'Human mastadenovirus B' genomes.

THIS PYTHON FILE WAS GENERATED BY A COMPUTER PROGRAM! DO NOT EDIT!
"""

import sys

from hybseldesign.datasets import GenomesDatasetSingleChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'


ds = GenomesDatasetSingleChrom(__name__, __file__, __spec__)
ds.add_fasta_path("data/mastadenovirus_b.fasta", relative=True)
sys.modules[__name__] = ds
