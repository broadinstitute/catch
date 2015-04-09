"""Dataset with Powassan sequences.

A dataset with 12 Powassan samples. This is a subset of the
genomes in the all_deer_tick_borne dataset.
"""

import sys

from hybseldesign.datasets import GenomesDatasetSingleChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'

ds = GenomesDatasetSingleChrom(__name__, __file__)
ds.set_fasta_path("data/powassan.fasta", relative=True)
sys.modules[__name__] = ds

