"""Dataset with MERS virus sequences.

A dataset with 204 MERS virus genomes. 152 of these come from
virological.org (part of MERS-CoV_22_15Apr15.fas.txt); some
of these 152 sequences contain gap characters ('-').
"""

import sys

from hybseldesign.datasets import GenomesDatasetSingleChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'

ds = GenomesDatasetSingleChrom(__name__, __file__)
ds.set_fasta_path("data/mers.fasta", relative=True)
sys.modules[__name__] = ds
