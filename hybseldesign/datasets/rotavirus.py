"""Dataset with Rotavirus sequences.

A dataset with 253 Rotavirus segments. (The labels for some
segments are unclear -- i.e., seem to be 'null'.)
"""

import sys

from hybseldesign.datasets import GenomesDatasetSingleChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'

ds = GenomesDatasetSingleChrom(__name__, __file__)
ds.set_fasta_path("data/rotavirus.fasta", relative=True)
sys.modules[__name__] = ds
