"""Dataset with Measles sequences.

A dataset with 53 Measles genomes.
"""

import sys

from hybseldesign.datasets import GenomesDatasetSingleChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'

ds = GenomesDatasetSingleChrom(__name__, __file__, __spec__)
ds.set_fasta_path("data/measles.fasta", relative=True)
sys.modules[__name__] = ds
