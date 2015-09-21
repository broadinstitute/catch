"""Dataset with Marburg sequences.

A dataset with 78 Marburg genomes.
"""

import sys

from hybseldesign.datasets import GenomesDatasetSingleChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'

ds = GenomesDatasetSingleChrom(__name__, __file__, __spec__)
ds.add_fasta_path("data/marburg.fasta", relative=True)
sys.modules[__name__] = ds
