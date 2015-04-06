"""A dataset with 12 Powassan samples. This is a subset of the
genomes in the all_deer_tick_borne dataset.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

import sys
from hybseldesign.datasets import GenomesDatasetSingleChrom


ds = GenomesDatasetSingleChrom(__name__, __file__)
ds.set_fasta_path("data/powassan.fasta", relative=True)
sys.modules[__name__] = ds

