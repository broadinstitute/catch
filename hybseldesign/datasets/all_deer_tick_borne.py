"""A dataset with the full genomes of all deer tick borne pathogens
(41 in total). These are mostly bacterial genomes with one genome
per bacteria (i.e., a sample of n=1 for each pathogen).
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

import sys
from hybseldesign.datasets import GenomesDatasetSingleChrom


ds = GenomesDatasetSingleChrom(__name__, __file__)
ds.set_fasta_path("data/all_deer_tick_borne.fasta", relative=True)
sys.modules[__name__] = ds
