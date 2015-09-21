"""Dataset with all deer tick borne pathogens.

A dataset with the full genomes of all deer tick borne pathogens
(41 in total). These are mostly bacterial genomes with one genome
per bacteria (i.e., a sample of n=1 for each pathogen).
"""

import sys

from hybseldesign.datasets import GenomesDatasetSingleChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'

ds = GenomesDatasetSingleChrom(__name__, __file__, __spec__)
ds.add_fasta_path("data/all_deer_tick_borne.fasta", relative=True)
sys.modules[__name__] = ds
