"""Dataset with HIV-1 sequences in which the LTRs are removed.

A dataset with 1749 HIV-1 genomes where the LTR is stripped from
both the 5' and 3' end of each genome.
"""

import sys

from hybseldesign.datasets import GenomesDatasetSingleChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'

ds = GenomesDatasetSingleChrom(__name__, __file__, __spec__)
ds.add_fasta_path("data/hiv1_without_ltr.fasta", relative=True)
sys.modules[__name__] = ds
