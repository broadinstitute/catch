"""Dataset with HIV-2 sequences in which the LTRs are removed.

A dataset with 30 HIV-2 genomes where the LTR is stripped from
both the 5' and 3' end of each genome.
"""

import sys

from hybseldesign.datasets import GenomesDatasetSingleChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'

ds = GenomesDatasetSingleChrom(__name__, __file__, __spec__)
ds.set_fasta_path("data/hiv2_without_ltr.fasta", relative=True)
sys.modules[__name__] = ds
