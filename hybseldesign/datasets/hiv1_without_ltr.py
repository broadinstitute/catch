"""Dataset with 'Human immunodeficiency virus 1' sequences without LTRs.

A dataset with 1962 'Human immunodeficiency virus 1' genomes where the
LTR is stripped from both the 5' and 3' ends of each genome. These
sequences are based on the sequences in the 'hiv1' dataset.
"""

import sys

from hybseldesign.datasets import GenomesDatasetSingleChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'


ds = GenomesDatasetSingleChrom(__name__, __file__, __spec__)
ds.add_fasta_path("data/hiv1_without_ltr.fasta", relative=True)
sys.modules[__name__] = ds
