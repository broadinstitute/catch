"""Dataset with 'Zika virus' partial sequences.

A dataset with 221 'Zika virus' partial sequences. These consist mostly
of UTRs and partial coding sequences (i.e., particular proteins). This
dataset is distinct from the 'zika' dataset, which only contains full
genomes (or complete coding sequences). The sequences in the 'zika'
dataset can be found on GenBank with a search for:
    "Zika virus"[porgn:__txid64320] AND complete
whereas the sequences in this dataset can be found with a search for:
    ("Zika virus"[porgn:__txid64320]) NOT complete

**why we keep them separate / why probe design is hard on partial
sequences**
"""

import sys

from hybseldesign.datasets import GenomesDatasetSingleChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'


ds = GenomesDatasetSingleChrom(__name__, __file__, __spec__)
ds.add_fasta_path("data/zika_partial.fasta", relative=True)
sys.modules[__name__] = ds
