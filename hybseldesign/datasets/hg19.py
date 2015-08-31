"""Dataset of the human genome 19.
"""

import sys

from hybseldesign.datasets import GenomesDatasetMultiChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'


def seq_header_to_chr(header):
    return header

chrs = [str(x) for x in range(1, 23)] + ["X", "Y", "MT"]
ds = GenomesDatasetMultiChrom(__name__, __file__, __spec__,
                              chrs, seq_header_to_chr)
sys.modules[__name__] = ds
