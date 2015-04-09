"""Dataset of the human genome 19.
"""

import sys

from hybseldesign.datasets import GenomesDatasetMultiChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'

chrs = [str(x) for x in range(1, 23)] + ["X", "Y", "MT"]
ds = GenomesDatasetMultiChrom(__name__, __file__, chrs)
sys.modules[__name__] = ds
