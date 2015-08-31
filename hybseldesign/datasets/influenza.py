"""Dataset with Influenza sequences.

A dataset with 1125 Influenza genomes, with a total of 8983 segments
(both Influenza A and Influenza B); some are incomplete.
"""

from os.path import dirname
from os.path import join
from os import listdir
import sys

from hybseldesign.datasets import GenomesDatasetMultiChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'


chrs = ["segment_" + str(i) for i in range(1, 9)]

def seq_header_to_chr(header):
    import re
    c = re.compile(r'Segment:([0-9]+)')
    m = c.search(header)
    if not m:
        raise ValueError("Unknown segment in header %s" % header)
    seg_num = m.group(1)
    valid_seg_nums = [str(x) for x in range(1, 9)]
    if seg_num not in valid_seg_nums:
        raise ValueError("Unknown segment %s" % seg_num)
    return "segment_" + seg_num


ds = GenomesDatasetMultiChrom(__name__, __file__, __spec__,
                              chrs, seq_header_to_chr)

for f in listdir(join(dirname(__file__), "data/influenza/")):
    ds.add_fasta_path("data/influenza/" + f, relative=True)

sys.modules[__name__] = ds
