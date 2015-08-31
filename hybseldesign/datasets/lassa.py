"""Dataset with Lassa virus sequences.

A dataset with 222 Lassa virus genomes. Each segment is either
a small (S) or large (L) segment of a genome. Most genomes have
both segments, but some are incomplete and only have one segment.
The total number of segments across the 222 genomes is 413.
"""

from os.path import dirname
from os.path import join
from os import listdir
import sys

from hybseldesign.datasets import GenomesDatasetMultiChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'


chrs = ["segment_S", "segment_L"]

def seq_header_to_chr(header):
    import re
    c = re.compile(r'segment ([SL])')
    m = c.search(header)
    if not m:
        raise ValueError("Unknown segment in header %s" % header)
    seg = m.group(1)
    valid_segs = ["S", "L"]
    if seg not in valid_segs:
        raise ValueError("Unknown segment %s" % seg)
    return "segment_" + seg


ds = GenomesDatasetMultiChrom(__name__, __file__, __spec__,
                              chrs, seq_header_to_chr)

for f in listdir(join(dirname(__file__), "data/lassa/")):
    ds.add_fasta_path("data/lassa/" + f, relative=True)

sys.modules[__name__] = ds
