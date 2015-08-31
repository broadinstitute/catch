"""Dataset with Rift Valley fever virus sequences.

A dataset with 158 Rift Valley fever virus genomes. Many
of these are complete (i.e., contain the three segments S,
L, and M), but some are incomplete (e.g., contain just one
segment).
"""

from os.path import dirname
from os.path import join
from os import listdir
import sys

from hybseldesign.datasets import GenomesDatasetMultiChrom

__author__ = 'Hayden Metsky <hayden@mit.edu>'


chrs = ["segment_S", "segment_L", "segment_M"]

def seq_header_to_chr(header):
    import re
    c = re.compile(r'segment ([SLM])')
    m = c.search(header)
    if not m:
        raise ValueError("Unknown segment in header %s" % header)
    seg = m.group(1)
    valid_segs = ["S", "L", "M"]
    if seg not in valid_segs:
        raise ValueError("Unknown segment %s" % seg)
    return "segment_" + seg


ds = GenomesDatasetMultiChrom(__name__, __file__, __spec__,
                              chrs, seq_header_to_chr)

for f in listdir(join(dirname(__file__), "data/rift_valley_fever/")):
    ds.add_fasta_path("data/rift_valley_fever/" + f, relative=True)

sys.modules[__name__] = ds
