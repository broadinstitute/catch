"""Script that runs hybseldesign to design probes from the Ebola
2014 dataset using a filter-based approach.
"""

# Author: Hayden Metsky <hayden@mit.edu>

from hybseldesign.datasets import ebola2014
from hybseldesign.utils import seq_io
from hybseldesign.filter import probe_designer
from hybseldesign.filter import reverse_complement_filter
from hybseldesign.filter import duplicate_filter


# For now, simply:
# 1) add reverse complements to the set of candidate probes
# 2) filter duplicates from the set of candidate probes
seqs = seq_io.read_fasta(ebola2014.fasta_path()).values()
rc = reverse_complement_filter.ReverseComplementFilter()
df = duplicate_filter.DuplicateFilter()
pb = probe_designer.ProbeDesigner(seqs, [rc, df])
pb.design()

print "Number of candidate probes:", len(pb.candidate_probes)
print "Number of filtered probes:", len(pb.final_probes)
