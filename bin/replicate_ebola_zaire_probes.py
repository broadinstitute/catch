"""Script that runs hybseldesign to replicate the probes generated
on the Ebola Zaire dataset by the prior version of code (Matlab
and RemoveSimilarSequences.jar).
"""

# Author: Hayden Metsky <hayden@mit.edu>

import argparse

from hybseldesign.datasets import ebola_zaire
from hybseldesign.utils import seq_io
from hybseldesign.filter import probe_designer
from hybseldesign.filter import reverse_complement_filter
from hybseldesign.filter import duplicate_filter
from hybseldesign.filter import naive_redundant_filter


def main(args):
  # Read the Ebola Zaire FASTA sequences
  seqs = seq_io.read_fasta(ebola_zaire.fasta_path()).values()

  # Setup the filters needed for replication
  # The filters needed for replication are, in order:
  #  1) Reverse complement (rc) -- add the reverse complement of
  #     each probe as a candidate
  rc = reverse_complement_filter.ReverseComplementFilter()
  #  2) Duplicate filter (df) -- condense all candidate probes that
  #     are identical down to one; this is not necessary for
  #     correctness, as the naive redundant filter achieves the same
  #     task implicitly, but it does significantly lower runtime by
  #     decreasing the input size to the naive redundant filter
  df = duplicate_filter.DuplicateFilter()
  #  3) Naive redundant filter (nf) -- execute a greedy algorithm to
  #     condense 'similar' probes down to one
  nf = naive_redundant_filter.NaiveRedundantFilter(
        naive_redundant_filter.redundant_shift_and_mismatch_count(
          shift=args.shift, mismatch_thres=args.mismatch_thres))
  filters = [ rc, df, nf ]

  # Design the probes
  # (including the replicate_first_version argument to execute
  #  bugs that were present in the code's prior version when
  #  generating the set of candidate probes)
  pb = probe_designer.ProbeDesigner(seqs, filters,
          replicate_first_version=True)
  pb.design()

  # Print the arguments and the number of final probes
  # (The final probes are stored in pb.final_probes if their
  #  sequences are desired)
  print args.shift, args.mismatch_thres, len(pb.final_probes)


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("-s", "--shift", required=True, type=int,
      help=("Shift one probe from -shift to +shift bp relative to the "
            "other"))
  parser.add_argument("-m", "--mismatch_thres", required=True, type=int,
      help=("Deem one probe redundant to another if, as one is "
            "shifted relative to the other, the minimum number of "
            "mismatches between them is <= 'mismatch_thres'"))
  args = parser.parse_args()

  main(args)
