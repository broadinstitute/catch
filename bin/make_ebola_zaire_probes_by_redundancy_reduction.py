#!/usr/bin/env python
"""Script that runs hybseldesign to design probes for the Ebola
Zaire dataset by using either a naive redundant filter or a
dominating set filter to remove redundancies in candidate probes
based on the longest common substring between two probes.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

import argparse

from hybseldesign.datasets import ebola_zaire
from hybseldesign.utils import seq_io, version
from hybseldesign.filter import probe_designer
from hybseldesign.filter import reverse_complement_filter
from hybseldesign.filter import duplicate_filter
from hybseldesign.filter import naive_redundant_filter
from hybseldesign.filter import dominating_set_filter


def main(args):
  # Read the Ebola Zaire FASTA sequences
  seqs = seq_io.read_fasta(ebola_zaire.fasta_path()).values()

  # Setup the filters
  # The filters we use are, in order:
  #  1) Duplicate filter (df) -- condense all candidate probes that
  #     are identical down to one; this is not necessary for
  #     correctness, as the naive redundant filter achieves the same
  #     task implicitly, but it does significantly lower runtime by
  #     decreasing the input size to the naive redundant filter
  df = duplicate_filter.DuplicateFilter()
  if args.main_filter == "ds":
    # 2) Dominating set filter (ds) -- treat the problem as an
    #    instance of the dominating set problem
    mf = dominating_set_filter.DominatingSetFilter(
          naive_redundant_filter.redundant_longest_common_substring(
            mismatches=args.mismatches, lcf_thres=args.lcf_thres))
  elif args.main_filter == "nr":
    #  2) Naive redundant filter (nr) -- execute a greedy algorithm
    #     to condense 'similar' probes down to one
    mf = naive_redundant_filter.NaiveRedundantFilter(
          naive_redundant_filter.redundant_longest_common_substring(
            mismatches=args.mismatches, lcf_thres=args.lcf_thres))
  else:
    raise ValueError("Unknown main filter " + args.main_filter)
  #  3) Reverse complement (rc) -- add the reverse complement of
  #     each probe that remains
  rc = reverse_complement_filter.ReverseComplementFilter()
  filters = [ df, mf, rc ]

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
  print args.mismatches, args.lcf_thres, len(pb.final_probes)


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("-m", "--mismatches", required=True, type=int,
      help=("Determine redundancy based on the longest common "
            "substring with at most 'mismatches' mismatches "
            "between two probes"))
  parser.add_argument("-l", "--lcf_thres", required=True, type=int,
      help=("Deem two probes to be redundant if the longest "
            "common substring with at most 'mismatches' mismatches "
            "between them has a length that is >= 'lcf_thres' bp"))
  parser.add_argument("-f", "--main_filter", type=str, default="ds",
      help=("The primary filter to use: 'ds' (dominating set filter) "
            "[default] or 'nr' (naive redundant filter)"))
  parser.add_argument('--version', '-V', action='version',
      version=version.get_version())
  args = parser.parse_args()

  main(args)
