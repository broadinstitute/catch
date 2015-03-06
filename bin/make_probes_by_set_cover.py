#!/usr/bin/env python
"""Script that runs hybseldesign to design probes for a given
dataset by treating the problem as an instance of set cover.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

import argparse

from hybseldesign.datasets import ebola_zaire
from hybseldesign.datasets import ebola2014
from hybseldesign.datasets import hg19
from hybseldesign.datasets import marburg
from hybseldesign.utils import seq_io, version
from hybseldesign.filter import probe_designer
from hybseldesign.filter import reverse_complement_filter
from hybseldesign.filter import duplicate_filter
from hybseldesign.filter import set_cover_filter

DATASETS = { "ebola_zaire": ebola_zaire,
             "ebola2014": ebola2014 }

# Set the path to the hg19 fasta file (assuming this is a Broad
# machine)
hg19.set_fasta_path("/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta")


def main(args):
  # Read the FASTA sequences
  if args.dataset not in DATASETS:
    raise ValueError("Unknown dataset %s" % args.dataset)
  fasta_path = DATASETS[args.dataset].fasta_path()
  seqs = seq_io.read_fasta(fasta_path).values()

  if args.limit_target_genomes:
    seqs = seqs[:args.limit_target_genomes]

  blacklisted_genomes = []
  if args.blacklist_hg19:
    blacklisted_genomes += [ hg19.fasta_path() ]
  if args.blacklist_marburg:
    blacklisted_genomes += [ marburg.fasta_path() ]

  # Setup the filters
  # The filters we use are, in order:
  #  1) Duplicate filter (df) -- condense all candidate probes that
  #     are identical down to one; this is not necessary for
  #     correctness, as the set cover filter achieves the same
  #     task implicitly, but it does significantly lower runtime by
  #     decreasing the input size to the set cover filter
  df = duplicate_filter.DuplicateFilter()
  #  2) Set cover filter (scf) -- solve the problem by treating it
  #     as an instance of the set cover problem
  scf = set_cover_filter.SetCoverFilter(
          mismatches=args.mismatches, lcf_thres=args.lcf_thres,
          blacklisted_genomes=blacklisted_genomes,
          coverage_frac=args.coverage_frac)
  #  3) Reverse complement (rc) -- add the reverse complement of
  #     each probe that remains
  rc = reverse_complement_filter.ReverseComplementFilter()
  filters = [ df, scf, rc ]

  # Don't apply the set cover filter if desired
  if args.skip_set_cover:
    filters.remove(scf)

  # Design the probes
  pb = probe_designer.ProbeDesigner(seqs, filters)
  pb.design()

  # Print the arguments and the number of final probes
  # (The final probes are stored in pb.final_probes if their
  #  sequences are desired)
  if args.limit_target_genomes:
    print args.mismatches, args.lcf_thres, args.coverage_frac, \
        args.limit_target_genomes, len(pb.final_probes)
  else:
    print args.mismatches, args.lcf_thres, args.coverage_frac, \
        len(pb.final_probes)


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("-m", "--mismatches", required=True, type=int,
      help=("Allow for this number of mismatches when determining "
            "whether a probe covers a sequence"))
  parser.add_argument("-l", "--lcf_thres", required=True, type=int,
      help=("Say that a portion of a probe covers a portion of a "
            "sequence if the two share a substring with at most "
            "'mismatches' mismatches that has length >= 'lcf_thres' "
            "bp"))
  parser.add_argument("-c", "--coverage_frac", type=float, default=1.0,
      help=("A float in [0,1] giving the fraction of each target "
            "genome that must be covered by the selected probes"))
  parser.add_argument("--skip_set_cover", dest="skip_set_cover",
      action="store_true",
      help=("Skip the set cover filter; this is useful when we "
            "wish to see the probes generated from only the "
            "duplicate and reverse complement filters, to gauge "
            "the effects of the set cover filter"))
  parser.add_argument("--blacklist_hg19", dest="blacklist_hg19",
      action="store_true",
      help=("Penalize probes based on how much of hg19 they cover"))
  parser.add_argument("--blacklist_marburg", dest="blacklist_marburg",
      action="store_true",
      help=("Penalize probes based on how much of Marburg they cover"))
  parser.add_argument("-d", "--dataset", type=str, default="ebola_zaire",
      help=("A label for the dataset to use"))
  parser.add_argument("--limit_target_genomes", type=int,
      help=("(Optional) Use only the first N target genomes in the "
            "dataset"))
  parser.add_argument('--version', '-V', action='version',
      version=version.get_version())
  args = parser.parse_args()

  main(args)
