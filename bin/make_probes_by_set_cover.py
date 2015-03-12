#!/usr/bin/env python
"""Script that runs hybseldesign to design probes for a given
dataset by treating the problem as an instance of set cover.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

import argparse
import logging
import importlib

from hybseldesign.datasets import hg19
from hybseldesign.filter import probe_designer
from hybseldesign.filter import reverse_complement_filter
from hybseldesign.filter import duplicate_filter
from hybseldesign.filter import set_cover_filter
from hybseldesign.utils import seq_io, version, log

# Set the path to the hg19 fasta file (assuming this is a Broad
# machine)
hg19.set_fasta_path("/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta")


def main(args):
  # Read the FASTA sequences
  try:
    dataset = importlib.import_module(
                 'hybseldesign.datasets.' + args.dataset)
  except ImportError:
    raise ValueError("Unknown dataset %s" % args.dataset)
  seqs = seq_io.read_fasta(dataset.fasta_path()).values()

  if args.limit_target_genomes:
    seqs = seqs[:args.limit_target_genomes]

  # Store the FASTA paths of blacklisted genomes
  blacklisted_genomes_fasta = []
  if args.blacklist_genomes:
    for bg in args.blacklist_genomes:
      try:
        dataset = importlib.import_module(
                     'hybseldesign.datasets.' + bg)
      except ImportError:
        raise ValueError("Unknown dataset %s" % bg)
      blacklisted_genomes_fasta += [ dataset.fasta_path() ]

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
          blacklisted_genomes=blacklisted_genomes_fasta,
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

  if args.output_probes:
    # Write the final probes to the file args.output_probes
    seq_io.write_probes(pb.final_probes, args.output_probes)

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
  parser.add_argument("--blacklist_genomes", nargs='+',
      help=("One or more blacklisted genomes; penalize probes based "
            "on how much of each of these genomes they cover; the "
            "label should be a dataset (e.g., 'hg19' or 'marburg')"))
  parser.add_argument("-d", "--dataset", type=str, default="ebola_zaire",
      help=("A label for the dataset to use"))
  parser.add_argument("--limit_target_genomes", type=int,
      help=("(Optional) Use only the first N target genomes in the "
            "dataset"))
  parser.add_argument("-o", "--output_probes",
      help=("(Optional) The file to which all final probes should be "
            "written; if not specified, the final probes are not "
            "written to a file"))
  parser.add_argument("--debug", dest="log_level",
      action="store_const", const=logging.DEBUG,
      default=logging.WARNING,
      help=("Debug output"))
  parser.add_argument("--verbose", dest="log_level",
      action="store_const", const=logging.INFO,
      help=("Verbose output"))
  parser.add_argument('--version', '-V', action='version',
      version=version.get_version())
  args = parser.parse_args()

  log.configure_logging(args.log_level)
  main(args)
