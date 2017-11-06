#!/usr/bin/env python3
"""Design probes in naive ways.

This offers a few options to design probes using naive methods.
This is mainly used for comparison with the probes generated
by bin/make_probes.py.
"""

import argparse
import importlib
import logging
import random

from hybseldesign import coverage_analysis
from hybseldesign.filter import duplicate_filter
from hybseldesign.filter import naive_redundant_filter
from hybseldesign.filter import dominating_set_filter
from hybseldesign.filter import probe_designer
from hybseldesign.filter import reverse_complement_filter
from hybseldesign.utils import seq_io, version, log

__author__ = 'Hayden Metsky <hayden@mit.edu>'


def main(args):
    # Read the FASTA sequences
    try:
        dataset = importlib.import_module(
            'hybseldesign.datasets.' + args.dataset)
    except ImportError:
        raise ValueError("Unknown dataset %s" % ds)
    seqs = [seq_io.read_dataset_genomes(dataset)]

    if (args.limit_target_genomes and
            args.limit_target_genomes_randomly_with_replacement):
        raise Exception(("Cannot --limit_target_genomes and "
                         "--limit_target_genomes_randomly_with_replacement at "
                         "the same time"))
    elif args.limit_target_genomes:
        seqs = [genomes[:args.limit_target_genomes] for genomes in seqs]
    elif args.limit_target_genomes_randomly_with_replacement:
        k = args.limit_target_genomes_randomly_with_replacement
        seqs = [random.choices(genomes, k=k) for genomes in seqs]

    # Setup the filters needed for replication
    filters = []
    # The filters we use are, in order:

    #  Duplicate filter (df) -- condense all candidate probes that
    #  are identical down to one; this is not necessary for
    #  correctness, as the naive redundant filter achieves the same
    #  task implicitly, but it does significantly lower runtime by
    #  decreasing the input size to the naive redundant filter
    df = duplicate_filter.DuplicateFilter()
    filters += [df]

    if args.naive_redundant_filter and args.dominating_set_filter:
        raise Exception(("Cannot use both 'naive_redundant_filter' and "
            "'dominating_set_filter' at the same time. (You could "
            "of course do one after the other, but it was probably "
            "a mistake to specify both.)"))
    elif args.naive_redundant_filter or args.dominating_set_filter:
        if args.naive_redundant_filter:
            # Naive redundant filter -- execute a greedy algorithm to
            # condense 'similar' probes down to one
            shift, mismatches = args.naive_redundant_filter
            filt_class = naive_redundant_filter.NaiveRedundantFilter
        if args.dominating_set_filter:
            # Dominating set filter (dsf) -- construct a graph where each
            # node is a probe and edges connect 'similar' probes; then
            # approximate the smallest dominating set
            shift, mismatches = args.dominating_set_filter
            filt_class = dominating_set_filter.DominatingSetFilter
        # Construct a function to determine whether two probes are
        # redundant, and then instantiate the appropriate filter
        redundant_fn = naive_redundant_filter.redundant_shift_and_mismatch_count(
            shift=shift,
            mismatch_thres=mismatches)
        filt = filt_class(redundant_fn)
        filters += [filt]

    if not args.skip_reverse_complements:
        # Reverse complement (rc) -- add the reverse complement of
        # each probe as a candidate
        rc = reverse_complement_filter.ReverseComplementFilter()
        filters += [rc]

    # Design the probes
    pb = probe_designer.ProbeDesigner(seqs, filters,
                                      probe_length=args.probe_length,
                                      probe_stride=args.probe_stride)
    pb.design()

    if args.print_analysis:
        if args.naive_redundant_filter or args.dominating_set_filter:
            mismatch_thres = mismatches
        else:
            mismatch_thres = 0
        analyzer = coverage_analysis.Analyzer(pb.final_probes,
                                              mismatch_thres,
                                              args.probe_length,
                                              seqs,
                                              [args.dataset])
        analyzer.run()
        analyzer.print_analysis()
    else:
        # Just print the number of probes
        print(len(pb.final_probes))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-pl", "--probe_length",
        type=int,
        default=100,
        help=("(Optional) The number of bp in each probe"))
    parser.add_argument(
        "-ps", "--probe_stride",
        type=int,
        default=50,
        help=("(Optional) Generate candidate probes from the input "
              "that are separated by this number of bp"))
    parser.add_argument("-d", "--dataset",
                        required=True,
                        help="Label for the target dataset")
    parser.add_argument(
        "-nrf", "--naive_redundant_filter",
        nargs=2,
        type=int,
        help=("Args: <shift> <mismatches>. Use the 'naive redundant "
              "filter' -- i.e., iterate through a list of probes and, "
              "for each probe p, remove the following probes that are "
              "redundant to p. Deem one probe redundant to another if, "
              "as it is shifted from -shift to +shift bp relative "
              "to the other, the minimum number of mismatches between "
              "them is <= mismatches. This is meant to replicate a "
              "prior approach for design."))
    parser.add_argument(
        "-dsf", "--dominating_set_filter",
        nargs=2,
        type=int,
        help=("Args: <shift> <mismatches>. Use the 'dominating set "
              "filter' -- i.e., filter redundant probes by constructing "
              "a graph connecting redundant probes and approximating "
              "the smallest dominating set. Deem one probe redundant "
              "to another if, as it is shifted from -shift to +shift bp "
              "relative to the other, the minimum number of mismatches "
              "between them is <= mismatches"))
    parser.add_argument("--skip_reverse_complements",
                        dest="skip_reverse_complements",
                        action="store_true",
                        help=("Do not add to the output the reverse "
                              "complement of each probe"))
    parser.add_argument(
        "--limit_target_genomes",
        type=int,
        help=("(Optional) Use only the first N target genomes in the "
              "dataset"))
    parser.add_argument(
        "--limit_target_genomes_randomly_with_replacement",
        type=int,
        help=("(Optional) Randomly select N target genomes in the "
              "dataset with replacement"))
    parser.add_argument("--print_analysis",
                        dest="print_analysis",
                        action="store_true",
                        help="Print analysis of the probe set's coverage")
    parser.add_argument("--debug",
                        dest="log_level",
                        action="store_const",
                        const=logging.DEBUG,
                        default=logging.WARNING,
                        help=("Debug output"))
    parser.add_argument("--verbose",
                        dest="log_level",
                        action="store_const",
                        const=logging.INFO,
                        help=("Verbose output"))
    parser.add_argument('--version', '-V',
                        action='version',
                        version=version.get_version())
    args = parser.parse_args()

    log.configure_logging(args.log_level)
    main(args)
