#!/usr/bin/env python3
"""Run the coverage analysis module on a provided list of probe sequences."""

import argparse
import importlib
import logging

from hybseldesign import coverage_analysis
from hybseldesign import probe
from hybseldesign.utils import seq_io, version, log

__author__ = 'Hayden Metsky <hayden@mit.edu>'


def main(args):
    # Read the genomes from FASTA sequences
    genomes_grouped = []
    genomes_grouped_names = []
    for ds in args.dataset:
        try:
            dataset = importlib.import_module('hybseldesign.datasets.' + ds)
        except ImportError:
            raise ValueError("Unknown dataset %s" % ds)
        genomes_grouped += [seq_io.read_dataset_genomes(dataset)]
        genomes_grouped_names += [ds]

    if args.limit_target_genomes:
        genomes_grouped = [genomes[:args.limit_target_genomes]
                           for genomes in genomes_grouped]

    # Set the maximum number of processes in multiprocessing pools
    if args.max_num_processes:
        probe.set_max_num_processes_for_probe_finding_pools(
            args.max_num_processes)

    # Read the FASTA file of probes
    fasta = seq_io.read_fasta(args.probes_fasta)
    probes = [probe.Probe.from_str(seq) for _, seq in fasta.items()]

    # Run the coverage analyzer
    analyzer = coverage_analysis.Analyzer(
        probes,
        args.mismatches,
        args.lcf_thres,
        genomes_grouped,
        genomes_grouped_names,
        island_of_exact_match=args.island_of_exact_match,
        cover_extension=args.cover_extension)
    analyzer.run()
    if args.write_analysis_to_tsv:
        analyzer.write_data_matrix_as_tsv(
            args.write_analysis_to_tsv)
    if args.write_sliding_window_coverage:
        analyzer.write_sliding_window_coverage(
            args.write_sliding_window_coverage)
    if args.print_analysis:
        analyzer.print_analysis()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-m", "--mismatches",
        required=True,
        type=int,
        help=("Allow for this number of mismatches when determining "
              "whether a probe covers a sequence"))
    parser.add_argument(
        "-l", "--lcf_thres",
        required=True,
        type=int,
        help=("Say that a portion of a probe covers a portion of a "
              "sequence if the two share a substring with at most "
              "'mismatches' mismatches that has length >= 'lcf_thres' "
              "bp"))
    parser.add_argument(
        "--island_of_exact_match",
        type=int,
        default=0,
        help=("(Optional) When determining whether a probe covers a "
              "sequence, require that there be an exact match (i.e., "
              "no mismatches) of length at least 'island_of_exact_"
              "match' bp between a portion of the probe and a portion "
              "of the sequence"))
    parser.add_argument(
        "-e", "--cover_extension",
        type=int,
        default=0,
        help=("Extend the coverage of each side of a probe by this number "
              "of bp. That is, a probe covers a region that consists of the "
              "portion of a sequence it hybridizes to, as well as this "
              "number of bp on each side of that portion. This is useful "
              "in modeling hybrid selection, where a probe hybridizes to"
              "a fragment that includes the region targeted by the probe, "
              "along with surrounding portions of the sequence. Increasing "
              "its value should reduce the number of probes required to "
              "achieve the desired coverage."))
    parser.add_argument(
        '-f', '--probes_fasta',
        required=True,
        help=("Path to a FASTA file that provides the probes (one per "
              "sequence) whose coverage should be analyzed against the "
              "genomes in the given datasets"))
    parser.add_argument("-d", "--dataset",
                        nargs='+',
                        required=True,
                        help=("Labels for one or more target datasets (e.g., "
                              "one label per species)"))
    parser.add_argument(
        "--limit_target_genomes",
        type=int,
        help=("(Optional) Use only the first N target genomes in the "
              "dataset"))
    parser.add_argument("--print_analysis",
                        dest="print_analysis",
                        action="store_true",
                        help="Print analysis of the probe set's coverage")
    parser.add_argument(
        "--write_analysis_to_tsv",
        help=("The file to which to write a TSV-formatted matrix of the "
              "probe set's coverage analysis"))
    parser.add_argument(
        "--write_sliding_window_coverage",
        help=("The file to which to write the average coverage achieved "
              "by the probe set within sliding windows of each target "
              "genome"))

    def check_max_num_processes(val):
        ival = int(val)
        if ival >= 1:
            return ival
        else:
            raise argparse.ArgumentTypeError(("max_num_processes must be "
                                              "an int >= 1"))

    parser.add_argument(
        "--max_num_processes",
        type=check_max_num_processes,
        help=("(Optional) An int >= 1 that gives the maximum number of "
              "processes to use in multiprocessing pools; uses min(number "
              "of CPUs in the system, max_num_processes) processes"))
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
