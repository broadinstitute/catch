#!/usr/bin/env python3
"""Run the coverage analysis module on a provided list of probe sequences."""

import argparse
import importlib
import logging
import multiprocessing
import os

from catch import coverage_analysis
from catch import probe
from catch.utils import ncbi_neighbors, seq_io, version, log

__author__ = 'Hayden Metsky <hayden@mit.edu>'


def main(args):
    # Read the genomes from FASTA sequences
    genomes_grouped = []
    genomes_grouped_names = []
    for ds in args.dataset:
        if ds.startswith('download:'):
            # Download a FASTA for an NCBI taxonomic ID
            taxid = ds[len('download:'):]
            ds_fasta_tf = ncbi_neighbors.construct_fasta_for_taxid(taxid)
            genomes_grouped += [seq_io.read_genomes_from_fasta(ds_fasta_tf.name)]
            genomes_grouped_names += ['taxid:' + str(taxid)]
            ds_fasta_tf.close()
        elif os.path.isfile(ds):
            # Process a custom fasta file with sequences
            genomes_grouped += [seq_io.read_genomes_from_fasta(ds)]
            genomes_grouped_names += [os.path.basename(ds)]
        else:
            # Process an individual dataset
            raise ValueError(("Dataset labels are no longer allowed as "
                "input. Please specify only NCBI taxonomy IDs to download "
                "(via 'download:taxid') or FASTA files. If you already "
                "specified a FASTA file, please check that the path to "
                f"'{ds}' is valid."))

    if args.limit_target_genomes:
        genomes_grouped = [genomes[:args.limit_target_genomes]
                           for genomes in genomes_grouped]

    # Set the maximum number of processes in multiprocessing pools
    if args.max_num_processes:
        probe.set_max_num_processes_for_probe_finding_pools(
            args.max_num_processes)

    # On macOS, starting with Python 3.8, new processes begin following the
    #   spawn behavior rather than fork; apparently, forking processes in
    #   macOS can cause crashes, but CATCH with older versions of
    #   Python has not experienced those issues. Parts of CATCH --
    #   especially the pools in the `probe` module -- are written to follow
    #   the fork behavior, where child processes inherit memory from their
    #   parent. When spawning processes, those child processes do not inherit
    #   memory, and CATCH crashes because the global variables are not
    #   accessible.
    # See https://github.com/python/cpython/issues/84112 for another user
    #   reporting the change in Python 3.8.
    # A quick fix is to simply change the behavior, on macOS, to fork
    #   processes. This is implemented below. A longer term, more appropriate
    #   fix might be to change the global variables accessed by children to
    #   use multiprocessing.shared_memory.SharedMemory objects.
    if os.uname().sysname == 'Darwin':
        multiprocessing.set_start_method('fork')

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
        cover_extension=args.cover_extension,
        kmer_probe_map_k=args.kmer_probe_map_k)
    analyzer.run()
    if args.write_analysis_to_tsv:
        analyzer.write_data_matrix_as_tsv(
            args.write_analysis_to_tsv)
    if args.write_sliding_window_coverage:
        analyzer.write_sliding_window_coverage(
            args.write_sliding_window_coverage)
    if args.write_probe_map_counts_to_tsv:
        analyzer.write_probe_map_counts(
                args.write_probe_map_counts_to_tsv)
    if args.print_analysis:
        analyzer.print_analysis()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    # Input data
    parser.add_argument('-d', '--dataset',
        nargs='+',
        required=True,
        help=("One or more target datasets; each can be a "
              "FASTA file or taxonomy ID to download. The "
              "format is as specified for --dataset in design.py."))
    parser.add_argument('-f', '--probes-fasta',
        required=True,
        help=("Path to a FASTA file that provides the probes (one per "
              "sequence) whose coverage should be analyzed against the "
              "genomes in the given datasets"))

    # Parameters governing hybridization
    parser.add_argument('-m', '--mismatches',
        required=True,
        type=int,
        help=("Allow for this number of mismatches when determining "
              "whether a probe covers a sequence"))
    parser.add_argument('-l', '--lcf-thres',
        required=True,
        type=int,
        help=("Say that a portion of a probe covers a portion of a "
              "sequence if the two share a substring with at most "
              "MISMATCHES mismatches that has length >= LCF_THRES "
              "bp"))
    parser.add_argument('--island-of-exact-match',
        type=int,
        default=0,
        help=("(Optional) When determining whether a probe covers a "
              "sequence, require that there be an exact match (i.e., "
              "no mismatches) of length at least ISLAND_OF_EXACT_"
              "MATCH bp between a portion of the probe and a portion "
              "of the sequence"))
    parser.add_argument('-e', '--cover-extension',
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

    # Limiting input
    parser.add_argument('--limit-target-genomes',
        type=int,
        help=("(Optional) Use only the first N target genomes in the "
              "dataset"))

    # Analysis output
    parser.add_argument('--print-analysis',
                        dest="print_analysis",
                        action="store_true",
                        help="Print analysis of the probe set's coverage")
    parser.add_argument('--write-analysis-to-tsv',
        help=("The file to which to write a TSV-formatted matrix of the "
              "probe set's coverage analysis"))
    parser.add_argument('--write-sliding-window-coverage',
        help=("The file to which to write the average coverage achieved "
              "by the probe set within sliding windows of each target "
              "genome"))
    parser.add_argument('--write-probe-map-counts-to-tsv',
        help=("The file to which to write a TSV-formatted list of the "
              "number of sequences each probe maps to. This explicitly "
              "does not count reverse complements."))

    # Technical adjustments
    def check_max_num_processes(val):
        ival = int(val)
        if ival >= 1:
            return ival
        else:
            raise argparse.ArgumentTypeError(("MAX_NUM_PROCESSES must be "
                                              "an int >= 1"))
    parser.add_argument('--max-num-processes',
        type=check_max_num_processes,
        help=("(Optional) An int >= 1 that gives the maximum number of "
              "processes to use in multiprocessing pools; uses min(number "
              "of CPUs in the system, MAX_NUM_PROCESSES) processes"))
    parser.add_argument('--kmer-probe-map-k',
        type=int,
        default=10,
        help=("(Optional) Use this value (KMER_PROBE_LENGTH_K) as the "
              "k-mer length when constructing a map of k-mers to the probes "
              "that contain these k-mers. This map is used when mapping "
              "the given probes to target sequences and the k-mers serve "
              "as seeds for calculating whether a probe 'covers' "
              "a subsequence. The value should be sufficiently less than "
              "the probe length (PROBE_LENGTH) so that it can find mappings "
              "even when the candidate probe and target sequence are "
              "divergent. In particular, CATCH will try to find a value k >= "
              "KMER_PROBE_LENGTH_K (by default, >=10) such that k divides "
              "PROBE_LENGTH and k < PROBE_LENGTH / MISMATCHES (if "
              "MISMATCHES=0, then k=PROBE_LENGTH). It will then use this "
              "k as the k-mer length in mappings; if no such k exists, it "
              "will use a randomized approach with KMER_PROBE_LENGTH_K as "
              "the k-mer length."))

    # Logging levels and version
    parser.add_argument('--debug',
        dest="log_level",
        action="store_const",
        const=logging.DEBUG,
        default=logging.WARNING,
        help=("Debug output"))
    parser.add_argument('--verbose',
        dest="log_level",
        action="store_const",
        const=logging.INFO,
        help=("Verbose output"))
    parser.add_argument('-V', '--version',
        action='version',
        version=version.get_version())

    args = parser.parse_args()

    log.configure_logging(args.log_level)
    main(args)
