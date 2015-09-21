#!/usr/bin/env python
"""Design probes by treating the problem as an instance of set cover."""

import argparse
import importlib
import logging

from hybseldesign import coverage_analysis
from hybseldesign import probe
from hybseldesign.datasets import hg19
from hybseldesign.filter import adapter_filter
from hybseldesign.filter import duplicate_filter
from hybseldesign.filter import fasta_filter
from hybseldesign.filter import probe_designer
from hybseldesign.filter import reverse_complement_filter
from hybseldesign.filter import set_cover_filter
from hybseldesign.utils import seq_io, version, log

__author__ = 'Hayden Metsky <hayden@mit.edu>'

# Set the path to the hg19 fasta file (assuming this is a Broad machine)
hg19.add_fasta_path(
    "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta")


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

    # Store the FASTA paths of blacklisted genomes
    blacklisted_genomes_fasta = []
    if args.blacklist_genomes:
        for bg in args.blacklist_genomes:
            try:
                dataset = importlib.import_module(
                    'hybseldesign.datasets.' + bg)
            except ImportError:
                raise ValueError("Unknown dataset %s" % bg)
            for fp in dataset.fasta_paths:
                blacklisted_genomes_fasta += [fp]

    # Set the maximum number of processes in multiprocessing pools
    if args.max_num_processes:
        probe.set_max_num_processes_for_probe_finding_pools(
            args.max_num_processes)

    # Setup the filters
    # The filters we use are, in order:
    #  1) Duplicate filter (df) -- condense all candidate probes that
    #     are identical down to one; this is not necessary for
    #     correctness, as the set cover filter achieves the same task
    #     implicitly, but it does significantly lower runtime by
    #     decreasing the input size to the set cover filter
    df = duplicate_filter.DuplicateFilter()
    #  2) Set cover filter (scf) -- solve the problem by treating it as
    #     an instance of the set cover problem
    scf = set_cover_filter.SetCoverFilter(
        mismatches=args.mismatches,
        lcf_thres=args.lcf_thres,
        mismatches_tolerant=args.mismatches_tolerant,
        lcf_thres_tolerant=args.lcf_thres_tolerant,
        identify=args.identify,
        blacklisted_genomes=blacklisted_genomes_fasta,
        coverage=args.coverage,
        cover_extension=args.cover_extension,
        cover_groupings_separately=args.cover_groupings_separately)
    #  3) Adapter filter (af) -- add adapters to both the 5' and 3' ends
    #     of each probe
    af = adapter_filter.AdapterFilter(tuple(args.adapter_a),
                                      tuple(args.adapter_b),
                                      mismatches=args.mismatches,
                                      lcf_thres=args.lcf_thres)
    #  4) Reverse complement (rc) -- add the reverse complement of each
    #     probe that remains
    rc = reverse_complement_filter.ReverseComplementFilter()
    filters = [df, scf, af, rc]

    # Add a FASTA filter to the beginning if desired
    if args.filter_from_fasta:
        ff = fasta_filter.FastaFilter(args.filter_from_fasta,
                                      skip_reverse_complements=True)
        filters.insert(0, ff)

    # Don't apply the set cover filter if desired
    if args.skip_set_cover:
        filters.remove(scf)

    # Don't add adapters if desired
    if args.skip_adapters:
        filters.remove(af)

    # Design the probes
    pb = probe_designer.ProbeDesigner(genomes_grouped, filters)
    pb.design()

    if args.output_probes:
        # Write the final probes to the file args.output_probes
        seq_io.write_probe_fasta(pb.final_probes, args.output_probes)

    if (args.print_analysis or args.write_analysis_to_tsv or
            args.write_sliding_window_coverage):
        analyzer = coverage_analysis.Analyzer(
            pb.final_probes,
            genomes_grouped,
            genomes_grouped_names,
            mismatches=args.mismatches,
            lcf_thres=args.lcf_thres,
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
    else:
        # Just print the number of probes
        print(len(pb.final_probes))


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
        "-mt", "--mismatches_tolerant",
        type=int,
        help=("(Optional) A more tolerant value for '--mismatches'; "
              "this should be greater than the value of '--mismatches'. "
              "Allows for capturing more possible hybridizations "
              "(i.e., more sensitivity) when designing probes for "
              "identification or when genomes are blacklisted."))
    parser.add_argument(
        "-lt", "--lcf_thres_tolerant",
        type=int,
        help=("(Optional) A more tolerant value for '--lcf_thres'; "
              "this should be less than the value of '--lcf_thres'. "
              "Allows for capturing more possible hybridizations "
              "(i.e., more sensitivity) when designing probes for "
              "identification or when genomes are blacklisted."))
    parser.add_argument(
        "-i", "--identify",
        dest="identify",
        action="store_true",
        help=("Design probes meant to make it possible to identify "
              "nucleic acid from a particular input dataset against "
              "the other datasets; when set, the coverage should "
              "generally be small"))

    def check_coverage(val):
        fval = float(val)
        ival = int(fval)
        if fval >= 0 and fval <= 1:
            # a float in [0,1] giving fractional coverage
            return fval
        elif fval > 1 and fval == ival:
            # an int > 1 giving number of bp to cover
            return ival
        else:
            raise argparse.ArgumentTypeError(("%s is an invalid coverage "
                                              "value") % val)

    parser.add_argument(
        "-c", "--coverage",
        type=check_coverage,
        default=1.0,
        help=("If this is a float in [0,1], it gives the fraction of "
              "each target genome that must be covered by the selected "
              "probes; if this is an int > 1, it gives the number of "
              "bp of each target genome that must be covered by the "
              "selected probes"))
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
        "--skip_set_cover",
        dest="skip_set_cover",
        action="store_true",
        help=("Skip the set cover filter; this is useful when we "
              "wish to see the probes generated from only the "
              "duplicate and reverse complement filters, to gauge "
              "the effects of the set cover filter"))
    parser.add_argument("--skip_adapters",
                        dest="skip_adapters",
                        action="store_true",
                        help=("Do not add adapters to the ends of probes"))
    parser.add_argument(
        "--filter_from_fasta",
        help=("(Optional) A FASTA file from which to select candidate probes. "
              "Before running any other filters, keep only the candidate "
              "probes that are equal to sequences in the file and remove "
              "all probes not equal to any of these sequences. This, by "
              "default, ignores sequences in the file whose header contains "
              "the string 'reverse complement'; that is, if there is some "
              "probe with sequence S, it may be filtered out (even if there "
              "is a sequence S in the file) if the header of S in the file "
              "contains 'reverse complement'. This is useful if we already "
              "have probes decided by the set cover filter, but simply "
              "want to process them further by, e.g., adding adapters or "
              "running a coverage analysis. For example, if we have already "
              "run the time-consuming set cover filter and have a FASTA "
              "containing those probes, we can provide a path to that "
              "FASTA file for this argument, and also provide the "
              "'--skip_set_cover' argument, in order to add adapters to "
              "those probes without having to re-run the set cover filter."))
    parser.add_argument(
        "--blacklist_genomes",
        nargs='+',
        help=("One or more blacklisted genomes; penalize probes based "
              "on how much of each of these genomes they cover; the "
              "label should be a dataset (e.g., 'hg19' or 'marburg')"))
    parser.add_argument(
        "--cover_groupings_separately",
        dest="cover_groupings_separately",
        action="store_true",
        help=("Run a separate instance of set cover with the target genomes "
              "from each grouping and pool (union) the resulting probes. "
              "When set, the software will run faster than when not set, but "
              "it may yield more probes than when it is not set."))
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
    parser.add_argument(
        "--adapter_a",
        nargs=2,
        default=['ATACGCCATGCTGGGTCTCC', 'CGTACTTGGGAGTCGGCCAT'],
        help=("(Optional) Custom A adapter to use; two ordered arguments x "
              "and y such that x is the A adapter sequence to place on the "
              "5' end of a probe and y is the A adapter sequence to place "
              "on the 3' end of a probe"))
    parser.add_argument(
        "--adapter_b",
        nargs=2,
        default=['AGGCCCTGGCTGCTGATATG', 'GACCTTTTGGGACAGCGGTG'],
        help=("(Optional) Custom B adapter to use; two ordered arguments x "
              "and y such that x is the B adapter sequence to place on the "
              "5' end of a probe and y is the B adapter sequence to place "
              "on the 3' end of a probe"))
    parser.add_argument(
        "-o", "--output_probes",
        help=("(Optional) The file to which all final probes should be "
              "written; if not specified, the final probes are not "
              "written to a file"))
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
