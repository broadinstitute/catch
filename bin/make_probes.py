#!/usr/bin/env python
"""Design probes by treating the problem as an instance of set cover."""

import argparse
import importlib
import logging

from hybseldesign import coverage_analysis
from hybseldesign.datasets import hg19
from hybseldesign.filter import adapter_filter
from hybseldesign.filter import duplicate_filter
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
            if dataset.is_multi_chr():
                for fp in dataset.fasta_paths:
                    blacklisted_genomes_fasta += [fp]
            else:
                blacklisted_genomes_fasta += [dataset.fasta_path]

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
        cover_groupings_separately=args.cover_groupings_separately)
    #  3) Adapter filter (af) -- add adapters to both the 5' and 3' ends
    #     of each probe
    af = adapter_filter.AdapterFilter(mismatches=args.mismatches,
                                      lcf_thres=args.lcf_thres)
    #  4) Reverse complement (rc) -- add the reverse complement of each
    #     probe that remains
    rc = reverse_complement_filter.ReverseComplementFilter()
    filters = [df, scf, af, rc]

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

    if args.print_analysis:
        analyzer = coverage_analysis.Analyzer(pb.final_probes,
                                              genomes_grouped,
                                              genomes_grouped_names,
                                              mismatches=args.mismatches,
                                              lcf_thres=args.lcf_thres)
        analyzer.run()
        analyzer.print_analysis()
    else:
        # Just print the number of probes
        print len(pb.final_probes)


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
        "-o", "--output_probes",
        help=("(Optional) The file to which all final probes should be "
              "written; if not specified, the final probes are not "
              "written to a file"))
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
