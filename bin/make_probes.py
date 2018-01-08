#!/usr/bin/env python3
"""Design probes by treating the problem as an instance of set cover."""

import argparse
import importlib
import logging
import os
import random

from hybseldesign import coverage_analysis
from hybseldesign import probe
from hybseldesign.datasets import hg19
from hybseldesign.filter import adapter_filter
from hybseldesign.filter import duplicate_filter
from hybseldesign.filter import fasta_filter
from hybseldesign.filter import n_expansion_filter
from hybseldesign.filter import probe_designer
from hybseldesign.filter import reverse_complement_filter
from hybseldesign.filter import set_cover_filter
from hybseldesign.utils import seq_io, version, log

__author__ = 'Hayden Metsky <hayden@mit.edu>'


def main(args):
    logger = logging.getLogger(__name__)

    # Read the genomes from FASTA sequences
    genomes_grouped = []
    genomes_grouped_names = []
    for ds in args.dataset:
        if ds.startswith('collection:'):
            # Process a collection of datasets
            collection_name = ds[len('collection:'):]
            try:
                collection = importlib.import_module(
                    'hybseldesign.datasets.collections.' + collection_name)
            except ImportError:
                raise ValueError("Unknown dataset collection %s" %
                                 collection_name)
            for name, dataset in collection.import_all():
                genomes_grouped += [seq_io.read_dataset_genomes(dataset)]
                genomes_grouped_names += [name]
        elif os.path.isfile(ds):
            # Process a custom fasta file with sequences
            genomes_grouped += [seq_io.read_genomes_from_fasta(ds)]
            genomes_grouped_names += [os.path.basename(ds)]
        else:
            # Process an individual dataset
            try:
                dataset = importlib.import_module(
                            'hybseldesign.datasets.' + ds)
            except ImportError:
                raise ValueError("Unknown file or dataset '%s'" % ds)
            genomes_grouped += [seq_io.read_dataset_genomes(dataset)]
            genomes_grouped_names += [ds]

    if (args.limit_target_genomes and
            args.limit_target_genomes_randomly_with_replacement):
        raise Exception(("Cannot --limit-target-genomes and "
                         "--limit-target-genomes-randomly-with-replacement at "
                         "the same time"))
    elif args.limit_target_genomes:
        genomes_grouped = [genomes[:args.limit_target_genomes]
                           for genomes in genomes_grouped]
    elif args.limit_target_genomes_randomly_with_replacement:
        k = args.limit_target_genomes_randomly_with_replacement
        genomes_grouped = [random.choices(genomes, k=k)
                           for genomes in genomes_grouped]

    # Store the FASTA paths of blacklisted genomes
    blacklisted_genomes_fasta = []
    if args.blacklist_genomes:
        for bg in args.blacklist_genomes:
            if os.path.isfile(bg):
                # Process a custom fasta file with sequences
                blacklisted_genomes_fasta += [bg]
            else:
                # Process an individual dataset
                try:
                    dataset = importlib.import_module(
                        'hybseldesign.datasets.' + bg)
                except ImportError:
                    raise ValueError("Unknown file or dataset '%s'" % bg)
                for fp in dataset.fasta_paths:
                    blacklisted_genomes_fasta += [fp]

    # Setup and verify parameters related to probe length
    if not args.lcf_thres:
        args.lcf_thres = args.probe_length
    if args.probe_stride > args.probe_length:
        logger.warning(("PROBE_STRIDE (%d) is greater than PROBE_LENGTH "
                        "(%d), which is usually undesirable and may lead "
                        "to undefined behavior"),
                        args.probe_stride, args.probe_length)
    if args.lcf_thres > args.probe_length:
        logger.warning(("LCF_THRES (%d) is greater than PROBE_LENGTH "
                        "(%d), which is usually undesirable and may lead "
                        "to undefined behavior"),
                        args.lcf_thres, args.probe_length)
    if args.island_of_exact_match > args.probe_length:
        logger.warning(("ISLAND_OF_EXACT_MATCH (%d) is greater than "
                        "PROBE_LENGTH (%d), which is usually undesirable "
                        "and may lead to undefined behavior"),
                        args.island_of_exact_match, args.probe_length)

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
        island_of_exact_match=args.island_of_exact_match,
        mismatches_tolerant=args.mismatches_tolerant,
        lcf_thres_tolerant=args.lcf_thres_tolerant,
        island_of_exact_match_tolerant=args.island_of_exact_match_tolerant,
        identify=args.identify,
        blacklisted_genomes=blacklisted_genomes_fasta,
        coverage=args.coverage,
        cover_extension=args.cover_extension,
        cover_groupings_separately=args.cover_groupings_separately,
        kmer_probe_map_use_native_dict=args.use_native_dict_when_finding_tolerant_coverage)
    #  3) Adapter filter (af) -- add adapters to both the 5' and 3' ends
    #     of each probe
    af = adapter_filter.AdapterFilter(tuple(args.adapter_a),
                                      tuple(args.adapter_b),
                                      mismatches=args.mismatches,
                                      lcf_thres=args.lcf_thres,
                                      island_of_exact_match=\
                                        args.island_of_exact_match)
    filters = [df, scf, af]

    #  4) Reverse complement (rc) -- add the reverse complement of each
    #     probe that remains (if requested)
    if args.add_reverse_complements:
        rc = reverse_complement_filter.ReverseComplementFilter()
        filters += [rc]

    # Add a FASTA filter to the beginning if desired
    if args.filter_from_fasta:
        ff = fasta_filter.FastaFilter(args.filter_from_fasta,
                                      skip_reverse_complements=True)
        filters.insert(0, ff)

    # Add an N expansion filter just before the reverse complement
    # filter if desired
    if args.expand_n:
        nef = n_expansion_filter.NExpansionFilter()
        if args.add_reverse_complements:
            rc_pos = filters.index(rc)
            filters.insert(rc_pos, nef)
        else:
            # There is no reverse complement filter, so that it
            # to the end
            filters += [nef]

    # Don't apply the set cover filter if desired
    if args.skip_set_cover:
        filters.remove(scf)

    # Don't add adapters if desired
    if args.skip_adapters:
        filters.remove(af)

    # Design the probes
    pb = probe_designer.ProbeDesigner(genomes_grouped, filters,
                                      probe_length=args.probe_length,
                                      probe_stride=args.probe_stride,
                                      allow_small_seqs=args.small_seq_min)
    pb.design()

    if args.output_probes:
        # Write the final probes to the file args.output_probes
        seq_io.write_probe_fasta(pb.final_probes, args.output_probes)

    if (args.print_analysis or args.write_analysis_to_tsv or
            args.write_sliding_window_coverage):
        rc_too = False if args.skip_reverse_complements else True
        analyzer = coverage_analysis.Analyzer(
            pb.final_probes,
            args.mismatches,
            args.lcf_thres,
            genomes_grouped,
            genomes_grouped_names,
            island_of_exact_match=args.island_of_exact_match,
            cover_extension=args.cover_extension,
            rc_too=rc_too)
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

    # Input data
    parser.add_argument('dataset',
        nargs='+',
        help=("One or more target datasets (e.g., one per species). If "
              "DATASET is a path to a file, then that file is treated as "
              "a FASTA file and its sequences are read. Otherwise, it "
              "is assumed that this is a label for a dataset included "
              "in this package (e.g., 'zika'). If the label starts with "
              "'collection:' (e.g., 'collection:viruses_with_human_host'), "
              "then this reads from an available collection of datasets."))

    # Parameters on probe length and stride
    parser.add_argument('-pl', '--probe-length',
        type=int,
        default=100,
        help=("(Optional) Make probes be PROBE_LENGTH nt long"))
    parser.add_argument('-ps', '--probe-stride',
        type=int,
        default=50,
        help=("(Optional) Generate candidate probes from the input "
              "that are separated by PROBE_STRIDE nt"))

    # Parameters governing probe hybridization
    parser.add_argument('-m', '--mismatches',
        type=int,
        default=0,
        help=("(Optional) Allow for MISMATCHES mismatches when determining "
              "whether a probe covers a sequence"))
    parser.add_argument('-l', '--lcf-thres',
        type=int,
        help=("(Optional) Say that a portion of a probe covers a portion "
              "of a sequence if the two share a substring with at most "
              "MISMATCHES mismatches that has length >= LCF_THRES "
              "nt; if unspecified, this is set to PROBE_LENGTH"))
    parser.add_argument('--island-of-exact-match',
        type=int,
        default=0,
        help=("(Optional) When determining whether a probe covers a "
              "sequence, require that there be an exact match (i.e., "
              "no mismatches) of length at least ISLAND_OF_EXACT_"
              "MATCH nt between a portion of the probe and a portion "
              "of the sequence"))

    # Desired coverage of target genomes
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
    parser.add_argument('-c', '--coverage',
        type=check_coverage,
        default=1.0,
        help=("If this is a float in [0,1], it gives the fraction of "
              "each target genome that must be covered by the selected "
              "probes; if this is an int > 1, it gives the number of "
              "bp of each target genome that must be covered by the "
              "selected probes"))

    # Amount of cover extension to assume
    parser.add_argument('-e', '--cover-extension',
        type=int,
        default=0,
        help=("Extend the coverage of each side of a probe by COVER_EXTENSION "
              "nt. That is, a probe covers a region that consists of the "
              "portion of a sequence it hybridizes to, as well as this "
              "number of nt on each side of that portion. This is useful "
              "in modeling hybrid selection, where a probe hybridizes to"
              "a fragment that includes the region targeted by the probe, "
              "along with surrounding portions of the sequence. Increasing "
              "its value should reduce the number of probes required to "
              "achieve the desired coverage."))

    # Differential identification and blacklisting
    parser.add_argument('-i', '--identify',
        dest="identify",
        action="store_true",
        help=("Design probes meant to make it possible to identify "
              "nucleic acid from a particular input dataset against "
              "the other datasets; when set, the coverage should "
              "generally be small"))
    parser.add_argument('--blacklist-genomes',
        nargs='+',
        help=("One or more blacklisted genomes; penalize probes based "
              "on how much of each of these genomes they cover. If "
              "the value is a path to a file, then that file is treated "
              "as a FASTA file and its sequences are read. Otherwise, "
              "it is assumed that this is a label for a dataset included "
              "in this package (e.g., 'zika')."))
    parser.add_argument('-mt', '--mismatches-tolerant',
        type=int,
        help=("(Optional) A more tolerant value for 'mismatches'; "
              "this should be greater than the value of MISMATCHES. "
              "Allows for capturing more possible hybridizations "
              "(i.e., more sensitivity) when designing probes for "
              "identification or when genomes are blacklisted."))
    parser.add_argument('-lt', '--lcf-thres-tolerant',
        type=int,
        help=("(Optional) A more tolerant value for 'lcf_thres'; "
              "this should be less than LCF_THRES. "
              "Allows for capturing more possible hybridizations "
              "(i.e., more sensitivity) when designing probes for "
              "identification or when genomes are blacklisted."))
    parser.add_argument('--island-of-exact-match-tolerant',
        type=int,
        default=0,
        help=("(Optional) A more tolerant value for 'island_of_"
              "exact_match'; this should be less than ISLAND_OF_ "
              "EXACT_MATCH. Allows for capturing more "
              "possible hybridizations (i.e., more sensitivity) "
              "when designing probes for identification or when "
              "genomes are blacklisted."))

    # Outputting probe sequences and coverage analyses
    parser.add_argument('-o', '--output-probes',
        help=("(Optional) The file to which all final probes should be "
              "written; if not specified, the final probes are not "
              "written to a file"))
    parser.add_argument('--print-analysis',
        dest="print_analysis",
        action="store_true",
        help="Print analysis of the probe set's coverage")
    parser.add_argument('--write-analysis-to-tsv',
        help=("(Optional) The file to which to write a TSV-formatted matrix "
              "of the probe set's coverage analysis"))
    parser.add_argument('--write-sliding-window-coverage',
        help=("(Optional) The file to which to write the average coverage "
              "achieved by the probe set within sliding windows of each "
              "target genome"))

    # Accepting probes as input and skipping set cover process
    parser.add_argument('--filter-from-fasta',
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
              "--skip-set-cover argument, in order to add adapters to "
              "those probes without having to re-run the set cover filter."))
    parser.add_argument('--skip-set-cover',
        dest="skip_set_cover",
        action="store_true",
        help=("Skip the set cover filter; this is useful when we "
              "wish to see the probes generated from only the "
              "duplicate and reverse complement filters, to gauge "
              "the effects of the set cover filter"))

    # Adding adapters
    parser.add_argument('--skip-adapters',
        dest="skip_adapters",
        action="store_true",
        help=("Do not add adapters to the ends of probes"))
    parser.add_argument('--adapter-a',
        nargs=2,
        default=['ATACGCCATGCTGGGTCTCC', 'CGTACTTGGGAGTCGGCCAT'],
        help=("(Optional) Args: <X> <Y>; Custom A adapter to use; two ordered "
              "where X is the A adapter sequence to place on the 5' end of "
              "a probe and Y is the A adapter sequence to place on the 3' "
              "end of a probe"))
    parser.add_argument('--adapter-b',
        nargs=2,
        default=['AGGCCCTGGCTGCTGATATG', 'GACCTTTTGGGACAGCGGTG'],
        help=("(Optional) Args: <X> <Y>; Custom B adapter to use; two ordered "
              "where X is the B adapter sequence to place on the 5' end of "
              "a probe and Y is the B adapter sequence to place on the 3' "
              "end of a probe"))

    # Adjusting probe output
    parser.add_argument('--add-reverse-complements',
        dest="add_reverse_complements",
        action="store_true",
        help=("Add to the output the reverse complement of each probe"))
    parser.add_argument('--expand-n',
        dest="expand_n",
        action="store_true",
        help=("Expand each probe so that 'N' bases are replaced by real "
              "bases; for example, the probe 'ANA' would be replaced "
              "with the probes 'AAA', 'ATA', 'ACA', and 'AGA'; this is "
              "done combinatorially across all 'N' bases in a probe, and "
              "thus the number of new probes grows exponentially with the "
              "number of 'N' bases in a probe"))

    # Limiting input
    parser.add_argument('--limit-target-genomes',
        type=int,
        help=("(Optional) Use only the first LIMIT_TARGET_GENOMES target "
              "genomes in the dataset"))
    parser.add_argument('--limit-target-genomes-randomly-with-replacement',
        type=int,
        help=("(Optional) Randomly select LIMIT_TARGET_GENOMES_RANDOMLY_"
              "WITH_REPLACMENT target genomes in the dataset with "
              "replacement"))

    # Technical adjustments
    parser.add_argument('--cover-groupings-separately',
        dest="cover_groupings_separately",
        action="store_true",
        help=("Run a separate instance of set cover with the target genomes "
              "from each grouping and pool (union) the resulting probes. "
              "When set, the software will run faster than when not set, but "
              "it may yield more probes than when it is not set."))
    parser.add_argument('--small-seq-min',
        type=int,
        help=("(Optional) If set, allow sequences as input that are "
              "shorter than PROBE_LENGTH (when not set, the program will "
              "error on such input). SMALL_SEQ_MIN is the "
              "minimum sequence length that should be accepted as input. "
              "When a sequence is less than PROBE_LENGTH, a candidate "
              "probe is created that is equal to the sequence; thus, "
              "the output probes may have different lengths. Note that, "
              "when this is set, it might be a good idea to also set "
              "LCF_THRES to be a value smaller than PROBE_LENGTH -- "
              "e.g., the length of the shortest input sequence; otherwise, "
              "when a probe of length p_l is mapped to a sequence of length "
              "s_l, then lcf_thres is treated as being min(LCF_THRES, p_l, "
              "s_l) so that a probe is able to 'cover' a sequence shorter "
              "than the probe and so that a probe shorter than lcf_thres "
              "is able to 'cover' a sequence"))
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
    parser.add_argument('--use-native-dict-when-finding-tolerant-coverage',
        dest="use_native_dict_when_finding_tolerant_coverage",
        action="store_true",
        help=("When finding probe coverage for blacklisting and "
              "identification (i.e., when using tolerant parameters), "
              "use a native Python dict as the kmer_probe_map across "
              "processes, rather than the primitives in SharedKmerProbeMap "
              "that are more suited to sharing across processes. Depending "
              "on the input (particularly if there are many candidate probes) "
              "this may result in substantial memory usage; but it may provide "
              "an improvement in runtime when there are relatively few "
              "candidate probes and a very large blacklisted input"))

    # Log levels and version
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
