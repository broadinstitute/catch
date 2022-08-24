#!/usr/bin/env python3
"""Design probes for genome capture.

This is the main executable of CATCH for probe design.
"""

import argparse
import importlib
import logging
import os
import random

from catch import coverage_analysis
from catch import probe
from catch.filter import adapter_filter
from catch.filter import duplicate_filter
from catch.filter import fasta_filter
from catch.filter import n_expansion_filter
from catch.filter import near_duplicate_filter
from catch.filter import polya_filter
from catch.filter import probe_designer
from catch.filter import reverse_complement_filter
from catch.filter import set_cover_filter
from catch.utils import cluster
from catch.utils import ncbi_neighbors
from catch.utils import seq_io, version, log

__author__ = 'Hayden Metsky <hayden@mit.edu>'


def main(args):
    logger = logging.getLogger(__name__)

    # Set NCBI API key
    if args.ncbi_api_key:
        ncbi_neighbors.ncbi_api_key = args.ncbi_api_key

    # Read the genomes from FASTA sequences
    genomes_grouped = []
    genomes_grouped_names = []
    for ds in args.dataset:
        if ds.startswith('collection:'):
            # Process a collection of datasets
            collection_name = ds[len('collection:'):]
            try:
                collection = importlib.import_module(
                    'catch.datasets.collections.' + collection_name)
            except ImportError:
                raise ValueError("Unknown dataset collection %s" %
                                 collection_name)
            for name, dataset in collection.import_all():
                genomes_grouped += [seq_io.read_dataset_genomes(dataset)]
                genomes_grouped_names += [name]
        elif ds.startswith('download:'):
            # Download a FASTA for an NCBI taxonomic ID
            taxid = ds[len('download:'):]
            if args.write_taxid_acc:
                taxid_fn = os.path.join(args.write_taxid_acc,
                        str(taxid) + '.txt')
            else:
                taxid_fn = None
            if '-' in taxid:
                taxid, segment = taxid.split('-')
            else:
                segment = None
            ds_fasta_tf = ncbi_neighbors.construct_fasta_for_taxid(taxid,
                    segment=segment, write_to=taxid_fn)
            genomes_grouped += [seq_io.read_genomes_from_fasta(ds_fasta_tf.name)]
            genomes_grouped_names += ['taxid:' + str(taxid)]
            ds_fasta_tf.close()
        elif os.path.isfile(ds):
            # Process a custom fasta file with sequences
            genomes_grouped += [seq_io.read_genomes_from_fasta(ds)]
            genomes_grouped_names += [os.path.basename(ds)]
        else:
            # Process an individual dataset
            try:
                dataset = importlib.import_module(
                            'catch.datasets.' + ds)
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
                        'catch.datasets.' + bg)
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

    # Setup and verify parameters related to k-mer length in probe map
    if args.kmer_probe_map_k:
        # Check that k is sufficiently small
        if args.kmer_probe_map_k > args.probe_length:
            raise Exception(("KMER_PROBE_MAP_K (%d) exceeds PROBE_LENGTH "
                             "(%d), which is not permitted") %
                            (args.kmer_probe_map_k, args.probe_length))

        # Use this value for the SetCoverFilter, AdapterFilter, and
        # the Analyzer
        kmer_probe_map_k_scf = args.kmer_probe_map_k
        kmer_probe_map_k_af = args.kmer_probe_map_k
        kmer_probe_map_k_analyzer = args.kmer_probe_map_k
    else:
        if args.probe_length <= 20:
            logger.warning(("PROBE_LENGTH (%d) is small; you may want to "
                            "consider setting --kmer-probe-map-k to be "
                            "small as well in order to be more sensitive "
                            "in mapping candidate probes to target sequence"),
                            args.probe_length)

        # Use a default k of 20 for the SetCoverFilter and AdapterFilter,
        # and 10 for the Analyzer since we would like to be more sensitive
        # (potentially at the cost of slower runtime) for the latter
        kmer_probe_map_k_scf = 20
        kmer_probe_map_k_af = 20
        kmer_probe_map_k_analyzer = 10

    # Set the maximum number of processes in multiprocessing pools
    if args.max_num_processes:
        probe.set_max_num_processes_for_probe_finding_pools(
            args.max_num_processes)
        cluster.set_max_num_processes_for_computing_distances(
            args.max_num_processes)

    # Raise exceptions or warn based on use of adapter arguments
    if args.add_adapters:
        if not (args.adapter_a or args.adapter_b):
            logger.warning(("Adapter sequences will be added, but default "
                            "sequences will be used; to provide adapter "
                            "sequences, use --adapter-a and --adapter-b"))
    else:
        if args.adapter_a or args.adapter_b:
            raise Exception(("Adapter sequences were provided with "
                "--adapter-a and --adapter-b, but --add-adapters is required "
                "to add adapter sequences onto the ends of probes"))

    # Do not allow both --small-seq-skip and --small-seq-min, since they
    # have different intentions
    if args.small_seq_skip is not None and args.small_seq_min is not None:
        raise Exception(("Both --small-seq-skip and --small-seq-min were "
            "specified, but both cannot be used together"))

    # Check arguments involving clustering
    if args.cluster_and_design_separately and args.identify:
        raise Exception(("Cannot use --cluster-and-design-separately with "
            "--identify, because clustering collapses genome groupings into "
            "one"))
    if args.cluster_from_fragments and not args.cluster_and_design_separately:
        raise Exception(("Cannot use --cluster-from-fragments without also "
            "setting --cluster-and-design-separately"))

    # Check for whether a custom hybridization function was provided
    if args.custom_hybridization_fn:
        custom_cover_range_fn = tuple(args.custom_hybridization_fn)
    else:
        custom_cover_range_fn = None
    if args.custom_hybridization_fn_tolerant:
        custom_cover_range_tolerant_fn = tuple(args.custom_hybridization_fn_tolerant)
    else:
        custom_cover_range_tolerant_fn = None

    # Setup the filters
    # The filters we use are, in order:
    filters = []

    # [Optional]
    # Fasta filter (ff) -- leave out candidate probes
    if args.filter_from_fasta:
        ff = fasta_filter.FastaFilter(args.filter_from_fasta,
                                      skip_reverse_complements=True)
        filters += [ff]

    # [Optional]
    # Poly(A) filter (paf) -- leave out probes with stretches of 'A' or 'T'
    if args.filter_polya:
        polya_length, polya_mismatches = args.filter_polya
        if polya_length > args.probe_length:
            logger.warning(("Length of poly(A) stretch to filter (%d) is "
                            "greater than PROBE_LENGTH (%d), which is usually "
                            "undesirable"), polya_length, args.probe_length)
        if polya_length < 10:
            logger.warning(("Length of poly(A) stretch to filter (%d) is "
                            "short, and may lead to many probes being "
                            "filtered"), polya_length)
        if polya_mismatches > 10:
            logger.warning(("Number of mismatches to tolerate when searching "
                            "for poly(A) stretches (%d) is high, and may "
                            "lead to many probes being filtered"),
                           polya_mismatches)
        paf = polya_filter.PolyAFilter(polya_length, polya_mismatches)
        filters += [paf]

    # Duplicate filter (df) -- condense all candidate probes that
    #     are identical down to one; this is not necessary for
    #     correctness, as the set cover filter achieves the same task
    #     implicitly, but it does significantly lower runtime by
    #     decreasing the input size to the set cover filter
    # Near duplicate filter (ndf) -- condense candidate probes that
    #     are near-duplicates down to one using locality-sensitive
    #     hashing; like the duplicate filter, this is not necessary
    #     but can significantly lower runtime and reduce memory usage
    #     (even more than the duplicate filter)
    if (args.filter_with_lsh_hamming is not None and
            args.filter_with_lsh_minhash is not None):
        raise Exception(("Cannot use both --filter-with-lsh-hamming "
            "and --filter-with-lsh-minhash"))
    if args.filter_with_lsh_hamming is not None:
        if args.filter_with_lsh_hamming > args.mismatches:
            logger.warning(("Setting FILTER_WITH_LSH_HAMMING (%d) to be greater "
                "than MISMATCHES (%d) may cause the probes to achieve less "
                "than the desired coverage"), args.filter_with_lsh_hamming,
                args.mismatches)
        ndf = near_duplicate_filter.NearDuplicateFilterWithHammingDistance(
            args.filter_with_lsh_hamming, args.probe_length)
        filters += [ndf]
    elif args.filter_with_lsh_minhash is not None:
        ndf = near_duplicate_filter.NearDuplicateFilterWithMinHash(
            args.filter_with_lsh_minhash)
        filters += [ndf]
    else:
        df = duplicate_filter.DuplicateFilter()
        filters += [df]

    # Set cover filter (scf) -- solve the problem by treating it as
    #     an instance of the set cover problem
    scf = set_cover_filter.SetCoverFilter(
        mismatches=args.mismatches,
        lcf_thres=args.lcf_thres,
        island_of_exact_match=args.island_of_exact_match,
        mismatches_tolerant=args.mismatches_tolerant,
        lcf_thres_tolerant=args.lcf_thres_tolerant,
        island_of_exact_match_tolerant=args.island_of_exact_match_tolerant,
        custom_cover_range_fn=custom_cover_range_fn,
        custom_cover_range_tolerant_fn=custom_cover_range_tolerant_fn,
        identify=args.identify,
        blacklisted_genomes=blacklisted_genomes_fasta,
        coverage=args.coverage,
        cover_extension=args.cover_extension,
        kmer_probe_map_k=kmer_probe_map_k_scf,
        kmer_probe_map_use_native_dict=args.use_native_dict_when_finding_tolerant_coverage)
    filters += [scf]

    # [Optional]
    # Adapter filter (af) -- add adapters to both the 5' and 3' ends
    #    of each probe
    if args.add_adapters:
        # Set default adapter sequences, if not provided
        if args.adapter_a:
            adapter_a = tuple(args.adapter_a)
        else:
            adapter_a = ('ATACGCCATGCTGGGTCTCC', 'CGTACTTGGGAGTCGGCCAT')
        if args.adapter_b:
            adapter_b = tuple(args.adapter_b)
        else:
            adapter_b = ('AGGCCCTGGCTGCTGATATG', 'GACCTTTTGGGACAGCGGTG')

        af = adapter_filter.AdapterFilter(adapter_a,
                                          adapter_b,
                                          mismatches=args.mismatches,
                                          lcf_thres=args.lcf_thres,
                                          island_of_exact_match=\
                                            args.island_of_exact_match,
                                          custom_cover_range_fn=\
                                            custom_cover_range_fn,
                                          kmer_probe_map_k=kmer_probe_map_k_af)
        filters += [af]

    # [Optional]
    # N expansion filter (nef) -- expand Ns in probe sequences
    # to avoid ambiguity
    if args.expand_n is not None:
        nef = n_expansion_filter.NExpansionFilter(
            limit_n_expansion_randomly=args.expand_n)
        filters += [nef]

    # [Optional]
    # Reverse complement (rc) -- add the reverse complement of each
    #    probe that remains
    if args.add_reverse_complements:
        rc = reverse_complement_filter.ReverseComplementFilter()
        filters += [rc]

    # If requested, don't apply the set cover filter
    if args.skip_set_cover:
        filter_before_scf = filters[filters.index(scf) - 1]
        filters.remove(scf)

    # Define parameters for clustering sequences
    if args.cluster_and_design_separately:
        cluster_threshold = args.cluster_and_design_separately
        if args.skip_set_cover:
            cluster_merge_after = filter_before_scf
        else:
            cluster_merge_after = scf
        cluster_method = args.cluster_and_design_separately_method
        cluster_fragment_length = args.cluster_from_fragments
    else:
        cluster_threshold = None
        cluster_merge_after = None
        cluster_method = None
        cluster_fragment_length = None

    # Design the probes
    pb = probe_designer.ProbeDesigner(genomes_grouped, filters,
                                      probe_length=args.probe_length,
                                      probe_stride=args.probe_stride,
                                      allow_small_seqs=args.small_seq_min,
                                      seq_length_to_skip=args.small_seq_skip,
                                      cluster_threshold=cluster_threshold,
                                      cluster_merge_after=cluster_merge_after,
                                      cluster_method=cluster_method,
                                      cluster_fragment_length=cluster_fragment_length)
    pb.design()

    # Write the final probes to the file args.output_probes
    seq_io.write_probe_fasta(pb.final_probes, args.output_probes)

    if (args.print_analysis or args.write_analysis_to_tsv or
            args.write_sliding_window_coverage or
            args.write_probe_map_counts_to_tsv):
        analyzer = coverage_analysis.Analyzer(
            pb.final_probes,
            args.mismatches,
            args.lcf_thres,
            genomes_grouped,
            genomes_grouped_names,
            island_of_exact_match=args.island_of_exact_match,
            custom_cover_range_fn=custom_cover_range_fn,
            cover_extension=args.cover_extension,
            kmer_probe_map_k=kmer_probe_map_k_analyzer,
            rc_too=args.add_reverse_complements)
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
    else:
        # Just print the number of probes
        print(len(pb.final_probes))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    # Input data
    parser.add_argument('dataset',
        nargs='+',
        help=("One or more target datasets (e.g., one per species). Each "
              "dataset can be specified in one of multiple ways. (a) If "
              "dataset is in the format 'download:TAXID', then CATCH downloads "
              "from NCBI all whole genomes for the NCBI taxonomy with id "
              "TAXID, and uses these sequences as input. (b) If dataset is "
              "a path to a FASTA file, then its sequences are read and used "
              "as input. (c) Otherwise, it is assumed that this is a label "
              "for a dataset included in this package (e.g., 'zika'). If "
              "the label starts with 'collection:' (e.g., 'collection:viruses"
              "_with_human_host'), then this reads from an available "
              "collection of datasets. For segmented viruses, the format "
              "for NCBI downloads can also be 'download:TAXID-SEGMENT'."))

    # Outputting probes
    parser.add_argument('-o', '--output-probes',
        required=True,
        help=("The file to which all final probes should be "
              "written; they are written in FASTA format"))

    # Outputting downloaed data
    parser.add_argument('--write-taxid-acc',
        help=("If 'download:' labels are used in datasets, write downloaded "
              "accessions to a file in this directory. Accessions are written "
              "to WRITE_TAXID_ACC/TAXID.txt"))

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

    # Custom function (dynamically loaded) to determine probe hybridization
    # When set, this makes values of the above arguments (--mismatches,
    # --lcf-thres, and --island-of-exact-match) meaningless
    parser.add_argument('--custom-hybridization-fn',
        nargs=2,
        help=("(Optional) Args: <PATH> <FUNC>; PATH is a path to a Python "
              "module (.py file) and FUNC is a string giving the name of "
              "a function in that module. FUNC provides a custom model of "
              "hybridization between a probe and target sequence to use in "
              "the probe set design. If this is set, the arguments "
              "--mismatches, --lcf-thres, and --island-of-exact-match are "
              "not used because these are meant for the default model of "
              "hybridization. The function FUNC in PATH is dynamically "
              "loaded to use when determining whether a probe hybridizes to "
              "a target sequence (and, if so, what portion). FUNC must "
              "accept the following arguments in order, though it "
              "may choose to ignore some values: (1) array giving sequence "
              "of a probe; (2) str giving subsequence of target sequence to "
              "which the probe may hybridize, of the same length as the "
              "given probe sequence; (3) int giving the position in the "
              "probe (equivalently, the target subsequence) of the start "
              "of a k-mer around which the probe and target subsequence "
              "are anchored (the probe and target subsequence are aligned "
              "using this k-mer as an anchor); (4) int giving the end "
              "position (exclusive) of the anchor k-mer; (5) int giving the "
              "full length of the probe (the probe provided in (1) may be "
              "cutoff on an end if it extends further than where the "
              "target sequence ends); (6) int giving the full length of the "
              "target sequence of which the subsequence in (2) is part. "
              "FUNC must return None if it deems that the probe does not "
              "hybridize to the target subsequence; otherwise, it must "
              "return a tuple (start, end) where start is an int giving "
              "the start position in the probe (equivalently, in the "
              "target subsequence) at which the probe will hybridize to "
              "the target subsequence, and end is an int (exclusive) giving "
              "the end position of the hybridization."))

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
    parser.add_argument('--custom-hybridization-fn-tolerant',
        nargs=2,
        help=("(Optional) A more tolerant model than the one "
              "implemented in custom_hybridization_fn. This should capture "
              "more possible hybridizations (i.e., be more sensitive) "
              "when designing probes for identification or when genomes "
              "are blacklisted. See --custom-hybridization-fn for details "
              "of how this function should be implemented and provided."))

    # Outputting coverage analyses
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
    parser.add_argument('--write-probe-map-counts-to-tsv',
        help=("(Optional) The file to which to write a TSV-formatted list of "
              "the number of sequences each probe maps to. This explicitly "
              "does not count reverse complements."))

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
    parser.add_argument('--add-adapters',
        dest="add_adapters",
        action="store_true",
        help=("Add adapters to the ends of probes; to specify adapter "
              "sequences, use --adapter-a and --adapter-b"))
    parser.add_argument('--adapter-a',
        nargs=2,
        help=("(Optional) Args: <X> <Y>; Custom A adapter to use; two ordered "
              "where X is the A adapter sequence to place on the 5' end of "
              "a probe and Y is the A adapter sequence to place on the 3' "
              "end of a probe"))
    parser.add_argument('--adapter-b',
        nargs=2,
        help=("(Optional) Args: <X> <Y>; Custom B adapter to use; two ordered "
              "where X is the B adapter sequence to place on the 5' end of "
              "a probe and Y is the B adapter sequence to place on the 3' "
              "end of a probe"))

    # Filtering poly(A) sequence from probes
    parser.add_argument('--filter-polya',
        nargs=2,
        type=int,
        help=("(Optional) Args: <X> <Y> (integers); do not output any probe "
              "that contains a stretch of X or more 'A' bases, tolerating "
              "up to Y mismatches (and likewise for 'T' bases)"))

    # Adjusting probe output
    parser.add_argument('--add-reverse-complements',
        dest="add_reverse_complements",
        action="store_true",
        help=("Add to the output the reverse complement of each probe"))
    parser.add_argument('--expand-n',
        nargs='?',
        type=int,
        default=None,
        const=3,
        help=("Expand each probe so that 'N' bases are replaced by real "
              "bases; for example, the probe 'ANA' would be replaced "
              "with the probes 'AAA', 'ATA', 'ACA', and 'AGA'; this is "
              "done combinatorially across all 'N' bases in a probe, and "
              "thus the number of new probes grows exponentially with the "
              "number of 'N' bases in a probe. If followed by a command- "
              "line argument (INT), this only expands at most INT randomly "
              "selected N bases, and the rest are replaced with random "
              "unambiguous bases (default INT is 3)."))

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

    # Clustering input sequences
    def check_cluster_and_design_separately(val):
        fval = float(val)
        if fval > 0 and fval <= 0.5:
            # a float in (0,0.5]
            return fval
        else:
            raise argparse.ArgumentTypeError(("%s is an invalid average "
                                              "nucleotide dissimilarity") % val)
    parser.add_argument('--cluster-and-design-separately',
        type=check_cluster_and_design_separately,
        help=("(Optional) If set, cluster all input sequences using their "
              "MinHash signatures, design probes separately on each cluster, "
              "and combine the resulting probes. This can significantly lower "
              "runtime and memory usage, but may lead to a suboptimal "
              "solution. The value CLUSTER_AND_DESIGN_SEPARATELY gives the "
              "distance threshold for determining clusters in terms of "
              "average nucleotide dissimilarity (1-ANI, where ANI is "
              "average nucleotide identity; see --cluster-and-design-"
              "separately-method for details); higher values "
              "result in fewer clusters, and thus longer runtime. Values "
              "must be in (0,0.5], and generally should be around 0.1 or "
              "0.2. When used, this creates a separate genome for each "
              "input sequence -- it collapses all sequences, across both "
              "groups and genomes, into one list of sequences in one group. "
              "Therefore, genomes will not be grouped as specified in the "
              "input and sequences will not be grouped by genome, and "
              "differential identification is not supported"))
    parser.add_argument('--cluster-and-design-separately-method',
        choices=['simple', 'hierarchical'], default='simple',
        help=("(Optional) Method for clustering input sequences, which is "
              "only used if --cluster-and-design-separately is set. If "
              "'simple', clusters are connected components of a graph in "
              "which each sequence is a vertex and two sequences are adjacent "
              "if their estimated nucleotide dissimilarity is within "
              "the value CLUSTER_AND_DESIGN_SEPARATELY. If 'hierarchical', "
              "clusters are determined by agglomerative hierarchical "
              "clustering and the the value CLUSTER_AND_DESIGN_SEPARATELY "
              "is the inter-cluster distance threshold to merge clusters."))
    parser.add_argument('--cluster-from-fragments',
        type=int,
        help=("(Optional) If set, break all sequences into sequences of "
              "length CLUSTER_FROM_FRAGMENTS nt, and cluster these fragments. "
              "This can be useful for improving runtime on input with "
              "especially large genomes, in which probes for different "
              "fragments can be designed separately. Values should generally "
              "be around 10,000. For this to be used, "
              "--cluster-and-design-separately must also be set."))

    # Filter candidate probes with LSH
    parser.add_argument('--filter-with-lsh-hamming',
        type=int,
        help=("(Optional) If set, filter candidate probes for near-"
              "duplicates using LSH with a family of hash functions that "
              "works with Hamming distance. FILTER_WITH_LSH_HAMMING gives "
              "the maximum Hamming distance at which to call near-"
              "duplicates; it should be commensurate with (but not greater "
              "than) MISMATCHES. Using this may significantly improve "
              "runtime and reduce memory usage by reducing the number of "
              "candidate probes to consider, but may lead to a slightly "
              "sub-optimal solution. It may also, particularly with "
              "relatively high values of FILTER_WITH_LSH_HAMMING, cause "
              "coverage obtained for each genome to be slightly less than "
              "the desired coverage (COVERAGE) when that desired coverage "
              "is the complete genome; it is recommended to also use "
              "--print-analysis or --write-analysis-to-tsv with this "
              "to see the coverage that is obtained."))
    def check_filter_with_lsh_minhash(val):
        fval = float(val)
        if fval >= 0.0 and fval <= 1.0:
            # a float in [0,1]
            return fval
        else:
            raise argparse.ArgumentTypeError(("%s is an invalid Jaccard "
                                              "distance") % val)
    parser.add_argument('--filter-with-lsh-minhash',
        type=check_filter_with_lsh_minhash,
        help=("(Optional) If set, filter candidate probes for near-"
              "duplicates using LSH with a MinHash family. "
              "FILTER_WITH_LSH_MINHASH gives the maximum Jaccard distance "
              "(1 minus Jaccard similarity) at which to call near-duplicates; "
              "the Jaccard similarity is calculated by treating each probe "
              "as a set of overlapping 10-mers. Its value should be "
              "commensurate with parameter values determining whether a probe "
              "hybridizes to a target sequence, but this can be difficult "
              "to measure compared to the input for --filter-with-lsh-hamming. "
              "However, this allows more sensitivity in near-duplicate "
              "detection than --filter-with-lsh-hamming (e.g., if near-"
              "duplicates should involve probes shifted relative to each "
              "other). The same caveats mentioned in help for "
              "--filter-with-lsh-hamming also apply here. Values of "
              "FILTER_WITH_LSH_MINHASH above ~0.7 may start to require "
              "significant memory and runtime for near-duplicate detection."))

    # Miscellaneous technical adjustments
    parser.add_argument('--small-seq-skip',
        type=int,
        help=("(Optional) Do not create candidate probes from sequences "
              "whose length is <= SMALL_SEQ_SKIP. If set to (PROBE_LENGTH - "
              "1), this avoids the error raised when sequences are less "
              "than the probe length"))
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
    parser.add_argument('--kmer-probe-map-k',
        type=int,
        help=("(Optional) Use this value (KMER_PROBE_LENGTH_K) as the "
              "k-mer length when constructing a map of k-mers to the probes "
              "that contain these k-mers. This map is used when mapping "
              "candidate probes to target sequences and the k-mers serve "
              "as seeds for calculating whether a candidate probe 'covers' "
              "a subsequence. The value should be sufficiently less than "
              "PROBE_LENGTH so that it can find mappings even when the "
              "candidate probe and target sequence are divergent. In "
              "particular, CATCH will try to find a value k >= "
              "KMER_PROBE_LENGTH_K (by default, >=20) such that k divides "
              "PROBE_LENGTH and k < PROBE_LENGTH / MISMATCHES (if "
              "MISMATCHES=0, then k=PROBE_LENGTH). It will then use this "
              "k as the k-mer length in mappings; if no such k exists, it "
              "will use a randomized approach with KMER_PROBE_LENGTH_K as "
              "the k-mer length. If --custom-hybridization-fn is set, "
              "it will always use the randomized approach with "
              "KMER_PROBE_LENGTH_K (by default, 20) as the k-mer length."))
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
    parser.add_argument('--ncbi-api-key',
        help=("API key to use for NCBI e-utils. Using this increases the "
              "limit on requests/second and may prevent an IP address "
              "from being block due to too many requests"))

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
