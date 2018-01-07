#!/usr/bin/env python3
"""Pool probes across datasets by searching for optimal parameters."""

import argparse
import logging

from hybseldesign.pool import param_search
from hybseldesign.utils import log, version
from hybseldesign.utils import pool_probes_io

__author__ = 'Hayden Metsky <hayden@mit.edu>'

logger = logging.getLogger(__name__)


def main(args):
    # Read the table of probe counts
    param_names, probe_counts = pool_probes_io.read_table_of_probe_counts(
        args.probe_count_tsv)

    if args.use_nd:
        # This does not round parameters after searching over the
        # dimensional space
        if args.round_params:
            raise Exception(("The arguments '--use-nd' and '--round-params' "
                "cannot both be used; this does not round parameters "
                "after searching over a space with n > 2"))

        # Perform a higher dimensional search for optimal values of
        # the parameters
        s_results = param_search.higher_dimensional_search(
            param_names, probe_counts, args.target_probe_count)
        write_type = 'float'
    else:
        # For the standard search, the only parameters must be (in order):
        #' mismatches' and 'cover_extension'. Verify this.
        if param_names != ('mismatches', 'cover_extension'):
            raise Exception(("For a standard search, the only parameters "
                "in the input table must be, in order: 'mismatches' and "
                "'cover_extension'. Consider using the '--use-nd' argument "
                "to search over additional parameters."))

        if args.round_params:
            mismatches_round, cover_extension_round = args.round_params
        else:
            mismatches_round, cover_extension_round = 1, 1

        # Perform a standard search for optimal values of mismatches and
        # cover extension
        s_results = param_search.standard_search(
            probe_counts, args.target_probe_count,
            mismatches_round=mismatches_round,
            cover_extension_round=cover_extension_round)
        write_type = 'int'

    opt_params, opt_params_count, opt_params_loss = s_results

    # Write a table of the optimal parameter values
    pool_probes_io.write_param_values_across_datasets(param_names, opt_params,
        args.param_vals_tsv, type=write_type)

    # Print the total number of probes and loss
    print("Number of probes: %d" % opt_params_count)
    print("Loss: %f" % opt_params_loss)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('probe_count_tsv',
        help=("Path to TSV file that contains probe counts for each "
              "dataset and combination of parameters; the first row "
              "must be a header, the first column must give a "
              "dataset ('dataset'), the last column must list a "
              "number of probes ('num_probes'), and the intermediary "
              "columns give parameter values"))
    parser.add_argument('target_probe_count', type=int,
        help=("Constraint on the total number of probes in the design; "
              "generally, parameters will be selected such that the "
              "number of probes, when pooled across datasets, is "
              "just below this number"))
    parser.add_argument('param_vals_tsv',
        help=("Path to TSV file in which to output optimal parameter "
              "values"))
    parser.add_argument('--round-params', nargs=2, type=int,
        help=("<m> <e>; round mismatches parameter to the nearest "
              "multiple of m and cover_extension parameter to the "
              "nearest multiple of e"))
    parser.add_argument('--use-nd', action='store_true',
        help=("Use the higher dimensional (n > 2) interpolation and "
              "search functions for optimizing parameters. This is "
              "required if the input table of probe counts contains "
              "more than 2 parameters. This does not round parameters "
              "to integers or to be placed on a grid -- i.e., they "
              "will be output as fractional values (from which probe "
              "counts were interpolated)."))
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
