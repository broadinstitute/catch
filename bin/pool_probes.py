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

    # Verify that the only parameters are, in order: 'mismatches' and
    # 'cover_extension'
    if param_names != ('mismatches', 'cover_extension'):
        raise Exception(("For a standard search, the only parameters in the "
                         "input table must be, in order: 'mismatches' and "
                         "cover_extension"))

    # Perform a standard search for optimal values of mismatches and
    # cover extension
    ss = param_search.standard_search(probe_counts, args.target_probe_count)
    opt_params, opt_params_count, opt_params_loss = ss

    # Write a table of the optimal parameter values
    pool_probes_io.write_param_values_across_datasets(param_names, opt_params,
        args.param_vals_tsv)

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
