"""Utilities for reading probe count tables and writing parameter values.
"""

import logging

__author__ = 'Hayden Metsky <hayden@mit.edu>'

logger = logging.getLogger(__name__)


def read_table_of_probe_counts(fn):
    """Read a table listing probe counts for combinations of parameters.

    In the input table, the first column must give the dataset (header:
    'dataset') and the last column must give the number of probes (header:
    'num_probes'). The intermediary columns give choices of parameters.

    Args:
        fn: path to tab-separated file to read, in the format described
            above

    Returns:
        tuple (x, y) where:
            - x is a tuple listing the names of parameters included in
              the table (according to the header)
            - y is a dict of the format:
                { dataset: { param_values: probe count } }
              where param_values is a tuple giving a choice of parameter
              values in the same order of parameters in x
    """
    d = {}
    with open(fn) as f:
        for i, line in enumerate(f):
            ls = line.rstrip().split('\t')
            if i == 0:
                # Read the header and verify it meets the required
                # order
                header = ls
                if header[0] != "dataset":
                    raise Exception(("First column in probe count table "
                        "must be 'dataset'"))
                if header[-1] != "num_probes":
                    raise Exception(("Last column in probe count table "
                        "must be 'num_probes'"))
                param_names = tuple(ls[1:-1])
                continue

            assert len(ls) == 2 + len(param_names)
            dataset = ls[0]
            num_probes = int(ls[-1])
            param_values = tuple([float(x) for x in ls[1:-1]])

            if dataset not in d:
                d[dataset] = {}
            if param_values in d[dataset]:
                raise Exception(("The same combination of dataset and "
                    "parameters is listed more than once in the probe "
                    "count table"))
            d[dataset][param_values] = num_probes
    return (param_names, d)


def read_table_of_dataset_weights(fn, datasets_to_check=None):
    """Read a table listing a weight for each dataset.

    In the input table, the first column must give the dataset (header:
    'dataset') and the second column must give the weight (header:
    'weight').

    Args:
        fn: path to tab-separated file to read, in the format described
            above
        datasets_to_check: if set, this verifies that each dataset
            present in this collection has a weight in the input
            table, and raises an Exception if one does not

    Returns:
        dict {dataset: weight}
    """
    d = {}
    with open(fn) as f:
        for i, line in enumerate(f):
            ls = line.rstrip().split('\t')
            if i == 0:
                # Read the header and verify it meets the required
                # order
                header = ls
                if header[0] != "dataset":
                    raise Exception(("First column in dataset weights table "
                        "must be 'dataset'"))
                if header[1] != "weight":
                    raise Exception(("Second column in dataset weights table "
                        "must be 'weight'"))
                if len(header) > 2:
                    raise Exception(("There can only be two columns in "
                        "the dataset weights table"))
                continue

            assert len(ls) == 2
            dataset = ls[0]
            weight = float(ls[1])

            if dataset in d:
                raise Exception(("The same dataset (%s) appears on more "
                    "than one row in the dataset weights table") % dataset)

            d[dataset] = weight

    if datasets_to_check is not None:
        for dataset in datasets_to_check:
            if dataset not in d:
                raise Exception(("dataset %s needs a weight, but one is "
                    "not given in the dataset weights table") % dataset)

    return d


def write_param_values_across_datasets(param_names, param_vals, out_tsv,
                                       type='int'):
    """Write parameter values for each dataset.

    This outputs a tab-separated file in which the first row is a
    header, the first column gives the name of a dataset, and the
    remaining columns give parameter values assigned to that dataset.

    Args:
        param_names: tuple of parameter names
        param_vals: dict {dataset: p} where p is a tuple consisting of
            parameter values, corresponding to the order in param_names
        out_tsv: path to output TSV file
        type: 'int' or 'float'; type to use when writing parameter
            values
    """
    header = '\t'.join(['dataset'] + list(param_names))
    lines = [header]
    for dataset in sorted(param_vals.keys()):
        vals = param_vals[dataset]
        if type == 'float':
            line = '\t'.join([dataset] + ['%f' % p for p in vals])
        elif type == 'int':
            line = '\t'.join([dataset] + ['%d' % p for p in vals])
        else:
            raise ValueError("Unknown type %s", type)
        lines += [line]

    with open(out_tsv, 'w') as f:
        for line in lines:
            f.write(line + '\n')

