# CATCH &nbsp;&middot;&nbsp; [![Build Status](https://travis-ci.org/broadinstitute/catch.svg?branch=master)](https://travis-ci.org/broadinstitute/catch) [![Coverage Status](https://coveralls.io/repos/broadinstitute/catch/badge.svg?branch=master)](https://coveralls.io/github/broadinstitute/catch?branch=master) [![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg)](https://github.com/broadinstitute/catch/pulls) [![MIT License](https://img.shields.io/badge/license-MIT-blue.svg)](./LICENSE)
#### Compact Aggregation of Targets for Comprehensive Hybridization

CATCH is a Python package for designing probe sets to use in hybrid capture experiments.

* **Comprehensive coverage**: CATCH accepts any collection of unaligned sequences &mdash; typically whole genomes of all known genetic diversity of one or more microbial species.
It designs oligo sequences that guarantee coverage of this diversity, enabling rapid design of exhaustive probe sets for customizable targets.
* **Compact designs**: CATCH can design with a specified constraint on the number of oligos (e.g., array size).
It searches a space of probe sets, which may pool many species, to find an optimal design.
This allows its designs to scale well with known genetic diversity, and also supports cost-effective experiments.
* **Flexibility**: CATCH supports applications beyond whole genome enrichment, such as differential identification of species.
It allows blacklisting sequence from the design (e.g., background in microbial enrichment), supports customized models of hybridization, enables weighting the sensitivity for different species, and more.
<br/>

### Table of contents

* [Setting up CATCH](#setting-up-catch)
  * [Python dependencies](#python-dependencies)
  * [Downloading and installing](#downloading-and-installing)
  * [Downloading viral sequence data](#downloading-viral-sequence-data)
  * [Testing](#testing)
* [Using CATCH](#using-catch)
  * [Designing with one choice of parameteters](#designing-with-one-choice-of-parameters-designpy)
  * [Pooling across many runs](#pooling-across-many-runs-poolpy)
* [Examples](#examples)
  * [Example of running design.py](#example-of-running-designpy)
  * [Example of running pool.py](#example-of-running-poolpy)
* [Contributing](#contributing)
* [License](#license)
<br/>

## Setting up CATCH

### Python dependencies

CATCH requires:
* [Python](https://www.python.org) &gt;= 3.5
* [NumPy](http://www.numpy.org) &gt;= 1.9.0
* [SciPy](https://www.scipy.org) &gt;= 1.0.0

Installing CATCH with `pip`, as described below, will install NumPy and SciPy if they are not already installed.

### Downloading and installing

The easiest way to setup CATCH is with `pip`:

```bash
git clone git@github.com:broadinstitute/catch.git
cd catch
pip install --user -e .
```

### Downloading viral sequence data

We distribute viral sequence data with CATCH, which can be used as input to probe design.
We use [Git LFS](https://git-lfs.github.com) to version and store this data.
If you wish to use this data, you'll need to install [Git LFS](https://git-lfs.github.com).
After installing it, you can download the viral sequence data by running:

```bash
git lfs install
git lfs pull
```

from inside the `catch` project directory.

### Testing

CATCH uses Python's `unittest` framework.
Some of these tests require you to have [downloaded](#downloading-viral-sequence-data) viral sequence data.
To execute all unit tests, run:

```bash
python -m unittest discover
```

## Using CATCH

### Designing with one choice of parameters ([`design.py`](./bin/design.py))

The main program to design probes is [`design.py`](./bin/design.py).
To see details on all the arguments that the program accepts, run:

```bash
design.py --help
```

[`design.py`](./bin/design.py) requires one or more `dataset`s that specify input sequence data to target:

```bash
design.py [dataset] [dataset ...]
```

Each `dataset` can be a path to a FASTA file. If you [downloaded](#downloading-viral-sequence-data) viral sequence data, it can also simply be a label for one of [350+ viral datasets](./catch/datasets/README.md) (e.g., `hiv1` or `zika`) distributed as part of this package.
Each of these datasets includes all available whole genomes (genome neighbors) in [NCBI's viral genome data](https://www.ncbi.nlm.nih.gov/genome/viruses/) that have human as a host, for one or more species, as of Sep. 2017.

Below is a summary of some useful arguments to `design.py`:

* `-pl PROBE_LENGTH`/`-ps PROBE_STRIDE`: Design probes to be PROBE_LENGTH nt long, and generate candidate probes using a stride of PROBE_STRIDE nt.
(Default: 100 and 50.)
* `-m MISMATCHES`: Tolerate up to MISMATCHES mismatches when determining whether a probe covers a target sequence.
(Also, see `-l/--lcf-thres` and `--island-of-exact-match` for adjusting hybridization criteria.)
Higher values lead to fewer probes.
(Default: 0.)
* `-c COVERAGE`: Guarantee that at least COVERAGE of each target genome is captured by probes, where COVERAGE is either a fraction of a genome or a number of nucleotides.
Higher values lead to more probes.
(Default: 1.0 &mdash; i.e., whole genome.)
* `-e COVER_EXTENSION`: Assume that a probe will capture both the region of the sequence to which it hybridizes, as well as COVER_EXTENSION nt on each side of that.
Higher values lead to fewer probes.
(Default: 0.)
* `--identify`: Design probes to perform differential identification.
This is typically used with small values of COVERAGE and &gt;1 specified `dataset`s.
Probes are designed such that each `dataset` should be captured by probes that are unlikely to hybridize to other `dataset`s.
* `--blacklist-genomes dataset [dataset ...]`: Design probes to be unlikely to hybridize to any of these datasets.
(Also, see `-mt/--mismatches-tolerant`, `-lt/--lcf-thres-tolerant`, and `--island-of-exact-match-tolerant` for this and for `--identify`.)
* `--add-adapters`: Add PCR adapters to the ends of each probe sequence.
This selects adapters to add to probe sequences so as to minimize overlap among probes that share an adapter, allowing probes with the same adapter to be amplified together.
(See `--adapter-a` and `--adapter-b` too.)
* `--custom-hybridization-fn PATH FN`: Specify a function, for CATCH to dynamically load, that implements a custom model of hybridization between a probe and target sequence.
See `design.py --help` for details on the expected input and output of this function.
If not set, CATCH uses its default model of hybridization based on `-m/--mismatches`, `-l/--lcf-thres`, and `--island-of-exact-match`.
(Relatedly, see `--custom-hybridization-fn-tolerant`.)
* `--filter-with-lsh-hamming FILTER_WITH_LSH_HAMMING`/`--filter-with-lsh-minhash FILTER_WITH_LSH_MINHASH`: Use locality-sensitive hashing to reduce the space of candidate probes.
This can significantly improve runtime and memory requirements when the input is especially large and diverse.
See `design.py --help` for details on using these options and downsides.
* `-o OUTPUT`: Write probe sequences in FASTA format to OUTPUT.

### Pooling across many runs ([`pool.py`](./bin/pool.py))

While [`design.py`](./bin/design.py) requires particular choices of parameter values, [`pool.py`](./bin/pool.py) is a program to find optimal hybridization parameters that can vary across many input, under a specified limit on the total number of probes (e.g., synthesis array size).
It does this by searching over a space of probe sets to solve a constrained optimization problem.
To see details on all the arguments that the program accepts, run:

```bash
pool.py --help
```

You need to run [`design.py`](./bin/design.py) on each dataset over a grid of parameters values that spans a reasonable domain.
Then, create a table that provides a probe count for each dataset and choice of parameters (TSV, in a format like [this](./catch/pool/tests/input/num-probes.V-WAfr.201506.tsv)).
Now, you can use this table as input:

```bash
pool.py INPUT_TSV TARGET_PROBE_COUNT OUTPUT_TSV
```
where INPUT_TSV is a path to the table described above, TARGET_PROBE_COUNT is a constraint on the number of probes to allow in the pool, and OUTPUT_TSV is a path to a file to which the program will write the optimal parameter values.

Below are two arguments that generalize the search:
* `--loss-coeffs COEFF [COEFF ...]`: Specify coefficients on parameters in the objective function.
This allows you to adjust how conservative each parameter is treated relative to others.
(Default: 1 for mismatches and 1/100 for cover extension.)
* `--dataset-weights WEIGHTS_TSV`: Assign a weight for each dataset to use in the objective function, where WEIGHTS_TSV is a path to a table that provides a weight for each dataset.
This allows you to seek that probes in the pooled design be more sensitive for some taxa than others.
(Default: 1 for all datasets.)

Each run of [`pool.py`](./bin/pool.py) may yield a different output based on the (random) initial guess.
We recommend running this multiple times and selecting the output that has the smallest loss, which is written to standard output at the end of the program.

## Examples

### Example of running [`design.py`](./bin/design.py)

Below is an example of designing probes to target a single taxon.

```bash
design.py zika -pl 75 -m 2 -l 60 -e 50 -o zika-probes.fasta
```

This will design probes that:
* target whole genomes of Zika virus (`zika`)
* are 75 nt long (`-pl 75`)
* capture the entirety of each genome under a model that a probe hybridizes to a region if the longest common substring, up to 2 mismatches (`-m 2`), between a probe and target is at least 60 nt (`-l 60`)
* assume 50 nt on each side of the hybridization is captured as well (`-e 50`)

and will save them to `zika-probes.fasta`.

It yields about 500 probes.
Note that the `zika` dataset distributed with CATCH contains 379 genomes, but the input can also be a path to any custom FASTA file.

### Example of running [`pool.py`](./bin/pool.py)

[Here](./catch/pool/tests/input/num-probes.V-WAfr.201506.tsv) is a table listing probe counts used in the design of the [V-WAfr probe set](./probe-designs).
It provides counts for each dataset and combination of two parameters (mismatches and cover extension) that were varied in the design.
Below is an example of designing that probe set using this table as input.

```bash
pool.py num-probes.V-WAfr.201506.tsv 90000 params.V-Wafr.201506.tsv --round-params 1 10
```

This will search for parameters that yield at most 90,000 probes across the datasets, and will output those to `params.V-Wafr.201506.tsv`.
Because the search is over a continuous space, here we use `--round-params 1 10` to set each value of the mismatches parameter to an integer and each value of the cover extension parameter to a multiple of 10 while still meeting the constraint on probe count.
The pooled design yields about 89,950 probes, depending on the initial guess.

## Contributing

We welcome contributions to CATCH. This can be in the form of an [issue](https://github.com/broadinstitute/catch/issues)  or [pull request](https://github.com/broadinstitute/catch/pulls). If you have questions, please create an [issue](https://github.com/broadinstitute/catch/issues) or email **Hayden Metsky** &lt;hayden@mit.edu&gt;.

## License

CATCH is licensed under the terms of the [MIT license](./LICENSE).

