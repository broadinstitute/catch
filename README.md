# CATCH &nbsp;&middot;&nbsp; [![Build Status](https://app.travis-ci.com/broadinstitute/catch.svg?branch=master)](https://app.travis-ci.com/broadinstitute/catch) [![Coverage Status](https://coveralls.io/repos/broadinstitute/catch/badge.svg?branch=master)](https://coveralls.io/github/broadinstitute/catch?branch=master) [![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg)](https://github.com/broadinstitute/catch/pulls) [![MIT License](https://img.shields.io/badge/license-MIT-blue.svg)](./LICENSE)
#### Compact Aggregation of Targets for Comprehensive Hybridization

CATCH is a Python package for designing probe sets to use for nucleic acid capture of diverse sequence.

* **Comprehensive coverage**: CATCH accepts any collection of unaligned sequences &mdash; typically whole genomes of all known genetic diversity of one or more microbial species.
It designs oligo sequences that guarantee coverage of this diversity, enabling rapid design of exhaustive probe sets for customizable targets.
* **Compact designs**: CATCH can design with a specified constraint on the number of oligos (e.g., array size).
It searches a space of probe sets, which may pool many species, to find an optimal design.
This allows its designs to scale well with known genetic diversity, and also supports cost-effective applications.
* **Flexibility**: CATCH supports applications beyond whole genome enrichment, such as differential identification of species.
It allows avoiding sequence from the design (e.g., background in microbial enrichment), supports customized models of hybridization, enables weighting the sensitivity for different species, and more.
<br/>

### Table of contents

* [Setting up CATCH](#setting-up-catch)
  * [Python dependencies](#python-dependencies)
  * [Setting up a conda environment](#setting-up-a-conda-environment)
  * [Downloading and installing](#downloading-and-installing)
  * [Downloading viral sequence data](#downloading-viral-sequence-data)
  * [Testing](#testing)
  * [Alternative approach: installing with conda](#alternative-approach-installing-with-conda)
* [Using CATCH](#using-catch)
  * [Cautious defaults and one setting of hybridization criteria](#option-1-cautious-defaults-and-one-setting-of-hybridization-criteria)
  * [Pragmatic defaults for large, diverse input and one setting of hybridization criteria](#option-2-pragmatic-defaults-for-large-diverse-input-and-one-setting-of-hybridization-criteria)
  * [Variable hybridization criteria across many taxa](#option-3-variable-hybridization-criteria-across-many-taxa)
* [Examples](#examples)
  * [Example of running design.py](#example-of-running-designpy)
  * [Example of running design_large.py](#example-of-running-designlargepy)
  * [Example of running pool.py](#example-of-running-poolpy)
* [Contributing](#contributing)
* [Citation](#citation)
* [License](#license)
<br/>

## Setting up CATCH

### Python dependencies

CATCH requires:
* [Python](https://www.python.org) &gt;= 3.8
* [NumPy](http://www.numpy.org) &gt;= 1.22
* [SciPy](https://www.scipy.org) &gt;= 1.8.0

CATCH may also work with older versions of Python, NumPy, and SciPy, but is only tested with the above versions.

Installing CATCH with `pip` (or conda), as described below, will install NumPy and SciPy if they are not already installed.

## Setting up a conda environment

_Note: This section is optional, but may be useful to users who are new to Python._

It is generally helpful to install and run Python packages inside of a [virtual environment](https://docs.python.org/3/glossary.html#term-virtual-environment), especially if you have multiple versions of Python installed or use multiple packages.
This can prevent problems when upgrading, conflicts between packages with different requirements, installation issues that arise from having different Python versions available, and more.

One option to manage packages and environments is to use [conda](https://conda.io/en/latest/).
A fast way to obtain conda is to install Miniconda: you can download it [here](https://conda.io/en/latest/miniconda.html) and find installation instructions for it [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html#regular-installation).
For example, on Linux you would run:
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Once you have conda, you can [create](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands) an environment for CATCH with Python 3.8:

```bash
conda create -n catch python=3.8
```

Then, you can activate the `catch` environment:

```bash
conda activate catch
```

After the environment is created and activated, you can install CATCH as described below.
You will need to activate the environment each time you use CATCH.

### Downloading and installing

An easy way to setup CATCH is to clone the repository and install with `pip`:

```bash
git clone https://github.com/broadinstitute/catch.git
cd catch
pip install -e .
```

If you do not have write permissions in the installation directory, you may need to supply `--user` to `pip install`.

### Downloading viral sequence data

_Note: This section is optional and not required for using CATCH._

CATCH's repository contains viral sequence data that can be used as input, but downloading this data is not necessary for using CATCH.
Instead, you can now [automatically download](#example-of-running-designpy) the latest genomes to use as input, so this data is no longer kept up-to-date and will eventually be removed.
Nevertheless, the data may still be helpful for some users and is still needed for some [unit tests](#testing).

We use [Git LFS](https://git-lfs.github.com) to version and store this data.
If you wish to use this data, you'll need to install [Git LFS](https://git-lfs.github.com).
After installing it, you can download the viral sequence data by running:

```bash
git lfs install
git lfs pull
```

from inside the `catch` project directory.

Depending on your setup, providing `-e` to `pip` during [installation](#downloading-and-installing) may be necessary for CATCH to access this data.

### Testing

_Note: This section is optional and not required for using CATCH._

CATCH uses Python's `unittest` framework.
Some of these tests, in order to pass, require you to have [downloaded](#downloading-viral-sequence-data) viral sequence data.
To execute all tests, run:

```bash
python -m unittest discover
```

### Alternative installation approach: conda

CATCH is also available through the [conda](https://conda.io/) package manager as part of the [bioconda](https://bioconda.github.io/) channel.
If you use conda, the easiest way to install CATCH is by running:

```bash
conda install -c bioconda catch
```

## Using CATCH

There are 3 ways of using CATCH to consider.
They are all related &mdash; in that they wrap around the command described in [Option #1](#option-1-cautious-defaults-and-one-setting-of-hybridization-criteria) &mdash; but they can have large differences in their required computational resources (runtime and memory) and the design criteria they use.
If you want a command that is quick to run without needing close familiarity with CATCH, we recommend [Option #2](#option-2-pragmatic-defaults-for-large-diverse-input-and-one-setting-of-hybridization-criteria).
If you have the time to design a complex probe set that works well in practice and/or carefully read documentation, we recommend [Option #1](#option-1-cautious-defaults-and-one-setting-of-hybridization-criteria) and [Option #3](#option-3-variable-hybridization-criteria-across-many-taxa).

### Option #1: Cautious defaults and one setting of hybridization criteria

This way of running CATCH uses [`design.py`](./bin/design.py).
To see details on all of the arguments that the program accepts, run:

```bash
design.py --help
```

[`design.py`](./bin/design.py) requires one or more `dataset`s that specify input sequence data to target, as well as a path to which the probe sequences are written:

```bash
design.py [dataset] [dataset ...] -o OUTPUT
```

Each `dataset` can be one of several input formats:
* A path to a FASTA file.
* An NCBI taxonomy ID, for which sequences will be automatically downloaded.
This is specified as `download:TAXID` where TAXID is the taxonomy ID.
CATCH will fetch all accessions (representing whole genomes) for this taxonomy and download the sequences.
For viruses, NCBI taxonomy IDs can be found via the [Taxonomy Browser](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Undef&id=10239).
* If you [downloaded](#downloading-viral-sequence-data) viral sequence data, a label for one of [550+ viral datasets](./catch/datasets/README.md) (e.g., `human_immunodeficiency_virus_1` or `zika`) distributed as part of this package.
Each of these datasets includes all available whole genomes (genome neighbors) in [NCBI's viral genome data](https://www.ncbi.nlm.nih.gov/genome/viruses/) for a species that has human as a host, as of Oct. 2018.

The probe sequences are written to OUTPUT in FASTA format.

**This program sets cautious hybridization criteria by default (e.g., tolerate no probe-target mismatches) and does not automatically enable options that, in practice, reduce runtime and memory usage.**
**The program was built and described in our [publication](https://www.nature.com/articles/s41587-018-0006-x) as a subroutine of CATCH applied to one viral species at a time.**
**Nevertheless, it can be applied to multiple taxa and to non-viral taxa with appropriate settings: if `dataset`s encompass multiple taxa and/or include large genomes, please (i) consider the arguments described below that reduce runtime and memory usage; or (ii) consider running CATCH as described in [Option #2](#option-2-pragmatic-defaults-for-large-diverse-input-and-one-setting-of-hybridization-criteria) or [Option #3](#option-3-variable-hybridization-criteria-across-many-taxa).**

Below are several important arguments to `design.py`:

* `-pl/--probe-length PROBE_LENGTH`/`-ps/--probe-stride PROBE_STRIDE`: Design probes to be PROBE_LENGTH nt long, and generate candidate probes using a stride of PROBE_STRIDE nt.
(Default: 100 and 50.)
* `-m/--mismatches MISMATCHES`: Tolerate up to MISMATCHES mismatches when determining whether a probe covers a target sequence.
Higher values lead to fewer probes.
This value can considerably affect runtime, with lower runtime at smaller values.
Also, see `-l/--lcf-thres` and `--island-of-exact-match` for adjusting additional hybridization criteria, as described by the output of `design.py --help`.
(Default: 0.)
* `-c/--coverage COVERAGE`: Guarantee that at least COVERAGE of each target genome is captured by probes, where COVERAGE is either a fraction of a genome or a number of nucleotides.
Higher values lead to more probes.
(Default: 1.0 &mdash; i.e., whole genome.)
* `-e/--cover-extension COVER_EXTENSION`: Assume that a probe will capture both the region of the sequence to which it hybridizes, as well as COVER_EXTENSION nt on each side of that.
This parameter is reasonable because library fragments are generally longer than the capture probes, and its value may depend on the library fragment length.
Higher values lead to fewer probes.
(Default: 0.)
* `-i/--identify`: Design probes to perform differential identification.
This is typically used with small values of COVERAGE and &gt;1 specified `dataset`s.
Probes are designed such that each `dataset` should be captured by probes that are unlikely to hybridize to other `dataset`s.
Also, see `-mt/--mismatches-tolerant`, `-lt/--lcf-thres-tolerant`, and `--island-of-exact-match-tolerant`, as described by the output of `design.py --help`.
* `--avoid-genomes dataset [dataset ...]`: Design probes to be unlikely to hybridize to any of these datasets.
Also, see `-mt/--mismatches-tolerant`, `-lt/--lcf-thres-tolerant`, and `--island-of-exact-match-tolerant`, as described by the output of `design.py --help`.
* `--add-adapters`: Add PCR adapters to the ends of each probe sequence.
This selects adapters to add to probe sequences so as to minimize overlap among probes that share an adapter, allowing probes with the same adapter to be amplified together.
(See `--adapter-a` and `--adapter-b` too.)
* `--custom-hybridization-fn PATH FN`: Specify a function, for CATCH to dynamically load, that implements a custom model of hybridization between a probe and target sequence.
See `design.py --help` for details on the expected input and output of this function.
If not set, CATCH uses its default model of hybridization based on `-m/--mismatches`, `-l/--lcf-thres`, and `--island-of-exact-match`.
Relatedly, see `--custom-hybridization-fn-tolerant` as described by the output of `design.py --help`.

Several arguments change the design process in a way that can considerably reduce the computational burden needed for design &mdash; especially for large and highly diverse inputs.
The trade-off of these arguments is an increase in the number of output probes, but this increase is usually small (<10%).
The arguments are:
* `--filter-with-lsh-minhash FILTER_WITH_LSH_MINHASH`: Use locality-sensitive hashing to reduce the space of candidate probes by filtering near-duplicate probes.
This can significantly improve runtime and memory requirements when the input is especially large and diverse.
See the output of `design.py --help` for detail on how this argument works and what its value means.
FILTER_WITH_LSH_MINHASH should generally be around 0.5 to 0.7, and 0.6 is a reasonable starting place.
One caveat: when specifying low values of MISMATCHES (~0, 1, or 2), using this argument can cause the coverage provided for each genome to be slightly less than desired with `--coverage` (e.g., 98% if COVERAGE is 100%).
* `--filter-with-lsh-hamming FILTER_WITH_LSH_HAMMING`: Similar to `--filter-with-lsh-minhash`, except uses a different technique for locality-sensitive hashing.
The improvement in runtime and memory usage is usually less than when using `--filter-with-lsh-minhash`, but the value to specify is more intuitive &mdash; namely, it filters near-duplicate probes, using FILTER_WITH_LSH_HAMMING as the maximum Hamming distance at which to call near-duplicates.
FILTER_WITH_LSH_HAMMING should be commensurate with (but not greater than) MISMATCHES; for example, if MISMATCHES is 5, we recommend setting this value to 3 or 4.
The same caveat regarding coverage noted for `--filter-with-lsh-minhash` also applies to this argument.
* `--cluster-and-design-separately CLUSTER_AND_DESIGN_SEPARATELY`: Cluster input sequences prior to design by computing their MinHash signatures and comparing them.
Then, design probes separately on each cluster and merge the resulting probes.
Like the above arguments, this is another option to improve runtime and memory requirements on large and diverse input.
CLUSTER_AND_DESIGN_SEPARATELY gives the inter-cluster distance threshold for merging clusters, expressed as average nucleotide dissimilarity (1-ANI, where ANI is average nucleotide identity).
Values must be in (0.1, 0.5] and generally should be around 0.1 to 0.2; we recommend 0.15 as a starting place (i.e., cluster at 15% nucleotide dissimilarity).
Also, see `--cluster-and-design-separately-method`, as described by the output of `design.py --help`, for more control over this argument and details on how clustering works.
* `--cluster-from-fragments CLUSTER_FROM_FRAGMENTS`: Break input sequences into fragments of length CLUSTER_FROM_FRAGMENTS nt, and cluster these fragments rather than whole input sequences.
This option requires setting `--cluster-and-design-separately` as well, and proceeds as described for that argument.
This option can be useful for large genomes (e.g., bacterial genomes) in which probes for different chunks can be designed independently.
The fragment length should be around 50,000 nt and 50000 is a reasonable starting place for CLUSTER_FROM_FRAGMENTS.

Note that these arguments may slightly increase runtime for small inputs so, when also considering the trade-offs noted above, we do not recommend them for every application.
Nevertheless, we recommend trying these arguments for large and diverse input (e.g., large genomes, many genomes, or highly divergent genomes), or in any case when runtime or memory usage presents an obstacle to running CATCH.

[`design.py`](./bin/design.py) accepts many other arguments that can be useful for particular design applications, to improve resource requirements, or for debugging; `design.py --help` describes all of them.

### Option #2: Pragmatic defaults for large, diverse input and one setting of hybridization criteria

This way of running CATCH uses [`design_large.py`](./bin/design_large.py).
This program wraps [`design.py`](./bin/design.py), setting default hybridization criteria that are reasonable in practice and enabling, by default, options that lower runtime and memory usage on large, highly diverse input.
The format for specifying inputs and outputs follows the description [above](#option-1-cautious-defaults-and-one-setting-of-hybridization-criteria).

In particular, [`design_large.py`](./bin/design_large.py):
* Sets default hybridization criteria that have worked well in practice: `--mismatches 5 --cover-extension 50`
* Turns on options to lower runtime and memory usage: `--filter-with-lsh-minhash 0.6 --cluster-and-design-separately 0.15 --cluster-from-fragments 50000`
* Uses all available CPUs to parallelize computations

The values of these arguments, as well as all of the arguments described [above](#option-1-cautious-defaults-and-one-setting-of-hybridization-criteria), can still be customized.

To see details on all of the arguments that the program accepts, including their default values, as well as more information about input and output, run:

```bash
design_large.py --help
```

### Option #3: Variable hybridization criteria across many taxa

This way of using CATCH combines multiple runs of [`design.py`](./bin/design.py) with one run of [`pool.py`](./bin/pool.py).
It is most true to how we built CATCH, described it in our [publication](https://www.nature.com/articles/s41587-018-0006-x), and apply it in practice.

[`design.py`](./bin/design.py) and [`design_large.py`](./bin/design_large.py) require a particular selection of hybridization criteria, but it may not be appropriate to use a single criterion uniformly across all taxa in a multi-taxa probe set.
Using [`pool.py`](./bin/pool.py), CATCH can find optimal hybridization settings that are allowed to vary across taxa, under a specified limit on the total number of probes (e.g., synthesis array size).
This process accounts for taxa having different degrees of variation, so that CATCH can balance a probe set's composition between the less diverse species needing few probes with the more diverse species that might otherwise dominate the probe set.
In other words, more diverse species are allowed lower enrichment efficiencies, in order to reduce the number of probes they require, relative to the less diverse species.
Species can also be weighted during this process if there is a greater likelihood of encountering some in a sample than others.

[`pool.py`](./bin/pool.py) searches over a space of potential probe sets to solve a constrained optimization problem.
To see details on all the arguments that the program accepts, run:

```bash
pool.py --help
```

In practice, the process follows a MapReduce-like model.
First, run [`design.py`](./bin/design.py) on each taxon (or other choice of dataset) over a grid of parameter values that spans a reasonable domain.
These parameters could encompass the number of mismatches (`--mismatches`) or cover extension (`--cover-extension`), or other parameters in a custom hybridization model.
Then, create a table that provides a probe count for each taxon and choice of parameters; the table must be a TSV, in a format like [this](./catch/pool/tests/input/num-probes.V-WAfr.201506.tsv).
Now, use this table as input to [`pool.py`](./bin/pool.py):

```bash
pool.py INPUT_TSV TARGET_PROBE_COUNT OUTPUT_TSV
```
where INPUT_TSV is a path to the table described above, TARGET_PROBE_COUNT is a constraint on the number of probes to allow in the pool, and OUTPUT_TSV is a path to a file to which the program will write the optimal parameter values.

Below are two arguments that generalize the search:
* `--loss-coeffs COEFF [COEFF ...]`: Specify coefficients on parameters in the objective function.
This allows you to adjust how conservative each parameter is treated relative to others.
(Default: 1 for mismatches and 1/100 for cover extension.)
* `--dataset-weights WEIGHTS_TSV`: Assign a weight for each dataset to use in the objective function, where WEIGHTS_TSV is a path to a table that provides a weight for each dataset.
This allows you to seek that probes in the pooled design be more sensitive for some taxa than others, for instance, if you believe you are more likely to encounter some taxa in a sample.
(Default: 1 for all datasets.)

Each run of [`pool.py`](./bin/pool.py) may yield a slightly different output based on the (random) initial guess.
You could run it multiple times and select the output that has the smallest loss, which is written to standard output at the end of the program.

## Examples

### Example of running [`design.py`](./bin/design.py)

Below is an example of designing probes to target a single taxon.

```bash
design.py download:64320 -pl 75 -m 2 -l 60 -e 50 -o zika-probes.fasta --verbose
```

This will download whole genomes of Zika virus (NCBI taxonomy ID [64320](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=64320)) and design probes that:
* are 75 nt long (`-pl 75`)
* capture the entirety of each genome under a model that a probe hybridizes to a region if the longest common substring, up to 2 mismatches (`-m 2`), between a probe and target is at least 60 nt (`-l 60`)
* assume 50 nt on each side of the hybridization is captured as well (`-e 50`)

and will save them to `zika-probes.fasta`.

It will provide detailed output during runtime (`--verbose`) and yield about 600 probes.
Note that using `-l 75` here (or, equivalently, leaving out `-l`) will run significantly faster, but result in more probes.
Also, note that the input can be `zika` to use the `zika` dataset distributed with CATCH, or a path to any custom FASTA file.

### Example of running [`design_large.py`](./bin/design_large.py)

Below is an example of designing probes with larger and more diverse input than that provided in the above example:

```bash
design_large.py download:64320 download:12637 -o zika-and-dengue-probes.fasta --verbose
```

This will download whole genomes of Zika virus (NCBI taxonomy ID [64320](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=64320)) and dengue virus (NCBI taxonomy ID [12637](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=12637)).
It will then design probes to enrich both species and save them to `zika-and-dengue-probes.fasta`.

The command will provide detailed output during runtime (`--verbose`) and yield about 3,200 probes.
It will take about 1 hour to run (with 8 CPUs).

For hybridization criteria, [`design_large.py`](./bin/design_large.py) will use default values (available in the output of `design_large.py --help`), which can be overridden with custom values if desired.

### Example of running [`pool.py`](./bin/pool.py)

[Here](./catch/pool/tests/input/num-probes.V-WAfr.201506.tsv) is a table listing probe counts used in the design of the V-WAfr probe set (provided [here](./probe-designs) and described in our [publication](https://www.nature.com/articles/s41587-018-0006-x)).
It provides counts for each dataset (species) and combination of two parameters (mismatches and cover extension) that were allowed to vary in CATCH's design; the table was produced by running [`design.py`](./bin/design.py) for each species over a grid of parameter values.
Below is an example of designing the V-WAfr probe set using this table as input.

```bash
pool.py num-probes.V-WAfr.201506.tsv 90000 params.V-Wafr.201506.tsv --round-params 1 10
```

This will search for parameters that yield at most 90,000 probes across the species, and will output those to `params.V-Wafr.201506.tsv`.
Because the search is over a continuous space, here we use `--round-params 1 10` to set each value of the mismatches parameter to an integer and each value of the cover extension parameter to a multiple of 10 while still meeting the constraint on probe count.
The pooled design yields about 89,950 probes, depending on the initial guess.

## Contributing

We welcome contributions to CATCH.
This can be in the form of an [issue](https://github.com/broadinstitute/catch/issues) or [pull request](https://github.com/broadinstitute/catch/pulls).
If you have questions, please create an [issue](https://github.com/broadinstitute/catch/issues) or email **Hayden Metsky** &lt;hayden@broadinstitute.org&gt;.

## Citation

For details on how CATCH works, please refer to our publication in [Nature Biotechnology](https://www.nature.com/articles/s41587-018-0006-x). If you find CATCH useful to your work, please cite our paper as:
  * Metsky HC and Siddle KJ _et al_. Capturing sequence diversity in metagenomes with comprehensive and scalable probe design. _Nature Biotechnology_, **37**(2), 160&ndash;168 (2019). doi: 10.1038/s41587-018-0006-x

## License

CATCH is licensed under the terms of the [MIT license](./LICENSE).

