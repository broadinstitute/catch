[![Build Status](https://magnum.travis-ci.com/broadinstitute/hybsel_design.svg?token=1e54w9HsqGg5ZQzQ9ruW&branch=master)](https://magnum.travis-ci.com/broadinstitute/hybsel_design)
[![Coverage Status](https://coveralls.io/repos/broadinstitute/hybsel_design/badge.svg?branch=master&t=b06RAq)](https://coveralls.io/r/broadinstitute/hybsel_design?branch=master)

hybsel_design
=============

hybsel_design is a Python package for designing probes to be used in hybrid selection.

## Dependencies

hybsel_design requires:
* Python >= 3.5
* NumPy >= 1.9.0

hybsel_design comes with viral sequence data to use as input.
If you wish to use this data, you'll need to install [Git LFS](https://git-lfs.github.com) so the data can be downloaded.
To install Git LFS on a Broad machine, simply run:

```bash
use .git-lfs-2.0.2
git lfs install
```

## Install

The easiest way to install hybsel_design is with `pip`:

```bash
git clone git@github.com:broadinstitute/hybsel_design.git
cd hybsel_design
pip install -e .
```

## Testing

hybsel_design uses Python's `unittest` framework. To execute all unit tests, run:

```bash
python -m unittest discover
```

## Running

The main script to design probes is `bin/make_probes.py`.
To see the arguments that the program accepts, run:

```bash
bin/make_probes.py -h
```

## Examples

Here is a command that uses some common options:

```bash
bin/make_probes.py -d lassa -pl 75 -m 2 -e 50 -o lassa-probes.fasta
```

This will design probes that:
* target Lassa virus (`-d lassa`)
* are 75 nt long (`-pl 75`)
* capture the entirety of each input genome, tolerating up to 2 mismatches (`-m 2`)
* assume 50 nt on each side of the hybridization is captured as well (`-e 50`)

and will save them to `lassa-probes.fasta`.

To use a custom input instead of a provided dataset, use `-d custom:[path-to-sequences.fasta]`.
