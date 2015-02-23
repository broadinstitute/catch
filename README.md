[![Build Status](https://magnum.travis-ci.com/broadinstitute/hybsel_design.svg?token=1e54w9HsqGg5ZQzQ9ruW&branch=master)](https://magnum.travis-ci.com/broadinstitute/hybsel_design)
[![Coverage Status](https://coveralls.io/repos/broadinstitute/hybsel_design/badge.svg?branch=master)](https://coveralls.io/r/broadinstitute/hybsel_design)

hybsel_design
=============

hybsel_design is a Python package for designing probes to be used in hybrid selection.

## Dependencies

hybsel_design is tested under Python 2.7, but should also work under Python 2.6 and newer versions. It requires NumPy; it is tested with NumPy 1.9.0, but should work with other versions as well.


## Install and Develop

### Install (not intended for development)
hybsel_design is installed with `setuptools`. To install it in your local user directory, run:

```
$ python setup.py install --user
```

### Install for development
To setup the package in your local user directory in development mode, in which changes are reflected immediately, run:

```
$ python setup.py develop --install-dir=~/.local/lib/python2.7/site-packages/
```

on a Broad machine. Or run:

```
$ python setup.py develop --install-dir=~/Library/Python/2.7/lib/python/site-packages/
```

on OS X.

## Testing

hybsel_design uses Python's `unittest` framework. To execute all unit tests, run:

```
$ python -m unittest discover
```

It is also possible to run unit tests on a particular package or module, or even run a particular test. For example:

```
$ python -m unittest hybseldesign.utils.tests.test_seq_io
```

runs tests for the `hybseldesign.utils.seq_io` module.
