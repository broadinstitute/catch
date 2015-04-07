"""Holds convenience classes and functions for working with datasets.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

import logging
from os.path import join
from os.path import dirname

logger = logging.getLogger(__name__)


"""Holds one or more individual (sampled) genomes, usually from
the same species.
"""
class GenomesDataset:

  """
  '__name__' is the name of this dataset (module).
  '__file__' is the path to the module corresponding to this dataset.
  """
  def __init__(self, __name__, __file__):
    self.__name__ = __name__
    self.__file__ = __file__

  def is_multi_chr(self):
    try:
      self.chrs
    except AttributeError:
      return False
    else:
      return True


"""This is used when the genome held by this dataset can be fully
specified by one sequence (i.e., there are not multiple chromosomes).
All of the sample/individual genomes should be in a single FASTA file
set by set_fasta_path().
"""
class GenomesDatasetSingleChrom(GenomesDataset):

  """
  '__name__' is the name of this dataset (module).
  '__file__' is the path to the module corresponding to this dataset.
  """
  def __init__(self, __name__, __file__):
    GenomesDataset.__init__(self, __name__, __file__)
    self.fasta_path = None

  """Set the FASTA file for the genome samples stored by this
  dataset.

  When 'relative' is True, path is assumed to be relative to the
  path of the dataset file (self.__file__).
  """
  def set_fasta_path(self, path, relative=False):
    if relative:
      path = join(dirname(self.__file__), path)
    self.fasta_path = path


"""This is used when the genome held by this dataset has more than one
chromosome.  Each FASTA file offers just one sample/individual genome
in which each sequence in the file corresponds to a chromosome; there
is typically more than one FASTA file for the dataset, and they should
each be added with add_fasta_path().
"""
class GenomesDatasetMultiChrom(GenomesDataset):

  """
  '__name__' is the name of this dataset (module).
  '__file__' is the path to the module corresponding to this dataset.
  'chrs' is a list of chromosomes for the genome held by this
  dataset; these should be the header of each sequence in the FASTA
  files.
  """
  def __init__(self, __name__, __file__, chrs):
    if len(chrs) == 1:
      logger.critical(("Typically when there is just one chromosome, "
                       "the GenomesDatasetSingleChrom class should "
                       "used. Even though there is just one "
                       "chromosome, GenomesDatasetMultiChrom will "
                       "assume that each sequence (chromosome) is in "
                       "a separate FASTA file."))

    GenomesDataset.__init__(self, __name__, __file__)
    self.chrs = chrs
    self.fasta_paths = []

  """Add a FASTA file giving one sample/individual genome for this
  dataset.

  When 'relative' is True, path is assumed to be relative to the
  path of the dataset file (self.__file__).
  """
  def add_fasta_path(self, path, relative=False):
    if relative:
      path = join(dirname(self.__file__), path)
    self.fasta_paths += [path]

