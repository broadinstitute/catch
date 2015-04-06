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
  'dataset_file' is the path to the module corresponding to this
  dataset. 
  """
  def __init__(self, dataset_name, dataset_file):
    self.dataset_name = dataset_name
    self.dataset_file = dataset_file

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
  'dataset_file' is the path to the module corresponding to this
  dataset. 
  """
  def __init__(self, dataset_name, dataset_file):
    GenomesDataset.__init__(self, dataset_name, dataset_file)
    self.fasta_path = None

  """Set the FASTA file for the genome samples stored by this
  dataset.

  When 'relative' is True, path is assumed to be relative to the
  path of the dataset file (self.dataset_file).
  """
  def set_fasta_path(self, path, relative=False):
    if relative:
      path = join(dirname(self.dataset_file), path)
    self.fasta_path = path

  def fasta_path(self):
    if self.fasta_path == None:
      raise ValueError("The FASTA path has not yet been set")
    return self.fasta_path


"""This is used when the genome held by this dataset has more than one
chromosome.  Each FASTA file offers just one sample/individual genome
in which each sequence in the file corresponds to a chromosome; there
is typically more than one FASTA file for the dataset, and they should
each be added with add_fasta_path().
"""
class GenomesDatasetMultiChrom(GenomesDataset):

  """
  'dataset_file' is the path to the module corresponding to this
  dataset. 
  'chrs' is a list of chromosomes for the genome held by this
  dataset; these should be the header of each sequence in the FASTA
  files.
  """
  def __init__(self, dataset_name, dataset_file, chrs):
    if len(chrs) == 1:
      logger.critical(("Typically when there is just one chromosome, "
                       "the GenomesDatasetSingleChrom class should "
                       "used. Even though there is just one "
                       "chromosome, GenomesDatasetMultiChrom will "
                       "assume that each sequence (chromosome) is in "
                       "a separate FASTA file."))

    GenomesDataset.__init__(self, dataset_name, dataset_file)
    self.dataset_file = dataset_file
    self.chrs = chrs
    self.fasta_paths = []

  """Add a FASTA file giving one sample/individual genome for this
  dataset.

  When 'relative' is True, path is assumed to be relative to the
  path of the dataset file (self.dataset_file).
  """
  def add_fasta_path(self, path, relative=False):
    if relative:
      path = join(dirname(self.dataset_file), path)
    self.fasta_paths += [path]

  def fasta_paths(self):
    if len(self.fasta_paths) == 0:
      raise ValueError("No FASTA paths have been added yet")
    return self.fasta_paths

