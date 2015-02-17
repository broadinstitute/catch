"""A dataset with 99 sequences of Ebola taken from the 2014
outbreak."""

# Author: Hayden Metsky <hayden@mit.edu>

from os.path import join
from os.path import dirname

FASTA_RELATIVE_PATH = "data/ebola2014.fasta"


def fasta_path():
  return join(dirname(__file__), FASTA_RELATIVE_PATH)
