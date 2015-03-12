"""A dataset with 53 Ebola samples taken from places that do not
include Zaire. These samples do not include those in the ebola2014
dataset.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

from os.path import join
from os.path import dirname

FASTA_RELATIVE_PATH = "data/ebola_nonzaire.fasta"


def fasta_path():
  return join(dirname(__file__), FASTA_RELATIVE_PATH)
