"""A dataset with 78 Marburg samples.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

from os.path import join
from os.path import dirname

FASTA_RELATIVE_PATH = "data/marburg.fasta"


def fasta_path():
  return join(dirname(__file__), FASTA_RELATIVE_PATH)
