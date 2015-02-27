"""A dataset with 140 Ebola samples taken from Zaire.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

from os.path import join
from os.path import dirname

FASTA_RELATIVE_PATH = "data/ebola_zaire.fasta"


def fasta_path():
  return join(dirname(__file__), FASTA_RELATIVE_PATH)
