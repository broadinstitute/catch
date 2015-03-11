"""A dataset with 12 Powassan samples. This is a subset of the
genomes in the all_deer_tick_borne dataset.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

from os.path import join
from os.path import dirname

FASTA_RELATIVE_PATH = "data/powassan.fasta"


def fasta_path():
  return join(dirname(__file__), FASTA_RELATIVE_PATH)
