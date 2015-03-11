"""A dataset with the full genomes of all deer tick borne pathogens
(41 in total). These are mostly bacterial genomes with one genome
per bacteria (i.e., a sample of n=1 for each pathogen).
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

from os.path import join
from os.path import dirname

FASTA_RELATIVE_PATH = "data/all_deer_tick_borne.fasta"


def fasta_path():
  return join(dirname(__file__), FASTA_RELATIVE_PATH)
