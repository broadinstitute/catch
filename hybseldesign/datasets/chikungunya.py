"""A dataset with 213 Chikungunya sequences. There are 181 complete
genomes from NCBI and 32 isolates from VBRC.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

from os.path import join
from os.path import dirname

FASTA_RELATIVE_PATH = "data/chikungunya.fasta"


def fasta_path():
  return join(dirname(__file__), FASTA_RELATIVE_PATH)
