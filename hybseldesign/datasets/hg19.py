"""The human genome 19 dataset.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'


def set_fasta_path(absolute_path):
  FASTA_ABSOLUTE_PATH = absolute_path


def fasta_path():
  try:
    FASTA_ABSOLUTE_PATH
  except NameError:
    raise NameError("FASTA path not set")
  else:
    return FASTA_ABSOLUTE_PATH
