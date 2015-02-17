"""Utilities for working with sequence i/o.
"""

# Author: Hayden Metsky <hayden@mit.edu>

import numpy as np


"""Reads the FASTA file fn and return a mapping from the name of
each sequence to the sequence itself, where the sequence is
stored as a native Python string ('str') or numpy array ('np')
as determined by data_type.
"""
def read_fasta(fn, data_type='str'):
  m = {}
  with open(fn) as f:
    curr_seq_name = ""
    for line in f:
      line = line.rstrip()
      if curr_seq_name == "":
        # Must encounter a new sequence
        assert line.startswith('>')
      if len(line) == 0:
        # Reset the sequence being read on an empty line
        curr_seq_name = ""
      elif line.startswith('>'):
        curr_seq_name = line[1:]
        m[curr_seq_name] = ''
      else:
        # Append the sequence
        m[curr_seq_name] += line

  if data_type == 'str':
    # Already stored sequence as string
    m_converted = m
  elif data_type == 'np':
    m_converted = { seq_name: np.fromstring(seq, dtype='S1') \
                    for seq_name, seq in m.iteritems() }
  else:
    raise ValueError("Unknown data_type " + data_type)

  return m_converted

