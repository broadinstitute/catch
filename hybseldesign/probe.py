"""Structure(s) storing probes and methods for directly working
with them.
"""

# Author: Hayden Metsky <hayden@mit.edu>

import numpy as np


"""Immutable sequence representing a probe/bait.
"""
class Probe:
  def __init__(self, seq):
    self.seq = seq
    self.hash = hash(seq.tostring())

  """Counts the number of mismatches between self and other.

  They must be the same length.
  """
  def mismatches(self, other):
    return self.mismatches_at_offset(other, 0)
  
  """Counts the number of mismatches between self and other after
  other is shifted by offset bp.

  Specifically, returns the minimum of the number of mismatches
  after shifting other to the left and to the right. The two
  sequences must be the same length. offset must be less than the
  length of the sequence.
  """
  def mismatches_at_offset(self, other, offset):
    if len(self.seq) != len(other.seq):
      raise ValueError("Sequences must be of same length")
    if offset < 0 or offset >= len(other.seq):
      raise ValueError("Invalid offset value " + str(offset))
    if offset == 0:
      left = np.sum(self.seq != other.seq)
      right = left
    else:
      left = np.sum(self.seq[:-offset] != other.seq[offset:])
      right = np.sum(self.seq[offset:] != other.seq[:-offset])
    return min(left, right)

  """Returns a probe that is the reverse-complement of this
  probe.

  Using rc_map.get(b, b) ensures that this can process bases
  like 'N'. It returns the base itself (e.g., 'N') if it is
  not either 'A', 'T', 'C', or 'G'.
  """
  def reverse_complement(self):
    rc_map = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    rc_seq = np.array([rc_map.get(b, b) for b in self.seq[::-1]],
                dtype='S1')
    return Probe(rc_seq)

  def __hash__(self):
    return self.hash

  def __eq__(self, other):
    return isinstance(other, Probe) and \
      np.array_equal(self.seq, other.seq)

  def __str__(self):
    return self.seq.tostring()

  def __repr__(self):
    return str(self.seq)

  """Make a Probe from a Python str.
  """
  @staticmethod
  def from_str(seqStr):
    return Probe(np.fromstring(seqStr, dtype='S1'))
