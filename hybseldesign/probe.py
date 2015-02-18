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

  offset can be negative (corresponding to other being shifted left)
  or positive (corresponding to other being shifted right).
  """
  def mismatches_at_offset(self, other, offset):
    if len(self.seq) != len(other.seq):
      raise ValueError("Sequences must be of same length")
    if abs(offset) >= len(other.seq):
      raise ValueError("Invalid offset value " + str(offset))
    if offset == 0:
      return np.sum(self.seq != other.seq)
    elif offset < 0:
      return np.sum(self.seq[:offset] != other.seq[-offset:])
    else:
      return np.sum(self.seq[offset:] != other.seq[:-offset])

  """Counts the minimum number of mismatches between self and other
  as other is shifted with an offset between -max_shift and
  +max_shift relative to self.
  """
  def min_mismatches_within_shift(self, other, max_shift):
    return min(self.mismatches_at_offset(other, offset) \
                for offset in xrange(-max_shift, max_shift+1))

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
