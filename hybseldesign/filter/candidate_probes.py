"""Functions for generating lists of many (likely redundant) probes,
termed candidate probes, from a sequence or list of sequences.
"""

# Author: Hayden Metsky <hayden@mit.edu>

import sys
import re
import numpy as np

from hybseldesign.utils import seq_io
from hybseldesign import probe
from hybseldesign.datasets import ebola2014


"""Generates and returns a list of candidate probes from the sequence
seq.

The probes are probe_length bp long and separated by probeStride bp.
Possible probes that would contains strings of min_n_string_length or
more N's are discarded and, instead, new probes flanking the string
are added.

It is possible (likely, in fact, when there are strings of
N's) that duplicate probes are returned.
"""
def make_candidate_probes_from_sequence(seq, probe_length=100,
    probe_stride=50, min_n_string_length=2):
  if probe_length > len(seq):
    raise ValueError("Invalid probe_length " + str(probe_length))

  if type(seq) is np.ndarray or type(seq) is np.array:
    seq = seq.tostring()
  n_string_query = re.compile('(N{'+str(min_n_string_length)+',})')

  # Make one or two probes based on the subsequence seq[start:end].
  # Namely, if that subsequence contains no string of N's, then it
  # is itself a probe to be added. But if it does contain a string
  # of N's, then add two flanking probes as needed.
  def add_probe_from_subsequence(start, end):
    subseq = seq[start:end]
    probes = []

    # Search for strings of min_n_string_length or more N's in subseq
    # (Add start onto match.start() and match.end() because the
    #  search is within subseq, but the coordinates should be
    #  relative to seq.)
    found_n_string = False
    n_string_starts, n_string_ends = [], []
    for match in n_string_query.finditer(subseq):
      found_n_string = True
      n_string_starts += [start + match.start()]
      n_string_ends += [start + match.end()]

    if found_n_string:
      # There is indeed at least one string of N's, so add flanking
      # probes. Specifically, add a probe to the left of the
      # leftmost string and a probe to the right of the rightmost
      # string.
      # (If we were to add a right probe to the leftmost string, a
      #  a left probe to the rightmost string, or flanking probes on
      #  middle strings, these would overlap strings of N's.)
      leftmost_start = min(n_string_starts)
      rightmost_end = max(n_string_ends)
      if leftmost_start - probe_length >= 0:
        left_subseq = seq[(leftmost_start-probe_length):leftmost_start]
        if not n_string_query.search(left_subseq):
          # Only add the flanking probe if it doesn't contain a
          # string of N's (i.e., don't recursively chase flanking
          # probes)
          probes += [left_subseq]
      if rightmost_end + probe_length <= len(seq):
        right_subseq = seq[rightmost_end:(rightmost_end+probe_length)]
        if not n_string_query.search(right_subseq):
          # Only add the flanking probe if it doesn't contain a
          # string of N's (i.e., don't recursively chase flanking
          # probes)
          probes += [right_subseq]
    else:
      # There's no string of N's, so this subsequence is a valid probe
      probes += [subseq]

    return probes

  # Populate a list of probes
  probes = []
  for start in np.arange(0, len(seq), probe_stride):
    if start + probe_length > len(seq):
      break
    probes += add_probe_from_subsequence(start, start+probe_length)
  if len(seq) % probe_stride != 0:
    # There are bases on the right that were never covered, so add
    # another probe for this
    probes += add_probe_from_subsequence(len(seq)-probe_length, len(seq))

  # Convert the probes from a Python list of Python strings to a
  # list of probe.Probe
  for i in xrange(len(probes)):
    probes[i] = probe.Probe.from_str(probes[i])

  return probes


"""Generates and returns a list of probes from a list of sequences
seqs.

It is possible (perhaps even likely depending on where
the sequences come from) that duplicate probes are returned.
"""
def make_candidate_probes_from_sequences(seqs, probe_length=100,        
    probe_stride=50, min_n_string_length=2):
  if type(seqs) != list:
    raise ValueError("seqs must be a list of sequences")
  if len(seqs) == 0:
    raise ValueError("seqs must have at least one sequence")
  for seq in seqs:
    if type(seq) != str:
      raise ValueError("seqs must be a list of Python strings")

  probes = []
  for seq in seqs:
    probes += make_candidate_probes_from_sequence(seq,
                probe_length=probe_length, probe_stride=probe_stride,
                min_n_string_length=min_n_string_length)
  return probes

