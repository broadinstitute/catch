"""Utilities for working with sequence i/o.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

import logging
import numpy as np
import re
from collections import OrderedDict

from hybseldesign import genome

logger = logging.getLogger(__name__)


"""Reads the genomes of the given dataset, and returns these
as a list of genome.Genome.
"""
def read_dataset_genomes(dataset):
  try:
    dataset.CHRS
    read_by_chrs = True
  except AttributeError:
    read_by_chrs = False

  genomes = []

  if read_by_chrs:
    # The genomes in this dataset have more than one chromosome.
    # The dataset should have one or more paths to FASTA files
    # (dataset.fasta_paths()), where each file gives one genome
    # (sample) and the sequences in that file correspond to the
    # different chromosomes. The header of each sequence should be
    # a chromosome in dataset.CHRS.
    logger.debug("Reading dataset %s broken up by chromosome",
        dataset.__name__)
    for fn in dataset.fasta_paths():
      seq_map = read_fasta(fn)
      seqs = OrderedDict()
      for chr in dataset.CHRS:
        if chr not in seq_map:
          raise ValueError(("Chromosome %s is not in the FASTA file "
                            "%s for dataset %s, but should be") %
                            (chr, fn, dataset.__name__))
        seqs[chr] = seq_map[chr]
      genomes += [ genome.Genome.from_chrs(seqs) ]
  else:
    # There is just one sequence (chromosome) for each genome in
    # this dataset. The dataset should have a path to just one FASTA
    # file (dataset.fasta_path()), and the sequences in this file
    # should correspond to the different genomes. The headers of each
    # sequence are ignored.
    logger.debug("Reading dataset %s with one chromosome per genome",
        dataset.__name__)
    seqs = read_fasta(dataset.fasta_path()).values()
    for seq in seqs:
      genomes += [ genome.Genome.from_one_seq(seq) ]

  return genomes


"""Reads the FASTA file fn and return a mapping from the name of
each sequence to the sequence itself, where the sequence is
stored as a native Python string ('str') or numpy array ('np')
as determined by data_type.

The mapping returned is ordered by the order in which the sequence
is encountered in the FASTA file. This helps in particular with
replicating past results, where the input order could affect the
output.

The degenerate bases ('Y','R','W','S','M','K') are replaced with
'N' iff replace_degenerate is True.
"""
def read_fasta(fn, data_type='str', replace_degenerate=True):
  logger.info("Reading fasta file %s", fn)

  degenerate_pattern = re.compile('[YRWSMK]')

  m = OrderedDict()
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
        if replace_degenerate:
          line = degenerate_pattern.sub('N', line)
        m[curr_seq_name] += line

  if data_type == 'str':
    # Already stored sequence as string
    m_converted = m
  elif data_type == 'np':
    m_converted = OrderedDict()
    for seq_name, seq in m.iteritems():
      m_converted[seq_name] = np.fromstring(seq, dtype='S1')
  else:
    raise ValueError("Unknown data_type " + data_type)

  return m_converted


"""A generator that scans through the FASTA file fn and, upon
completing the read of a sequence, yields that sequence, where
the sequence is stored as a native Python string ('str') or numpy
array ('np') as determined by data_type.

The degenerate bases ('Y','R','W','S','M','K') are replaced with
'N' iff replace_degenerate is True.
"""
def iterate_fasta(fn, data_type='str', replace_degenerate=True):
  degenerate_pattern = re.compile('[YRWSMK]')

  def format_seq(seq):
    if data_type == 'str':
      # Already stored as str
      return seq
    elif data_type == 'np':
      return np.fromstring(seq, dtype='S1')
    else:
      raise ValueError("Unknown data_type " + data_tyoe)

  with open(fn) as f:
    curr_seq = ''
    for line in f:
      line = line.rstrip()
      if line.startswith('>'):
        # Yield the current sequence (if there is one) and reset the
        # sequence being read
        if len(curr_seq) > 0:
          yield format_seq(curr_seq)
        curr_seq = ''
      else:
        # Append the sequence
        if replace_degenerate:
          line = degenerate_pattern.sub('N', line)
        curr_seq += line
    if len(curr_seq) > 0:
      yield format_seq(curr_seq)


"""Write the sequences in 'probes' to the file 'out_fn' in FASTA
format.

This writes one probe sequence per line, with a header immediately
preceding the sequence. If set, the header written is the one in
probe.Probe.header. If not set, the probe.Probe.identifier() is used.
"""
def write_probe_fasta(probes, out_fn):
  with open(out_fn, 'w') as f:
    for p in probes:
      if p.header:
        f.write('>' + p.header + '\n')
      else:
        f.write('>probe_%s\n' % p.identifier())
      f.write(p.seq_str + '\n')

