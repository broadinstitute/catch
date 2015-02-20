"""Designs probes by generating a set of candidate probes and
passing these through a set of one or more given filters.
"""

# Author: Hayden Metsky <hayden@mit.edu>

from hybseldesign.filter import candidate_probes


class ProbeDesigner:

  """Creates a ProbeDesigner from a list 'seqs' of sequences and an
  (ordered) list 'filters' of filters.

  Each filter should be an instance of a subclass of Filter.

  When replicate_first_version is True, candidate probes are
  explicitly designed in a (buggy) way intended to replicate the
  first version (Matlab code) of probe design.
  """
  def __init__(self, seqs, filters, replicate_first_version=False):
    self.seqs = seqs
    self.filters = filters
    self.replicate_first_version = replicate_first_version

  """Generates the set of candidate probes and runs these through
  the provided filters.

  Stores the candidate probes in self.candidate_probes and the
  probes processed by the filters in self.final_probes.
  """
  def design(self):
    if self.replicate_first_version:
      replicate_args = { 'insert_bugs': True }
    else:
      replicate_args = {}
    self.candidate_probes = candidate_probes.\
        make_candidate_probes_from_sequences(self.seqs,
          **replicate_args)

    probes = self.candidate_probes
    for f in self.filters:
      probes = f.filter(probes)
    self.final_probes = probes

