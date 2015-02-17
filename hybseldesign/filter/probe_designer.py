"""Designs probes by generating a set of candidate probes and
filtering these through a set of one or more given filters.
"""

# Author: Hayden Metsky <hayden@mit.edu>

from hybseldesign.filter import candidate_probes


class ProbeDesigner:

  """Creates a ProbeDesigner from a list 'seqs' of sequences and an
  (ordered) list 'filters' of filters.

  Each filter should be an instance of a subclass of Filter.
  """
  def __init__(self, seqs, filters):
    self.seqs = seqs
    self.filters = filters

  """Generates the set of candidate probes and runs these through
  the provided filters.

  Stores the candidate probes in self.candidate_probes and the
  probes that passed the filters in self.final_probes.
  """
  def design(self):
    self.candidate_probes = candidate_probes.\
        make_candidate_probes_from_sequences(self.seqs)
    probes = self.candidate_probes
    for f in self.filters:
      probes = f.filter(probes)
    self.final_probes = probes

