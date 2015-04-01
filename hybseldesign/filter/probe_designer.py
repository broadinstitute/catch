"""Designs probes by generating a set of candidate probes and
passing these through a set of one or more given filters.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

import logging

from hybseldesign.filter import candidate_probes

logger = logging.getLogger(__name__)


class ProbeDesigner:

  """Creates a ProbeDesigner from a collection 'genomes' of grouped
  genomes and an (ordered) list 'filters' of filters.

  genomes is a list [g_1, g_2, g_m] of m groupings of genomes, where
  each g_i is a list of genome.Genomes belonging to group i. For
  example, a group may be a species and each g_i would be a list of
  the target genomes of species i. Each filter should be an instance
  of a subclass of Filter.

  When replicate_first_version is True, candidate probes are
  explicitly designed in a (buggy) way intended to replicate the
  first version (Matlab code) of probe design.
  """
  def __init__(self, genomes, filters, replicate_first_version=False):
    self.genomes = genomes
    self.filters = filters
    self.replicate_first_version = replicate_first_version

  """Generates the set of candidate probes and runs these through
  the provided filters.

  Stores the candidate probes in self.candidate_probes and the
  probes processed by the filters in self.final_probes.
  """
  def design(self):
    if self.replicate_first_version:
      replicate_args = {
                'insert_bugs': True,
                'move_all_n_string_flanking_probes_to_end': True
      }
    else:
      replicate_args = {}

    logger.info("Building candidate probes from target sequences")
    self.candidate_probes = []
    for genomes_from_group in self.genomes:
      for g in genomes_from_group:
        self.candidate_probes += candidate_probes.\
            make_candidate_probes_from_sequences(g.seqs,
              **replicate_args)

    probes = self.candidate_probes
    for f in self.filters:
      logger.info("Starting filter %s", f.__class__.__name__)
      f.target_genomes = self.genomes
      probes = f.filter(probes)
    self.final_probes = probes

