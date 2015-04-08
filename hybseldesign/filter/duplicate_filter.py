"""Removes duplicates from an input list of probes.

This acts as a filter on the probes by returning a list of
non-duplicate probes from the given input probes. This only removes
probes that are identical to another. It preserves the order of
the input.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

from collections import OrderedDict

from hybseldesign.filter.base_filter import BaseFilter


class DuplicateFilter(BaseFilter):

  """Filter that removes duplicates.
  """

  def _filter(self, input):
    """Return a subset of the input probes.
    """
    # `return list(set(input))` would be a short way to produce
    # non-duplicate probes, but would not preserve the input
    # order. Instead, preserve the order with an OrderedDict.
    return list(OrderedDict.fromkeys(input))

