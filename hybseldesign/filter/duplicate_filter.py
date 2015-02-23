"""Filters probes by returning a list of non-duplicate probes from
the given input probes.

This only removes probes that are identical and makes no
guarantees about the order of what is returned (i.e., it is
not stable with respect to the input order).
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

from hybseldesign.filter.base_filter import BaseFilter


class DuplicateFilter(BaseFilter):

  def _filter(self, input):
    return list(set(input))
