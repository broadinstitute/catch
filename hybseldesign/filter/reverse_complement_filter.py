"""Filters input probes by returning each input probe as well as its
reverse complement.

The number of output probes is twice the number of input probes.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

from hybseldesign.filter.base_filter import BaseFilter


class ReverseComplementFilter(BaseFilter):

  def _filter(self, input):
    reverse_complements = [p.reverse_complement() for p in input]
    return input + reverse_complements
