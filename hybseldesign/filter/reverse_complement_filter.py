"""Filters input probes by returning each input probe as well as its
reverse complement.

The number of output probes is twice the number of input probes.
In the output, the reverse complements are interspersed with
the original input probes --- i.e., if the input probes are
[P1, P2, P3], the output is [P1, P1r, P2, P2r, P3, P3r] where
Pnr is the reverse complement of Pn.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

from hybseldesign.filter.base_filter import BaseFilter


class ReverseComplementFilter(BaseFilter):

  def _filter(self, input):
    output = []
    for p in input:
      output += [ p ]
      output += [ p.reverse_complement() ]
    return output
