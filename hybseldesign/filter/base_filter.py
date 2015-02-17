"""An abstract class representing a filter through which candidate
probes are reduced.

All subclasses must implement a _filter(..) method that returns a
list of probes filtered from the given input list. This saves the
input probes in self.input_probes and the output probes in
self.output_probes.
"""

# Author: Hayden Metsky <hayden@mit.edu>


class BaseFilter:

  def filter(self, input):
    self.input_probes = input
    filtered = self._filter(input)
    self.output_probes = filtered
    return filtered

  def _filter(self, input):
    raise Exception(("A subclass of BaseFilter must implement "
                     "_filter(..)"))

