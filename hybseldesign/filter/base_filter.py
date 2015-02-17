"""An abstract class representing a filter through which candidate
probes are processed.

A filter may operate in the more traditional sense of reducing the
input probes to a smaller number (i.e., only allowing certain probes
to pass through); an example would be a filter the eliminates
duplicates. A filter may also operate in the more computing sense
of simply processing the input probes and returning the list of
processed probes; the output/processed probes may be altered
versions of the input or there may even be more output probes than
input probes.

All subclasses must implement a _filter(..) method that returns a
list of probes after processing from the given input list. This
saves the input probes in self.input_probes and the output probes in
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

