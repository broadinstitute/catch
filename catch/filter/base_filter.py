"""Wrappers around a filter, meant to abstract away common tasks.
"""

import inspect

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class BaseFilter:
    """Abstract class representing a filter for processing candidate probes.

    A filter may operate in the more traditional sense of reducing the
    input probes to a smaller number (i.e., only allowing certain probes
    to pass through); an example would be a filter the eliminates
    duplicates. A filter may also operate in the more computing sense
    of simply processing the input probes and returning the list of
    processed probes; the output/processed probes may be altered
    versions of the input or there may even be more output probes than
    input probes.

    All subclasses must implement a _filter(..) method that returns a
    list of probes after processing from the given input list.
    """

    def filter(self, input, target_genomes=None):
        """Perform the filtering.

        Args:
            input: list of candidate probes
            target_genomes: list [g_1, g_2, g_m] of m groupings of genomes,
                where each g_i is a list of genome.Genomes belonging to group
                i, that should be targeted by the probes; for example a
                group may be a species and each g_i would be a list of the
                target genomes of species i

        Returns:
            list of probes after applying a filter to the input
        """
        _filter_params = inspect.signature(self._filter).parameters
        if len(_filter_params) == 2:
            # _filter() should accept both probes and target genomes
            return self._filter(input, target_genomes)
        else:
            # _filter() may not need target genomes, and does not accept it
            return self._filter(input)

    def _filter(self, input):
        raise Exception(("A subclass of BaseFilter must implement "
                         "_filter(..)"))
