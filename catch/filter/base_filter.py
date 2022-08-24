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

    def filter(self, input, target_genomes=None, input_is_grouped=False):
        """Perform the filtering.

        Args:
            input: candidate probes from which to filter; see input_is_grouped
                for details
            target_genomes: list [g_1, g_2, ..., g_m] of m groupings of genomes,
                where each g_i is a list of genome.Genomes belonging to group
                i, that should be targeted by the probes; for example a
                group may be a species and each g_i would be a list of the
                target genomes of species i
            input_is_grouped: if True, input is list [p_1, p_2, ..., p_m] of
                m groupings of genomes, where each p_i is a list of candidate
                probes for group i; if False, input is a single list of
                candidate probes (ungrouped)

        Returns:
            if input_is_grouped is True:
                list [q_1, q_2, q_m] where each q_i is a list of probes after
                applying a filter to the corresponding input
            else:
                list of probes after applying a filter to the input
        """
        _filter_params = inspect.signature(self._filter).parameters

        # Determine whether self._filter() requires probes being
        #   split into groupings, or whether each group must be passed
        #   separately
        if (hasattr(self, 'requires_probe_groupings') and
                self.requires_probe_groupings is True):
            pass_groupings = True
        else:
            pass_groupings = False

        if pass_groupings:
            # Input must already be grouped
            assert input_is_grouped is True

            if len(_filter_params) == 2:
                # self._filter() should accept both probes and target genomes
                return self._filter(input, target_genomes)
            else:
                # self._filter() may not need target genomes, and does not
                # accept it
                return self._filter(input)
        else:
            if input_is_grouped:
                # Call _filter() separately for each group

                if len(_filter_params) == 2:
                    # self._filter() should accept both probes and target genomes
                    return [self._filter(i, target_genomes) for i in input]
                else:
                    # self._filter() may not need target genomes, and does not
                    # accept it
                    return [self._filter(i) for i in input]
            else:
                # Input is not grouped and there is no need to pass it grouped

                if len(_filter_params) == 2:
                    # self._filter() should accept both probes and target genomes
                    return self._filter(input, target_genomes)
                else:
                    # self._filter() may not need target genomes, and does not
                    # accept it
                    return self._filter(input)

    def _filter(self, input):
        raise Exception(("A subclass of BaseFilter must implement "
                         "_filter(..)"))
