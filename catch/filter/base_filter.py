"""Wrappers around a filter, meant to abstract away common tasks.
"""

import inspect
import multiprocessing

from catch.utils import fix_spawn_behavior

__author__ = 'Hayden Metsky <hayden@mit.edu>'


def set_max_num_processes_for_filter_over_groupings(max_num_processes=8):
    """Set the maximum number of processes to use for parallelizing calls
    to _filter() across groupings.

    Note that parallelization defined in this module does not always occur.
    See the `num_processes` arg in BaseFilter.filter() for when it does
    apply.

    Args:
        max_num_processes: an int (>= 1) specifying the maximum number of
            processes to use in a multiprocessing.Pool when parallelizing
            over groupings, i.e., the maximum number of target groupings
            to filter in parallel; it uses min(the number of CPUs
            in the system, max_num_processes) processes
    """
    global _fg_max_num_processes
    _fg_max_num_processes = max_num_processes
set_max_num_processes_for_filter_over_groupings()

# Define filter function to use in multiprocessing Pool; this must be
# top-level in the module, and we will ensure only one can be set at a time
global _global_filter_fn
_global_filter_fn = None


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

    For information about parallelization over groupings, see the
    `num_processes` argument below.

    All subclasses must implement a _filter(..) method that returns a
    list of probes after processing from the given input list.
    """

    def filter(self, input, target_genomes=None, input_is_grouped=False,
            num_processes=None):
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
            num_processes: number of processes to use when parallelizing over
                groupings; if None, this determines a number based on the
                maximum specified and the number of CPUs. Note that
                parallelization only happens when input_is_grouped is True
                *and* self.requires_probe_groupings is not set or is False; if
                that parameter is True and input_is_grouped is True, then
                all groupings are passed to the subclass's filter and it
                is up to self._filter() to parallelize over groupings

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
                # Call _filter() separately for each group, and parallelize
                # calls across groupings

                fix_spawn_behavior.fix_spawn_behavior()

                global _fg_max_num_processes
                if num_processes is None:
                    num_processes = min(multiprocessing.cpu_count(),
                                        _fg_max_num_processes)
                pool = multiprocessing.Pool(num_processes)

                # Order groupings in descending order
                #   by the number of possible probes (input size) in the group.
                #   The number is an indication of how long the grouping may
                #   take to filter, and we want to start the slower groupings
                #   first in the pool
                input_lens = list(enumerate([len(x) for x in input]))
                input_idx_ordered = [x[0] for x in sorted(input_lens,
                    key=lambda y: y[1], reverse=True)]
                input_idx_revert = {y: x for x, y in
                        enumerate(input_idx_ordered)}
                # Note that the reordered input is:
                #   [input[i] for i in input_idx_ordered]

                # The function called by a multiprocessing Pool must be
                #   top-level
                global _global_filter_fn
                if _global_filter_fn is not None:
                    raise Exception(("Only one filter() function can be "
                        "called in parallel at a time"))
                _global_filter_fn = self._filter

                # Construct args to _filter()
                if len(_filter_params) == 2:
                    # self._filter() should accept both probes and target genomes
                    pool_args = [(input[i], target_genomes)
                            for i in input_idx_ordered]
                else:
                    # self._filter() may not need target genomes, and does not
                    # accept it
                    pool_args = [tuple([input[i]])
                            for i in input_idx_ordered]

                # Run the pool, giving 1 grouping (chunksize=1) at a time
                pool_out = pool.starmap(_global_filter_fn, pool_args,
                        chunksize=1)
                pool.close()
                _global_filter_fn = None

                # Revert the order of the output to go back to the original
                #   ordering of the input
                pool_out_reordered = [pool_out[input_idx_revert[i]]
                        for i in range(len(pool_out))]
                return pool_out_reordered
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
