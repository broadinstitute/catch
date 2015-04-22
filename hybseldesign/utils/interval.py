"""Structures and functions for working with intervals.
"""

import bisect

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class IntervalSet:
    """Immutable collection of intervals.

    Each interval is a tuple of the form (start, end) where start
    is inclusive and end is exclusive. The intervals that are
    stored by this structure are sorted and non-overlapping.
    """

    def __init__(self, intervals):
        """
        Args:
            intervals: collection of intervals, each of the form
                    (start, end) where start is inclusive and
                    end is exclusive
        """
        # Save intervals in a sorted, non-overlapping tuple
        self.intervals = tuple(merge_overlapping(list(intervals)))
        # Save starting endpoints, which is useful for binary search over
        # the intervals
        self.starting_endpoints = [x[0] for x in self.intervals]
        # Find the very beginning and end of these intervals
        if len(self.intervals) > 0:
            self.first_start = self.intervals[0][0]
            self.last_end = self.intervals[-1][1]
        else:
            self.first_start = None
            self.last_end = None
        self.len_cached = None

    def _merge(self, other, include_in_merge, optimize_intersection=False):
        """Merge this interval set with another.

        Args:
            other: instance of IntervalSet
            include_in_merge: function f(x,y) such that x is a boolean
                indicating whether some element is contained within an
                interval in self.intervals and y is a boolean indicating
                whether that same element is contained within an interval
                in other.intervals; f(x,y) returns whether that element
                is to be contained within an interval in the merge
            optimize_intersection: When True, makes optimizations specific
                for performing an intersection

        Returns:
            list of intervals (tuples) such that those in self.intervals
            and other.intervals are merged according to the function
            include_in_merge
        """
        merge_intervals = []
        curr_merge_interval_start = None

        self_endpoint_index = 0
        other_endpoint_index = 0

        if optimize_intersection:
            if len(self.intervals) == 0 or len(other.intervals) == 0:
                return []
            if len(self.intervals) < len(other.intervals):
                smaller, larger = self, other
            else:
                smaller, larger = other, self
            # Start considering intervals in larger.intervals at an index
            # determined by the first interval in smaller.intervals
            # This is especially helpful when smaller.intervals contains
            # just one interval
            first_interval = smaller.intervals[0]
            insert_point = bisect.bisect_left(larger.starting_endpoints,
                                              first_interval[0])
            # Start looking at larger.intervals at its (insert_point - 1)'th
            # interval, since the first interval in smaller.intervals falls
            # in this interval or before it
            insert_point = max(0, insert_point - 1)
            # Multiply by 2 because the endpoint indices correspond to
            # endpoints of intervals, and there are two per interval
            if len(self.intervals) < len(other.intervals):
                other_endpoint_index = 2 * insert_point
            else:
                self_endpoint_index = 2 * insert_point

        # Scan through all the endpoints of the intervals in self.intervals
        # and in other.intervals simultaneously from left to right
        # (an interval has two endpoints)
        intervals_left_to_scan = True
        while intervals_left_to_scan:
            # The endpoint indices correspond to endpoints of intervals (two
            # per interval), so use these indices to fetch the corresponding
            # interval. Then determine whether curr_endpoint is contained
            # within each of the "current" intervals.
            if self_endpoint_index < 2 * len(self.intervals):
                curr_self_interval = self.intervals[self_endpoint_index / 2]
            else:
                curr_self_interval = None
            if other_endpoint_index < 2 * len(other.intervals):
                curr_other_interval = other.intervals[other_endpoint_index / 2]
            else:
                curr_other_interval = None

            if not curr_self_interval and not curr_other_interval:
                # There are no more intervals left to consider
                intervals_left_to_scan = False
            elif curr_self_interval and not curr_other_interval:
                if optimize_intersection:
                    # No further intersections are possible
                    intervals_left_to_scan = False
                # Consider the current endpoint of curr_self_interval
                curr_endpoint = curr_self_interval[self_endpoint_index % 2]
            elif not curr_self_interval and curr_other_interval:
                if optimize_intersection:
                    # No further intersections are possible
                    intervals_left_to_scan = False
                # Consider the current endpoint of curr_other_interval
                curr_endpoint = curr_other_interval[other_endpoint_index % 2]
            else:
                # Consider the smaller of the current endpoints of the two
                # intervals
                self_endpoint = curr_self_interval[self_endpoint_index % 2]
                other_endpoint = curr_other_interval[other_endpoint_index % 2]
                curr_endpoint = min(self_endpoint, other_endpoint)

            # Determine whether curr_endpoint is contained within each of
            # the "current" intervals
            if curr_self_interval:
                within_curr_self = (curr_endpoint >= curr_self_interval[0] and
                                    curr_endpoint < curr_self_interval[1])
            else:
                within_curr_self = False
            if curr_other_interval:
                within_other_self = (curr_endpoint >= curr_other_interval[0] and
                                     curr_endpoint < curr_other_interval[1])
            else:
                within_other_self = False

            # Determine whether curr_endpoint should be within an interval
            # in the merge and update merge_intervals accordingly
            if include_in_merge(within_curr_self, within_other_self):
                if curr_merge_interval_start is None:
                    # Explicitly check 'is None' in case
                    #  curr_merge_interval_start is 0
                    # curr_endpoint should be in a merge interval but we
                    # are not currently within the middle of a merge interval,
                    # so start one here
                    curr_merge_interval_start = curr_endpoint
                # Else, curr_endpoint should be in a merge interval but we
                # are already within the middle of one, so do nothing
            else:
                if curr_merge_interval_start is not None:
                    # Explicitly check 'is not None' in case
                    #  curr_merge_interval_start is 0
                    # curr_endpoint should not be in a merge interval but
                    # we currently are within the middle of a merge interval,
                    # so end it here
                    merge_intervals += [(curr_merge_interval_start,
                                         curr_endpoint)]
                    curr_merge_interval_start = None
                # Else, curr_endpoint should not be in a merge interval and
                # we are not currently within the middle of one, so do
                # nothing

            if (curr_self_interval and
                    (curr_endpoint == curr_self_interval[0] or
                     curr_endpoint == curr_self_interval[1])):
                # curr_endpoint is an endpoint of the "current" interval from
                # self.intervals, so advance the index for this
                self_endpoint_index += 1
            if (curr_other_interval and
                    (curr_endpoint == curr_other_interval[0] or
                     curr_endpoint == curr_other_interval[1])):
                # curr_endpoint is an endpoint of the "current" interval from
                # other.intervals, so advance the index for this
                other_endpoint_index += 1

        return merge_intervals

    def intersection(self, other):
        """Make the intersection of this interval set and other.

        Args:
            other: instance of IntervalSet

        Returns:
            instance of IntervalSet whose intervals consist of the
            intersection of self.intervals with other.intervals
        """
        def include_in_merge(within_self, within_other):
            return within_self and within_other
        return IntervalSet(self._merge(other, include_in_merge,
                                       optimize_intersection=True))

    def union(self, other):
        """Make the union of this interval set and other.

        Args:
            other: instance of IntervalSet

        Returns:
            instance of interval set whose intervals consist of the
            union of self.intervals with other.intervals
        """
        def include_in_merge(within_self, within_other):
            return within_self or within_other
        return IntervalSet(self._merge(other, include_in_merge))

    def difference(self, other):
        """Make the difference of this interval set with other.

        Args:
            other: instance of IntervalSet

        Returns:
            instance of IntervalSet whose intervals contain all elements
            within intervals of self.intervals but not within intervals
            of other.intervals (i.e., self.intervals - other.intervals
            where the minus operation acts on elements contained within
            intervals)
        """
        def include_in_merge(within_self, within_other):
            return within_self and not within_other
        return IntervalSet(self._merge(other, include_in_merge))

    def __len__(self):
        """Count the number of integer elements contained within all intervals.

        This assumes that interval endpoints are integers and that an
        element is an integer. For example, the interval (3,6) contains
        three elements (3, 4, and 5). The collection of intervals
        [(3,6), (8,10)] contains five elements (3, 4, 5, 8, and 9).

        Returns:
            number of integer elements contained within all intervals
        """
        if self.len_cached is None:
            self.len_cached = sum(x[1] - x[0] for x in self.intervals)
        return self.len_cached

    def __hash__(self):
        return hash(self.intervals)

    def __eq__(self, other):
        return isinstance(other, IntervalSet) and \
            self.intervals == other.intervals

    def __str__(self):
        return str(self.intervals)

    def __repr__(self):
        return str(self.intervals)


def merge_overlapping(intervals):
    """Merge a list of possibly overlapping intervals.

    Args:
        intervals: list of intervals, each of which is a tuple of the
            form (start, end); start is inclusive, end is exclusive

    Returns:
        list of intervals in which overlapping ones are merged. The
        list of intervals is sorted. For example, the input
        [(1,5), (3,7), (9,12)] yields [(1,7), (9,12)]. Intervals that
        are touching (e.g., (1,3) and (3,5)) are merged into one.
    """
    if len(intervals) == 0:
        return []

    intervals = sorted(intervals)
    intervals_merged = []
    curr_start, curr_end = intervals[0][0], intervals[0][1]
    for start, end in intervals:
        if start <= curr_end:
            curr_end = max(curr_end, end)
        else:
            intervals_merged += [(curr_start, curr_end)]
            curr_start = start
            curr_end = end
    intervals_merged += [(curr_start, curr_end)]

    return intervals_merged


def schedule(intervals):
    """Schedule the maximum number of compatible intervals.

    This uses the greedy interval scheduling algorithm to find (schedule)
    the maximum number of compatible (non-overlapping) intervals.

    Args:
        intervals: list of intervals, each of which is a tuple (x,y)
            in which x is a tuple of the form (start, end) and y is a
            reference to an object represented by the interval.

    Returns:
        list of the objects corresponding to the chosen intervals
        (i.e., the 'y' for each chosen element)
    """
    # Sort all intervals by their endpoint (the "finishing time")
    # x[0] gives the interval and x[0][1] gives the endpoint
    intervals = sorted(intervals, key=lambda x: x[0][1])

    # Scan through the intervals in sorted order and choose
    # compatible ones with the earliest endpoint
    last_chosen_interval = None
    chosen_objects = []
    for interval, obj in intervals:
        is_compatible = False
        if last_chosen_interval is None:
            # No intervals have been chosen yet, so interval is
            # of course compatible
            is_compatible = True
        else:
            # interval is compatible with the chosen intervals iff
            # its start comes after the finish of the last chosen
            # interval
            if interval[0] >= last_chosen_interval[1]:
                is_compatible = True
        if is_compatible:
            last_chosen_interval = interval
            chosen_objects += [obj]

    return chosen_objects
