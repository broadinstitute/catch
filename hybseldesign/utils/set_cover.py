"""Functions for working with instances of the set cover problem.
"""

from collections import defaultdict
import logging

from hybseldesign.utils import interval

__author__ = 'Hayden Metsky <hayden@mit.edu>'

logger = logging.getLogger(__name__)


def approx(sets, costs=None, p=1.0):
    """Approximates the solution to an instance of the set cover problem.

    This solves the problem using the well-known greedy algorithm.
    The following is a description of the problem and solution to the
    most basic case, in which sets are unweighted and we seek to cover
    the entire universe:
    We are given are universe U of (hashable) objects and a collection
    of m subsets S_1, S_2, ..., S_m of U whose union equals U. We wish
    to approximate the smallest number of these subets whose union is U
    (i.e., "covers" the universe U). Pseudocode of the greedy algorithm
    is as follows:
      C <- {}
      while the universe is not covered (U =/= {}):
        Pick the set S_i that covers the most of U (i.e., maximizes
          | S_i _intersection_ U |)
        C <- C _union_ { S_i }
        U <- U - S_i
      return C
    The collection of subsets C is the approximate set cover and a
    valid set cover is always returned. The loop goes through min(|U|,m)
    iterations and picking S_i at each iteration takes O(m*D) time where
    D is the cardinality of the largest subset. The returned solution
    is a ceil(ln(D))-approximation. In the worst-case, where D=|U|, this
    is a ceil(ln(|U|))-approximation. Inapproximability results show
    that it is NP-hard to approximate the problem to within c*ln(|U|)
    for any 0 < c < 1. Thus, this is a very good approximation given
    what is possible.

    This generalizes the above case to support weighted, partial set
    cover. There are two primary changes to the above pseudocode that
    support this:
     - Rather than looping while the universe is not covered, we instead
       loop until we have covered as much as the universe as we desire
       to cover (where that amount is determined by the parameter
       indicating the fraction of the universe to cover). This change
       alone supports 'partial cover'.
     - Rather than picking the set S_i that covers the most of U, we
       instead pick the set S_i that minimizes the ratio of the cost
       of S_i (call it c_i) to the amount of the remaining universe, U,
       that S_i covers. Let r be the number of elements that have yet
       to be covered in order to obtain the desired partial cover; there
       is no reason to favor a set that covers more than r. Thus, in
       particular, we pick the set S_i that minimizes the quotient
       [ c_i / min(r, | S_i \intersect U |) ].
    As shown by Petr Slavik in the paper "Improved performance of the
    greedy algorithm for the minimum set cover and minimum partial cover
    problems", the obtained solution is a ceil(ln(D))-approximation,
    where D is the cardinality of the largest subset. It also obtains a
    ceil(ln(p*|U|))-approximation where p is the fraction of the universe
    we wish to cover. When p=1, D <= |U| and therefore the
    ceil(ln(D))-approximation is the stronger one. But when p<1, it is
    possible that D > p*|U| and hence the ceil(ln(p*|U|))-approximation
    would be stronger. Note that in the weighted case, the solution
    seeks to minimize the sum of the weights of the chosen sets and
    the approximation is with respect to this.

    Args:
        sets: dict mapping set identifiers to sets
        costs: dict mapping set identifiers to the costs (or weights)
            of the set; the default is for every set to have a cost
            of 1, making this equivalent to the unweighted problem
        p: float in [0,1] that specifies the fraction of the universe
            we must cover; the default is p=1, making this equivalent
            to the problem in which the entire universe is covered

    Returns:
        a set consisting of the identifiers of the sets chosen to be
        in the set cover. For example, if the sets input is
          { 0: S_1, 1: S_2, 2: S_3 }
        where each S_i is a set, and the sets S_1 and S_3 are
        chosen to be in the set cover, then the output is {0,2}.
    """
    if p < 0 or p > 1:
        raise ValueError("p must be in [0,1]")
    if costs is None:
        # Give each set a default cost of 1
        costs = {set_id: 1 for set_id in sets.keys()}
    else:
        for c in costs.values():
            if c < 0:
                raise ValueError("All costs must be nonnegative")

    # Create the universe as the union of all input sets
    universe = set()
    for s in sets.values():
        universe.update(s)

    num_that_can_be_uncovered = int(len(universe) - p * len(universe))
    # Above, use int(..) to take the floor. Also, expand out
    # len(universe) rather than use int((1.0-p)*len(universe)) due to
    # precision errors in Python -- e.g., int((1.0-0.8)*5) yields 0 on
    # some machines.
    num_left_to_cover = len(universe) - num_that_can_be_uncovered

    set_ids_not_in_cover = set(sets.keys())
    set_ids_in_cover = set()
    # Keep iterating until desired partial cover is obtained
    while num_left_to_cover > 0:
        # Find the set that minimizes the ratio of its cost to the
        # number of uncovered elements (that need to be covered) that
        # it covers
        id_min_ratio, min_ratio = None, float('inf')
        for id in set_ids_not_in_cover:
            # There is no strict need to keep track of set_ids_not_in_cover
            # and iterate through these; we could have iterated over
            # all input sets. However, sets that are already put into the
            # set cover have zero intersection with universe (i.e., cover
            # none of it), so there is no reason to iterate over the sets
            # already placed in the cover. Not doing so should yield some
            # runtime improvements.
            s = sets[id]
            num_covered = len(s.intersection(universe))
            # There is no need to cover more than num_left_to_cover
            # elements
            num_needed_covered = min(num_left_to_cover, num_covered)
            if num_needed_covered == 0:
                # s covers no elements that need to be covered, so it should
                # not be in the set cover
                continue
            ratio = float(costs[id]) / num_needed_covered
            if ratio < min_ratio:
                id_min_ratio = id
                min_ratio = ratio
        # id_min_ratio goes into the set cover
        set_ids_in_cover.add(id_min_ratio)
        set_ids_not_in_cover.remove(id_min_ratio)
        universe.difference_update(sets[id_min_ratio])
        num_left_to_cover = max(0, len(universe) - num_that_can_be_uncovered)

    return set_ids_in_cover


def approx_multiuniverse(sets,
                         costs=None,
                         universe_p=None,
                         ranks=None,
                         use_arrays=False,
                         use_intervalsets=False):
    """Approximates the solution to a "multiuniverse" set problem.

    We define the "multiuniverse" set problem to be a version of the
    set cover problem in which there are multiple universes and we
    seek to find a collection of sets whose union covers a specified
    fraction of each given universe. The sets can contain elements
    from different universes; the input sets provide, for each set,
    which elements in which universes the set covers.

    Note that elements with the same value from different universes are
    effectively treated as different elements. For example, covering the
    element "42" from universe 2 does not necessarily cover the element
    "42" from universe 1 (unless a chosen set contains "42" from
    universe 2 as well as "42" from universe 1).

    Args:
        sets: dict mapping set identifiers to other dicts, which then
            map universe identifiers to sets of elements. That is,
            consider all the elements in some set set_id. sets[set_id]
            is a dict in which all the elements in set_id are split up
            by the universe they are a part of;
            sets[set_id][universe_id] is a set containing the elements
            in set_id that come from universe universe_id
        costs: dict mapping set identifiers to the costs (or weights)
            of the set; the default is for every set to have a cost of
            1, making this equivalent to the unweighted problem
        universe_p: dict mapping universe identifiers to floats in
            [0,1] that specifiy the fraction of the corresponding
            universe we must cover; the default is for the coverage
            fraction of each universe to be 1.0, making this equivalent
            to the problem in which there is a single universe (the
            union of all given universes) that must be completely
            covered
        ranks: dict mapping set identifiers to a rank (integer) for the
            set. The rank of a set makes it easy to define different
            levels of penalties on sets. If there are two sets A and B
            such that ranks[A] < ranks[B], then A will always be
            considered before B -- i.e., if coverage is needed and A
            provides that coverage, A will be chosen before B even if A
            provides less coverage than B would provide. When two sets
            have the same rank, then the costs of those sets are
            considered. In a sense, the ranks of two sets is given
            higher priority than the costs of those sets or the number
            of elements they cover; the ranks define different groups
            such that as much coverage as possible is sought from sets
            with lesser rank before considering sets with higher rank.
            [See implementation note below.]
        use_arrays: when True, the values inside the input 'sets' are
            expected to be stored in a Python array rather than in a
            Python set. Python sets require a considerable amount of
            overhead and, though the use of array is less time-efficient
            here, it can lead to substantial memory savings depending on
            the type of data.
        use_intervalsets: when True, the values inside the input 'sets'
            are expected to be stored as an instance of IntervalSet
            (one instance per set). For sets that contain many adjacent
            elements (i.e., stretches of contiguous integers), this
            should require less memory and less runtime than the use of
            Python sets or arrays. This cannot be true when 'use_arrays'
            is also True. This should only be used when the elements
            (usually integers) are well represented by intervals --
            i.e., they fit into contiguous stretches -- because otherwise
            this will run slowly compared to when the option is not set.
            If this is true and the value inside the input 'sets' at some
            entry is a tuple, then it is assumed that the set at that entry
            has just one interval (i.e., the one specified by the tuple),
            which is useful for saving space, and the interval is converted
            into an instance of IntervalSet as needed.

    Returns:
        a set consisting of the identifiers of the sets chosen to be
        in the set cover. For example, if the sets input is
          { 0: S_1, 1: S_2, 2: S_3 }
        where each S_i is a set, and the sets S_1 and S_3 are chosen
        to be in the set cover, then the output is {0,2}.

    Implementation notes:
      - The rank of a set can thought of as giving a sufficiently high
        cost to the set.  In fact, ranks are not at all necessary, as
        they can be emulated entirely using costs. We choose some
        sufficiently large constant C; C should be at least as large as
        the total number of elements across all the universes (e.g.,
        2^64).  Then, if we wish to assign some set a nonnegative rank
        r, we could instead take the cost of the set and multiply that
        cost by C^r prior to calling this function. That would have the
        same effect of delaying the consideration of the set until all
        sets with smaller ranks have been considered. However, while
        correct in theory, the method of using costs to emulate ranks is
        not practical. For example, if C=2^64 and r=100 for some set,
        then the cost of that set with r=100 would be its previous cost
        multiplied by (2^64)^100. While Python can compute this number,
        it cannot do floating point arithmetic with it and thus cannot
        consider its ratio in the weighted set cover approximation
        algorithm.
      - This is a generalization of the partial, weighted set cover
        problem whose solution is approximated by approx(..). (For any
        input to that function, simply create a single universe (e.g.,
        with id 0) whose elements consist of the union of all the input
        sets, and then create a new input sets to this function in which
        each original set set_id is instead found in sets[set_id][0].)
        Indeed, this function and approx(..) are very similar, with this
        one simply generalizing to allow for multiple universes.
        However, whereas approx(..) achieves a provable guarantee on the
        approximation factor, it is not clear whether the
        straightforward generalization implemented in this function
        achieves the same factor; therefore, there may be room for
        improvement in this function. Furthermore, the generalization to
        multiple universes will yield a slight increase in runtime
        (constant factors) compared to working directly with the more
        traditional notion of a single universe. For these reasons --
        although the approx(..) function could be greatly shortened by
        simply calling this function -- we choose to implement the
        versions separately.
      - In practice, often the computed minimum ratio in one iteration of
        the set cover algorithm is equal to the minimum ratio from the
        previous iteration. We exploit this fact to improve the running
        time of the greedy algorithm without affecting the output. We store
        the minimum ratio from the previous iteration (last_min_ratio) as
        well as all the sets whose ratio, as computed in that previous
        iteration, was equal to last_min_ratio
        (set_ids_with_same_ratio_as_last_min). On each iteration, we first
        check if there is a set whose ratio is equal to the last minimum
        ratio and is valid for being put into the set cover. Because the
        minimum ratio is nondecreasing as we step through iterations, such
        a set would have the minimum ratio on this iteration; thus, we
        simply add it to the set cover and progress to the next iteration
        of the greedy algorithm. If there is no such set, we must take
        the usual approach that would not use this heuristic -- i.e.,
        compute the ratios for all sets not yet in the set cover and find
        the set with the minimum ratio.
    """
    if use_arrays and use_intervalsets:
        raise ValueError("Cannot use both arrays and IntervalSets")

    if costs is None:
        # Give each set a default cost of 1
        costs = {set_id: 1 for set_id in sets.keys()}
    else:
        for c in costs.values():
            if c < 0:
                raise ValueError("All costs must be nonnegative")
        for set_id in sets.keys():
            if set_id not in costs:
                raise ValueError("costs is missing a value for set %d" %
                                 set_id)

    # Create the universes from given sets
    if use_intervalsets:
        # Store the elements of each universe in an IntervalSet
        universes = defaultdict(lambda: interval.IntervalSet([]))
    else:
        # Store the elements of each universe in a set
        universes = defaultdict(set)
    for sets_by_universe in sets.values():
        for universe_id, s in sets_by_universe.items():
            if use_intervalsets:
                if isinstance(s, tuple):
                    # s is a single interval
                    s = interval.IntervalSet([s])
                universes[universe_id] = universes[universe_id].union(s)
            elif use_arrays:
                for v in s:
                    universes[universe_id].add(v)
            else:
                universes[universe_id].update(s)
    universes = dict(universes)

    if universe_p is None:
        # Give each universe a coverage fraction of 1.0 (i.e., cover
        # all of it)
        universe_p = {universe_id: 1 for universe_id in universes.keys()}
    else:
        for p in universe_p.values():
            if p < 0 or p > 1:
                raise ValueError(("The coverage fraction (p) of each "
                                  "universe must be in [0,1]"))
        for universe_id in universes.keys():
            if universe_id not in universe_p:
                raise ValueError(("universe_p is missing a value for "
                                  "universe %d" % universe_id))

    if ranks is None:
        # Give each set a default rank of 1; since all sets have the
        # same rank, this is effectively the same as not using ranks
        # at all
        ranks = {set_id: 1 for set_id in sets.keys()}
    else:
        for set_id in sets.keys():
            if set_id not in ranks:
                raise ValueError("ranks is missing a value for set %d" %
                                 set_id)
    # Track the current index (curr_rank_index) as we step through all
    # of the (distinct) ranks (rank_vals)
    rank_vals = sorted(set(ranks.values()))
    curr_rank_index = 0

    num_that_can_be_uncovered = {}
    num_left_to_cover = {}
    for universe_id, universe in universes.items():
        p = universe_p[universe_id]
        num_that_can_be_uncovered[universe_id] = \
            int(len(universe) - p * len(universe))
        # Above, use int(..) to take the floor. Also, expand out
        # len(universe) rather than use int((1.0-p)*len(universe)) due to
        # precision errors in Python -- e.g., int((1.0-0.8)*5) yields 0 on
        # some machines.
        num_left_to_cover[universe_id] = \
            len(universe) - num_that_can_be_uncovered[universe_id]

    # Computing intersections between a particular set and a particular
    # universe is the bottleneck of this function. However, a lot of
    # time the intersection may be computed even though it has already
    # been computed between the same set and (unchanged) universe.
    # (While the contents of sets are not modified, elements in
    # universes are discarded. But not every universe is modified at
    # every iteration, so when one is not modified an intersection with
    # that universe is unnecessarily computed again.) To avoid this,
    # memoize the sizes of the intersections of universes with sets
    # and discard memoized values for a particular universe U_i when
    # U_i is updated. When using interval sets, it is possible to further
    # optimize this so that, on an update of U_i, the only values that
    # are discarded are ones that might change due to the update.
    memoized_intersect_counts = {
        universe_id: {}
        for universe_id in universes.keys()
    }

    def compute_ratio_for_set(set_id):
        # By iterating over all universes covered by set_id, compute the ratio
        # for set_id that will be used to determine whether it should be placed
        # in the set cover
        num_needed_covered_across_universes = 0
        for universe_id in sets[set_id].keys():
            if set_id in memoized_intersect_counts[universe_id]:
                # We have num_covered memoized
                num_covered = memoized_intersect_counts[universe_id][set_id]
            else:
                s = sets[set_id][universe_id]
                universe = universes[universe_id]
                if use_arrays:
                    # It may seem faster to compute, in the case where s is
                    # an array, num_covered as
                    # sum([1 for v in s if v in universe])
                    # in order to avoid converting s to a set. However,
                    # it appears that, in practice, converting s to a set
                    # and using set.intersection is faster.
                    s = set(s)
                if use_intervalsets and isinstance(s, tuple):
                    # s is a single interval
                    s = interval.IntervalSet([s])
                # If use_intervalsets, then s and universe should already
                # be IntervalSets, and the intersection method is defined
                # for these
                num_covered = len(s.intersection(universe))
                # Memoize num_covered
                memoized_intersect_counts[universe_id][set_id] = num_covered
            # There is no need to cover more than num_left_to_cover
            # elements for this universe
            num_needed_covered = min(num_left_to_cover[universe_id],
                                     num_covered)
            num_needed_covered_across_universes += num_needed_covered
        if num_needed_covered_across_universes == 0:
            # s covers no elements that need to be covered, so it should
            # not be in the set cover; give it an infinite ratio
            ratio = float('inf')
        else:
            ratio = float(costs[set_id]) / num_needed_covered_across_universes
        return ratio

    # Track the computed ratio of the set that was last added to the set cover
    last_min_ratio = None
    # When ratios must be computed to find the minimum, this is done across
    # all sets not in the set cover; track the sets whose computed ratio is
    # equal to last_min_ratio. Note that this collection does not include the
    # set with ratio last_min_ratio that was initially put into the set
    # cover when last_min_ratio was first computed.
    set_ids_with_same_ratio_as_last_min = []

    set_ids_not_in_cover = set(sets.keys())
    set_ids_in_cover = set()
    # Keep iterating until desired partial cover of each universe
    # is obtained (note that [] evaluates to False)
    while [True for universe_id in universes.keys()
           if num_left_to_cover[universe_id] > 0]:
        if len(set_ids_in_cover) % 10 == 0:
            logger.info(("Selected %d sets with a total of %d elements "
                         "remaining to be covered"), len(set_ids_in_cover),
                        sum(num_left_to_cover.values()))

        # Find the set that minimizes the ratio of its cost to the
        # number of uncovered elements (that need to be covered) that
        # it covers
        id_min_ratio = None

        # First, look among all sets whose ratio equals the last minimum
        # ratio. Because the minimum ratio is nondecreasing across iterations,
        # if one set's ratio equals the last minimum ratio, this one must also
        # be a minimum on this iteration.
        for set_id in set_ids_with_same_ratio_as_last_min:
            # Check that set_id was not yet chosen (i.e., is in
            # set_ids_not_in_cover); it may be the case that
            # set_ids_with_same_ratio_as_last_min was the same on a previous
            # iteration and set_id was chosen to be id_min_ratio on that
            # previous iteration.
            # Also, re-compute the ratio for set_id and check that it is
            # still equal to last_min_ratio; it may be the case that, on a
            # previous iteration, choosing a set with ratio last_min_ratio
            # altered the ratio of set_id. For example, that set may cover
            # some of the same elements as set_id covers, which would have
            # caused an invalidation to memoized_intersect_counts and would
            # increase the ratio for set_id.
            if (set_id in set_ids_not_in_cover and
                    compute_ratio_for_set(set_id) == last_min_ratio):
                id_min_ratio = set_id
                break

        if id_min_ratio == None:
            # The above heuristic -- which looks for sets whose ratio equals
            # last_min_ratio -- failed to find a set for this iteration. So
            # simply iterate over all sets in set_ids_not_in_cover, compute
            # the ratio for each set, and find the one with the minimum ratio.
            min_ratio = float('inf')
            for set_id in set_ids_not_in_cover:
                # There is no strict need to keep track of set_ids_not_in_cover
                # and iterate through these; we could have iterated over
                # all input sets. However, sets that are already put into the
                # set cover have zero intersection with any universe (i.e.,
                # cover none of it), so there is no reason to iterate over the
                # sets already placed in the cover. Not doing so should yield
                # some runtime improvements.
                if ranks[set_id] != rank_vals[curr_rank_index]:
                    # Skip this set because its rank is not the current rank
                    # being considered.
                    # We could (correctly) use '>' instead of '!='. But doing so
                    # would also consider any sets whose rank is less than the
                    # current rank being considered (rank_vals[curr_rank_index]).
                    # Because the rank being considered is strictly increasing,
                    # these sets should have already been considered at a
                    # previous point. Since the rank increased without these
                    # sets having been put into the cover, it must have been the
                    # case that they did not cover any elements that needed to
                    # be covered; this would still hold true for these sets, so
                    # there is no reason to consider them again.
                    continue
                ratio = compute_ratio_for_set(set_id)
                if ratio < min_ratio:
                    id_min_ratio = set_id
                    min_ratio = ratio
                    # Since there is a new min_ratio, reset the collection
                    # of sets whose ratio equals last_min_ratio
                    set_ids_with_same_ratio_as_last_min = []
                elif ratio == min_ratio:
                    set_ids_with_same_ratio_as_last_min += [set_id]
            last_min_ratio = min_ratio

        if id_min_ratio is None:
            # Increase the rank being considered and try again
            curr_rank_index += 1
            set_ids_with_same_ratio_as_last_min = []
            continue

        # id_min_ratio goes into the set cover
        set_ids_in_cover.add(id_min_ratio)
        set_ids_not_in_cover.remove(id_min_ratio)
        for universe_id, universe in universes.items():
            if universe_id not in sets[id_min_ratio]:
                # id_min_ratio covers nothing in this universe
                continue
            s = sets[id_min_ratio][universe_id]
            prev_universe_size = len(universe)
            # Remove s from universe
            if use_intervalsets:
                if isinstance(s, tuple):
                    # s is a single interval
                    s = interval.IntervalSet([s])
                universe = universe.difference(s)
                universes[universe_id] = universe
            elif use_arrays:
                for v in s:
                    universe.discard(v)
            else:
                universe.difference_update(s)
            num_left_to_cover[universe_id] = max(
                0, len(universe) - num_that_can_be_uncovered[universe_id])
            # Discard memoized values
            if len(universe) != prev_universe_size:
                if use_intervalsets:
                    # The universe was modified and since we are using interval
                    # sets we can optimize what values we choose to discard
                    # (i.e., invalidate). In particular, only invalidate
                    # sets whose interval set overlaps the interval set
                    # deleted from universe (s). Consider some set x that does
                    # not overlap s (i.e., lies entirely ouside of s); the
                    # intersection between x and universe cannot have possibly
                    # changed due to the modification of universe (in which
                    # s is removed from universe) and hence it would be
                    # unnecessary to discard/invalidate x's intersection count.
                    for set_id in list(memoized_intersect_counts[universe_id].keys()):
                        memoized_set = sets[set_id][universe_id]
                        if isinstance(memoized_set, tuple):
                            # memoized_set is a single interval
                            memoized_set_start, memoized_set_end = memoized_set
                        else:
                            memoized_set_start = memoized_set.first_start
                            memoized_set_end = memoized_set.last_end
                        if (memoized_set_start >= s.last_end or
                                memoized_set_end <= s.first_start):
                            # memoized_set lies entirely outside of s, so
                            # there is no need to invalidate its memoized value
                            continue
                        if (isinstance(memoized_set, interval.IntervalSet) and
                            not memoized_set.overlaps_interval(s.first_start,
                                                               s.last_end)):
                            # since memoized_set does not overlap
                            # (s.first_start, s.last_end), it is also true
                            # that memoized_set does not overlap s -- so there
                            # is not need to invalidate its memoized value
                            continue
                        # memoized_set might overlap s, so invalidate it
                        # ('might' because we only checked against
                        #  (s.first_start, s.last_end))
                        # (we could fully check with a call to
                        #  memoized_set.intersection(s) and seeing if the
                        #  result is nonempty, but this is just as expensive
                        #  as the consequence of invalidating memoized_set)
                        del memoized_intersect_counts[universe_id][set_id]
                else:
                    # The universe was modified. Since we are not using
                    # interval sets, there are no obvious optimizations --
                    # simply discard all memoized intersection counts that
                    # involve this universe.
                    memoized_intersect_counts[universe_id] = {}
            else:
                # Although the universe was not modified, the memoized
                # intersection count between id_min_ratio and universe_id
                # is no longer needed (and will never be accessed) because
                # id_min_ratio will never be considered again, as it was
                # placed in the set cover. Because the universe did not
                # change, the count is still valid and therefore there is
                # no strict need to invalidate (i.e., discard) the memoized
                # value. However, doing so can yield some modest runtime
                # improvements when using interval sets because the loop
                # in the optimization above (over
                # memoized_intersect_counts[universe_id].keys()) would no
                # longer needlessly iterate over id_min_ratio.
                if id_min_ratio in memoized_intersect_counts[universe_id]:
                    del memoized_intersect_counts[universe_id][id_min_ratio]

    return set_ids_in_cover
