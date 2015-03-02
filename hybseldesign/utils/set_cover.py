"""Functions for working with instances of the set cover problem.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

import math


"""Approximates the solution to an instance of the set cover problem
using the well-known greedy algorithm.

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
      | S_i \intersect U |)
    C <- C \union { S_i }
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

The input 'sets' is a dict mapping set identifiers to sets. The
optional input 'costs' is a dict mapping set identifiers to the
costs (or weights) of the set; the default is for every set to have
a cost of 1, making this equivalent to the unweighted problem. The
optional input 'p' is a float in [0,1] that specifies the fraction
of the universe we must cover; the default is p=1, making this
equivalent to the problem in which the entire universe is covered.
The output is a set consisting of the identifiers of the sets chosen
to be in the set cover. For example, if the sets input is
  { 0: S_1, 1: S_2, 2: S_3 }
where each S_i is a set, and the sets S_1 and S_3 are chosen to be
in the set cover, then the output is {0,2}.
"""
def approx(sets, costs=None, p=1.0):
  if p < 0 or p > 1:
    raise ValueError("p must be in [0,1]")
  if costs == None:
    # Give each set a default cost of 1
    costs = { set_id: 1 for set_id in sets.keys() }
  else:
    for c in costs.values():
      if c < 0:
        raise ValueError("All costs must be nonnegative")

  # Create the universe as the union of all input sets
  universe = set()
  for s in sets.values():
    universe.update(s)

  num_that_can_be_uncovered = int(len(universe) - p*len(universe))
  # Above, use int(..) to take the floor. Also, expand out
  # len(universe) rather than use int((1.0-p)*len(universe)) due to
  # precision errors in Python -- e.g., int((1.0-0.8)*5) yields 0 on
  # some machines.
  num_left_to_cover = len(universe) - num_that_can_be_uncovered

  set_ids_not_in_cover = sets.keys()
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
    num_left_to_cover = max(0,
        len(universe) - num_that_can_be_uncovered)

  return set_ids_in_cover


"""Approximates the solution to an instance of a version of the
set cover problem in which there are multiple universes and we
seek to find a collection of sets whose union covers a specified
fraction of each given universe.

The input 'sets' is a dict mapping set identifiers to sets. The
input 'universes' is a dict mapping universe identifiers to sets
consisting of the elements in the corresponding universe. The
optional input 'costs' is a dict mapping set identifiers to the
costs (or weights) of the set; the default is for every set to have
a cost of 1, making this equivalent to the unweighted problem. The
optional input 'universe_p' is a dict mapping universe identifiers
to floats in [0,1] that specifiy the fraction of the corresponding
universe we must cover; the default is for the coverage fraction of
each universe to be 1.0, making this equivalent to the problem in
which there is a single universe (the union of all the given
universes) that must be completely covered.  The output is a set
consisting of the identifiers of the sets chosen to be in the set
cover. For example, if the sets input is
  { 0: S_1, 1: S_2, 2: S_3 }
where each S_i is a set, and the sets S_1 and S_3 are chosen to be
in the set cover, then the output is {0,2}.

This is a generalization of the partial, weighted set cover problem
whose solution is approximated by approx(..). (For any input to
that function, simply create a single universe whose elements
consist of the union of all the input sets, and that supply that
universe as input to this function.) Indeed, this function and
approx(..) are very similar, with this one simply generalizing to
allow for multiple universes. However, whereas approx(..) achieves
a provable guarantee on the approximation factor, it is not clear
whether the straightforward generalization implemented in this
function achieves the same factor; therefore, there may be room for
improvement in this function. Furthermore, the generalization to
multiple universes will yield a slight increase in runtime (constant
factors) compared to working directly with the more traditional
notion of a single universe. For these reasons -- although the
approx(..) function could be greatly shortened by simply calling
this function -- we choose to implement the versions separately.
"""
def approx_multiuniverse(sets, universes, costs=None,
    universe_p=None):
  if costs == None:
    # Give each set a default cost of 1
    costs = { set_id: 1 for set_id in sets.keys() }
  else:
    for c in costs.values():
      if c < 0:
        raise ValueError("All costs must be nonnegative")
  if universe_p == None:
    # Give each universe a coverage fraction of 1.0 (i.e., cover
    # all of it)
    universe_p = { universe_id: 1 for universe_id in universes.keys() }
  else:
    for p in universe_p.values():
      if p < 0 or p > 1:
        raise ValueError(("The coverage fraction (p) of each "
                          "universe must be in [0,1]"))
    for universe_id in universes.keys():
      if universe_id not in universe_p:
        raise ValueError(("universe_p is missing a value for "
                          "universe %d" % universe_id))

  # Check that it is possible to achieve the desired coverage of
  # each universe (i.e., by taking all the elements in the union of
  # all sets)
  sets_union = set()
  for s in sets.values():
    sets_union.update(s)
  for universe_id, universe in universes.iteritems():
    p = universe_p[universe_id]
    num_possible = len(universe.intersection(sets_union))
    num_desired = math.ceil(p*len(universe))
    if num_possible < num_desired:
      raise ValueError(("%d elements of universe %d must be "
                        "covered, but it is only possible to cover "
                        "%d elements" % (num_desired, universe_id,
                        num_possible)))

  # Make new sets of all the universes so they can be modified
  # without disrupting the state of this function's caller
  universes = { k: set(universes[k]) for k in universes.keys() }

  num_that_can_be_uncovered = {}
  num_left_to_cover = {}
  for universe_id, universe in universes.iteritems():
    p = universe_p[universe_id]
    num_that_can_be_uncovered[universe_id] = \
        int(len(universe) - p*len(universe))
    # Above, use int(..) to take the floor. Also, expand out
    # len(universe) rather than use int((1.0-p)*len(universe)) due to
    # precision errors in Python -- e.g., int((1.0-0.8)*5) yields 0 on
    # some machines.
    num_left_to_cover[universe_id] = \
        len(universe) - num_that_can_be_uncovered[universe_id]

  set_ids_not_in_cover = sets.keys()
  set_ids_in_cover = set()
  # Keep iterating until desired partial cover of each universe
  # is obtained (note that [] evaluates to False)
  while [True for universe_id in universes.keys() \
        if num_left_to_cover[universe_id] > 0]:
    # Find the set that minimizes the ratio of its cost to the
    # number of uncovered elements (that need to be covered) that
    # it covers
    id_min_ratio, min_ratio = None, float('inf')
    for id in set_ids_not_in_cover:
      # There is no strict need to keep track of set_ids_not_in_cover
      # and iterate through these; we could have iterated over
      # all input sets. However, sets that are already put into the
      # set cover have zero intersection with any universe (i.e., cover
      # none of it), so there is no reason to iterate over the sets
      # already placed in the cover. Not doing so should yield some
      # runtime improvements.
      s = sets[id]
      num_needed_covered_across_universes = 0
      for universe_id, universe in universes.iteritems():
        num_covered = len(s.intersection(universe))
        # There is no need to cover more than num_left_to_cover
        # elements for this universe
        num_needed_covered = min(num_left_to_cover[universe_id],
                                  num_covered)
        num_needed_covered_across_universes += num_needed_covered
      if num_needed_covered_across_universes == 0:
        # s covers no elements that need to be covered, so it should
        # not be in the set cover
        continue
      ratio = float(costs[id]) / num_needed_covered_across_universes
      if ratio < min_ratio:
        id_min_ratio = id
        min_ratio = ratio
    # id_min_ratio goes into the set cover
    set_ids_in_cover.add(id_min_ratio)
    set_ids_not_in_cover.remove(id_min_ratio)
    for universe_id, universe in universes.iteritems():
      universe.difference_update(sets[id_min_ratio])
      num_left_to_cover[universe_id] = max(0,
          len(universe) - num_that_can_be_uncovered[universe_id])

  return set_ids_in_cover
