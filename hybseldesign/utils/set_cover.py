"""Functions for working with instances of the set cover problem.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'


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
    num_left_to_cover = len(universe) - num_that_can_be_uncovered

  return set_ids_in_cover

