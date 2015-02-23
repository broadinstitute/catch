"""Approximates the solution to an instance of the set cover problem
using the well-known greedy algorithm.

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

The input 'sets' is a dict mapping set identifiers to sets.
The output is a set consisting of the identifiers of the sets chosen
to be in the set cover. For example, if the input is
  { 0: S_1, 1: S_2, 2: S_3 }
where each S_i is a set, and the sets S_1 and S_3 are chosen to be
in the set cover, then the output is {0,2}.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'


def set_cover_approx(sets):
  # Create the universe as the union of all input sets
  universe = set()
  for s in sets.values():
    universe.update(s)

  set_ids_not_in_cover = sets.keys()
  set_ids_in_cover = set()
  # Keep iterating until universe is covered
  while len(universe) > 0:
    # Find the set the covers the most of universe
    id_covering_most, most_num_covered = None, 0
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
      if num_covered > most_num_covered:
        id_covering_most = id
        most_num_covered = num_covered
    # id_covering_most goes into the set cover
    set_ids_in_cover.add(id_covering_most)
    set_ids_not_in_cover.remove(id_covering_most)
    universe.difference_update(sets[id_covering_most])

  return set_ids_in_cover

