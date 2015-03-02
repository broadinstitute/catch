"""Tests for set_cover module.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

import unittest
import numpy as np

from hybseldesign.utils import set_cover as sc


"""Tests approx function.
"""
class TestSetCoverApprox(unittest.TestCase):

  def test_complete_unweighted(self):
    input = { 0: set([1,2]),
              1: set([1,2,4]),
              2: set([2,4]),
              3: set([4,5]),
              4: set([3]) }
    desired_output = set([1,3,4])
    self.assertEqual(sc.approx(input), desired_output)

  def test_partial_unweighted1(self):
    input = { 0: set([1,2]),
              1: set([1,2,4]),
              2: set([2,4]),
              3: set([4,5]),
              4: set([3]) }
    desired_output = set([1])
    self.assertEqual(sc.approx(input, p=0.6), desired_output)

  def test_partial_unweighted2(self):
    input = { 0: set([1,2]),
              1: set([1,2,4]),
              2: set([2,4]),
              3: set([4,5]),
              4: set([2,3,6]) }
    desired_output = set([1,4])
    self.assertEqual(sc.approx(input, p=0.81), desired_output)

  def test_complete_weighted1(self):
    input = { 0: set([1,2]),
              1: set([1,2,4]),
              2: set([2,4]),
              3: set([4,5]),
              4: set([3]) }
    costs = { 0: 2, 1: 1000, 2: 3, 3: 1, 4: 10 }
    desired_output = set([0,3,4])
    self.assertEqual(sc.approx(input, costs=costs), desired_output)

  def test_complete_weighted2(self):
    input = { 0: set([1,2]),
              1: set([1,2,3,4,5]),
              2: set([4]),
              3: set([5]),
              4: set([3]) }
    costs = { 0: 2, 1: 1000, 2: 3, 3: 1, 4: 10 }
    desired_output = set([0,2,3,4])
    self.assertEqual(sc.approx(input, costs=costs), desired_output)

  def test_partial_weighted1(self):
    input = { 0: set([1,2]),
              1: set([1,2,3,4,5]),
              2: set([4]),
              3: set([5]),
              4: set([3]) }
    costs = { 0: 2, 1: 1000, 2: 3, 3: 1, 4: 10 }
    desired_output = set([3])
    self.assertEqual(sc.approx(input, costs=costs, p=0.1),
        desired_output)

  def test_partial_weighted2(self):
    input = { 0: set([1,2]),
              1: set([2,3]),
              2: set([4,5]),
              3: set([5]),
              4: set([4]) }
    costs = { 0: 2, 1: 1000, 2: 100, 3: 10, 4: 10 }
    desired_output = set([0,3,4])
    self.assertEqual(sc.approx(input, costs=costs, p=0.7),
        desired_output)

  def test_partial_weighted3(self):
    input = { 0: set([1,2]),
              1: set([3]),
              2: set([4]),
              3: set([2,5]),
              4: set([1]) }
    costs = { 0: 2, 1: 1000, 2: 999, 3: 10, 4: 10 }
    desired_output = set([0,2,3])
    self.assertEqual(sc.approx(input, costs=costs, p=0.8),
        desired_output)

  def test_partial_weighted4(self):
    input = { 0: set([1,2]),
              1: set([3,4,5]),
              2: set([3]),
              3: set([4]),
              4: set([5]) }
    costs = { 0: 2.1, 1: 3, 2: 2, 3: 2, 4: 2 }
    desired_output = set([1])
    self.assertEqual(sc.approx(input, costs=costs, p=0.6),
        desired_output)

  def test_partial_weighted5(self):
    input = { 0: set([1,2]),
              1: set([2,3,4,5]),
              2: set([3]),
              3: set([4]),
              4: set([5]) }

    costs = { 0: 3, 1: 4, 2: 1, 3: 1, 4: 2 }
    desired_output = set([1])
    self.assertEqual(sc.approx(input, costs=costs, p=0.8),
        desired_output)

    costs = { 0: 3, 1: 4.1, 2: 1, 3: 1, 4: 2 }
    desired_output = set([0,2,3])
    # The optimal solution is [1], but the approximation fails to
    # find it
    self.assertEqual(sc.approx(input, costs=costs, p=0.8),
        desired_output)

  def test_no_elements(self):
    input = { }
    self.assertEqual(sc.approx(input), set([]))
    inptut = { 0: set([]) }
    self.assertEqual(sc.approx(input), set([]))

  def test_one_element(self):
    input = { 0: set([1]) }
    self.assertEqual(sc.approx(input), set([0]))


"""Tests approx_multiuniverse function.
"""
class TestSetCoverApproxMultiuniverse(unittest.TestCase):

  def test_one_universe_complete_unweighted(self):
    sets = { 0: set([1,2]),
             1: set([1,2,4]),
             2: set([2,4]),
             3: set([4,5]),
             4: set([3]) }
    universes = { 0: set([1,2,3,4,5]) }
    desired_output = set([1,3,4])
    self.assertEqual(sc.approx_multiuniverse(sets, universes),
        desired_output)

  def test_two_universes_complete_unweighted(self):
    sets = { 0: set([1,2]),
             1: set([1,2,4]),
             2: set([2,4]),
             3: set([4,5]),
             4: set([3]) }

    universes = { 0: set([1,2,4]), 1: set([3,5]) }
    desired_output = set([1,3,4])
    self.assertEqual(sc.approx_multiuniverse(sets, universes),
        desired_output)

    universes = { 0: set([1,2,4]) }
    desired_output = set([1])
    self.assertEqual(sc.approx_multiuniverse(sets, universes),
        desired_output)

  def test_impossible_to_cover(self):
    sets = { 0: set([1,2]),
             1: set([1,2,4]),
             2: set([2,4]),
             3: set([4,5]),
             4: set([3]) }
    universes = { 0: set([1,2,6]), 1: set([3,5]) }
    self.assertRaises(ValueError,
        sc.approx_multiuniverse, sets, universes)

  def test_one_universe_partial_unweighted(self):
    sets = { 0: set([1,2]),
             1: set([1,2,4]),
             2: set([2,4]),
             3: set([4,5]),
             4: set([3]) }
    universes = { 0: set([1,2,3,4,5]) }
    universe_p = { 0: 0.6 }
    desired_output = set([1])
    self.assertEqual(sc.approx_multiuniverse(sets, universes,
        universe_p=universe_p), desired_output)

  def test_two_universes_partial_unweighted(self):
    sets = { 0: set([1,2]),
             1: set([1,2,4]),
             2: set([2,4]),
             3: set([4,5]),
             4: set([3]) }

    universes = { 0: set([3,5]), 1: set([1,2,4]) }
    universe_p = { 0: 1.0, 1: 0.3 }
    desired_output = set([3,4])
    self.assertEqual(sc.approx_multiuniverse(sets, universes,
        universe_p=universe_p), desired_output)

    universes = { 0: set([2,3,4]), 1: set([1,5]) }
    universe_p = { 0: 1.0, 1: 0.5 }
    desired_output = set([1,4])
    self.assertEqual(sc.approx_multiuniverse(sets, universes,
        universe_p=universe_p), desired_output)

  def test_two_universes_partial_weighted1(self):
    sets = { 0: set([1,2]),
             1: set([1,2,3,4,5]),
             2: set([4]),
             3: set([5]),
             4: set([3]) }

    universes = { 0: set([1,2]), 1: set([3,4,5]) }
    costs = { 0: 2, 1: 1000, 2: 3, 3: 1, 4: 10 }
    universe_p = { 0: 0.1, 1: 0.1 }
    desired_output = set([0,3])
    self.assertEqual(sc.approx_multiuniverse(sets, universes,
        costs, universe_p), desired_output)

    universes = { 0: set([1,2]), 1: set([3,4,5]) }
    costs = { 0: 2, 1: 1000, 2: 3, 3: 1, 4: 10 }
    universe_p = { 0: 0.0, 1: 0.1 }
    desired_output = set([3])
    self.assertEqual(sc.approx_multiuniverse(sets, universes,
        costs, universe_p), desired_output)

    universes = { 0: set([1,2]), 1: set([3,4,5]) }
    costs = { 0: 2, 1: 1000, 2: 3, 3: 1, 4: 10 }
    universe_p = { 0: 0.5, 1: 0.5 }
    desired_output = set([0,2,3])
    self.assertEqual(sc.approx_multiuniverse(sets, universes,
        costs, universe_p), desired_output)

  def test_two_universes_partial_weighted2(self):
    sets = { 0: set([1,2]),
             1: set([2,3,4,5]),
             2: set([3]),
             3: set([4]),
             4: set([5]) }

    universes = { 0: set([1,2,3]), 1: set([4,5]) }
    costs = { 0: 3, 1: 4, 2: 1, 3: 1, 4: 2 }
    universe_p = { 0: 1.0, 1: 0.5 }
    desired_output = set([0,2,3])
    self.assertEqual(sc.approx_multiuniverse(sets, universes,
        costs, universe_p), desired_output)

    universes = { 0: set([3,4,5]), 1: set([1,2]) }
    costs = { 0: 1000, 1: 4, 2: 1, 3: 1, 4: 2 }
    universe_p = { 0: 0.6, 1: 0.5 }
    desired_output = set([1,2,3])
    # The optimal solution is [1] but the approximation fails to
    # find it
    self.assertEqual(sc.approx_multiuniverse(sets, universes,
        costs, universe_p), desired_output)

    universes = { 0: set([3,4,5]), 1: set([1,2]) }
    costs = { 0: 1000, 1: 4, 2: 1.5, 3: 1.5, 4: 2 }
    universe_p = { 0: 0.6, 1: 0.5 }
    desired_output = set([1])
    self.assertEqual(sc.approx_multiuniverse(sets, universes,
        costs, universe_p), desired_output)

  def test_three_universes_partial_weighted(self):
    sets = { 0: set([1,2]),
             1: set([2,3,4]),
             2: set([3]),
             3: set([4,6]),
             4: set([5]) }
    universes = { 0: set([1,2]), 1: set([3,4]), 2: set([5,6]) }
    costs = { 0: 3, 1: 4, 2: 1, 3: 1, 4: 1000 }
    universe_p = { 0: 0.5, 1: 0.5, 2: 1.0 }
    desired_output = set([0,3,4])
    self.assertEqual(sc.approx_multiuniverse(sets, universes,
        costs, universe_p), desired_output)

  def test_empty_universe(self):
    sets = { 0: set([1,2]) }
    universes = { 0: set([1]), 1: set() }
    universe_p = { 0: 1.0, 1: 1.0 }
    self.assertEqual(sc.approx_multiuniverse(sets, universes,
        universe_p=universe_p), set([0]))

  """Verifies that a computed set cover achieves the desired
  coverage of each universe.
  """
  def verify_partial_cover(self, sets, universes, universe_p,
      output):
    covered_elements = set()
    for set_id in output:
      covered_elements.update(sets[set_id])
    for universe_id, universe in universes.iteritems():
      desired_p = universe_p[universe_id]
      desired_num_covered = desired_p * len(universe)
      num_covered = len(covered_elements.intersection(universe))
      self.assertGreaterEqual(num_covered, desired_num_covered)

  """Returns the ratio of the sum of the weights of the sets in a
  computed set cover to the sum of the weights of all the sets.
  """
  def weight_frac(self, costs, output):
    sum_of_all_weights = np.sum(costs.values())
    sum_of_weights_of_selected_sets = \
        np.sum([costs[set_id] for set_id in output])
    return float(sum_of_weights_of_selected_sets) / sum_of_all_weights

  """Generates random instances of set cover, computes the solution,
  and verifies that the solution achieves the desired coverage. Also
  verifies that, on average, the solution achieves a reasonable
  reduction in the sum of weights of chosen sets versus choosing
  all sets.
  """
  def test_random(self):
    np.random.seed(1)
    weight_fracs = []
    for n in xrange(25):
      # Generate the universes
      num_universes = np.random.randint(1, 10)
      universes = {}
      universe_p = {}
      universe_union = set()
      max_of_previous_universe = 0
      for universe_id in xrange(num_universes):
        universe_size = np.random.randint(100, 500)
        els = set(np.random.randint(max_of_previous_universe+1,
                  max_of_previous_universe+500, size=universe_size))
        universes[universe_id] = els
        max_of_previous_universe = max(els)
        universe_union.update(els)
      universe_union_list = list(universe_union)
      # Generate the sets
      num_sets = np.random.randint(500, 1000)
      sets = {}
      sets_union = set()
      for set_id in xrange(num_sets):
        set_size = np.random.randint(10, 250)
        els = set(np.random.choice(universe_union_list,
                  size=set_size, replace=False))
        sets[set_id] = els
        sets_union.update(els)
      # Remove from all universes any elements that don't show
      # up in a set in order to ensure that all universes can be
      # covered fully
      for universe in universes.values():
        universe.intersection_update(sets_union)
      # Generate random set costs and random coverage fractions
      costs = { set_id: 1.0+10.0*np.random.random() \
                for set_id in xrange(num_sets) }
      universe_p = { universe_id: np.random.random() \
                      for universe_id in xrange(num_universes) }
      # Compute the set cover
      output = sc.approx_multiuniverse(sets, universes, costs,
                universe_p)
      self.verify_partial_cover(sets, universes, universe_p,
          output)
      weight_fracs += [self.weight_frac(costs, output)]
    # There's no guarantee that the average weight_frac should be
    # small, but in the average case it should be so test it anyway
    # (e.g., test that it's less than 0.01)
    self.assertLess(np.median(weight_fracs), 0.01)
