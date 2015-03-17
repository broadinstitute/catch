"""Tests for set_cover module.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

import unittest
import logging
import numpy as np
from collections import defaultdict

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

  def setUp(self):
    # Disable logging
    logging.disable(logging.INFO)

  def test_one_universe_complete_unweighted(self):
    sets = { 0: { 0: set([1,2]) },
             1: { 0: set([1,2,4]) },
             2: { 0: set([2,4]) },
             3: { 0: set([4,5]) },
             4: { 0: set([3]) } }
    desired_output = set([1,3,4])
    self.assertEqual(sc.approx_multiuniverse(sets),
        desired_output)

  def test_two_universes_complete_unweighted(self):
    sets = { 0: { 0: set([1,2]) },
             1: { 0: set([1,2,4]) },
             2: { 0: set([2,4]) },
             3: { 0: set([4]), 1: set([5]) },
             4: { 1: set([3]) } }
    desired_output = set([1,3,4])
    self.assertEqual(sc.approx_multiuniverse(sets),
        desired_output)

  def test_one_universe_partial_unweighted(self):
    sets = { 0: { 0: set([1,2]) },
             1: { 0: set([1,2,4]) },
             2: { 0: set([2,4]) },
             3: { 0: set([4,5]) },
             4: { 0: set([3]) } }
    universe_p = { 0: 0.6 }
    desired_output = set([1])
    self.assertEqual(sc.approx_multiuniverse(sets,
        universe_p=universe_p), desired_output)

  def test_two_universes_partial_unweighted1(self):
    sets = { 0: { 1: set([1,2]) },
             1: { 1: set([1,2,4]) },
             2: { 1: set([2,4]) },
             3: { 0: set([5]), 1: set([4]) },
             4: { 0: set([3]) } }
    universe_p = { 0: 1.0, 1: 0.3 }
    desired_output = set([3,4])
    self.assertEqual(sc.approx_multiuniverse(sets,
        universe_p=universe_p), desired_output)

  def test_two_universes_partial_unweighted2(self):
    sets = { 0: { 0: set([2]), 1: set([1]) },
             1: { 0: set([2,4]), 1: set([1]) },
             2: { 0: set([2,4]) },
             3: { 0: set([4]), 1: set([5]) },
             4: { 0: set([3]) } }
    universe_p = { 0: 1.0, 1: 0.5 }
    desired_output = set([1,4])
    self.assertEqual(sc.approx_multiuniverse(sets,
        universe_p=universe_p), desired_output)

  def test_two_universes_partial_weighted1(self):
    sets = { 0: { 0: set([1,2]) },
             1: { 0: set([1,2]), 1: set([3,4,5]) },
             2: { 1: set([4]) },
             3: { 1: set([5]) },
             4: { 1: set([3]) } }

    costs = { 0: 2, 1: 1000, 2: 3, 3: 1, 4: 10 }
    universe_p = { 0: 0.1, 1: 0.1 }
    desired_output = set([0,3])
    self.assertEqual(sc.approx_multiuniverse(sets,
        costs, universe_p), desired_output)

    costs = { 0: 2, 1: 1000, 2: 3, 3: 1, 4: 10 }
    universe_p = { 0: 0.0, 1: 0.1 }
    desired_output = set([3])
    self.assertEqual(sc.approx_multiuniverse(sets,
        costs, universe_p), desired_output)

    costs = { 0: 2, 1: 1000, 2: 3, 3: 1, 4: 10 }
    universe_p = { 0: 0.5, 1: 0.5 }
    desired_output = set([0,2,3])
    self.assertEqual(sc.approx_multiuniverse(sets,
        costs, universe_p), desired_output)

  def test_two_universes_partial_weighted2(self):
    sets = { 0: { 0: set([1,2]) },
             1: { 0: set([2,3]), 1: set([4,5]) },
             2: { 0: set([3]) },
             3: { 1: set([4]) },
             4: { 1: set([5]) } }
    costs = { 0: 3, 1: 4, 2: 1, 3: 1, 4: 2 }
    universe_p = { 0: 1.0, 1: 0.5 }
    desired_output = set([0,2,3])
    self.assertEqual(sc.approx_multiuniverse(sets,
        costs, universe_p), desired_output)

  def test_two_universes_partial_weighted3(self):
    sets = { 0: { 1: set([1,2]) },
             1: { 0: set([3,4,5]), 1: set([2]) },
             2: { 0: set([3]) },
             3: { 0: set([4]) },
             4: { 0: set([5]) } }

    costs = { 0: 1000, 1: 4, 2: 1, 3: 1, 4: 2 }
    universe_p = { 0: 0.6, 1: 0.5 }
    desired_output = set([1,2,3])
    # The optimal solution is [1] but the approximation fails to
    # find it
    self.assertEqual(sc.approx_multiuniverse(sets,
        costs, universe_p), desired_output)

    costs = { 0: 1000, 1: 4, 2: 1.5, 3: 1.5, 4: 2 }
    universe_p = { 0: 0.6, 1: 0.5 }
    desired_output = set([1])
    self.assertEqual(sc.approx_multiuniverse(sets,
        costs, universe_p), desired_output)

  def test_three_universes_partial_weighted(self):
    sets = { 0: { 0: set([1,2]) },
             1: { 0: set([2]), 1: set([3,4]) },
             2: { 1: set([3]) },
             3: { 1: set([4]), 2: set([6]) },
             4: { 2: set([5]) } }
    costs = { 0: 3, 1: 4, 2: 1, 3: 1, 4: 1000 }
    universe_p = { 0: 0.5, 1: 0.5, 2: 1.0 }
    desired_output = set([0,3,4])
    self.assertEqual(sc.approx_multiuniverse(sets,
        costs, universe_p), desired_output)

  def test_same_value_different_universe1(self):
    sets = { 0: { 0: set([1,2]) },
             1: { 1: set([1]) } }
    universe_p = { 0: 1.0, 1: 1.0 }
    desired_output = set([0,1])
    self.assertEqual(sc.approx_multiuniverse(sets,
        universe_p=universe_p), desired_output)

  def test_same_value_different_universe2(self):
    sets = { 0: { 0: set([1,2]), 1: set([1]) },
             1: { 1: set([1]) } }
    universe_p = { 0: 1.0, 1: 1.0 }
    desired_output = set([0])
    self.assertEqual(sc.approx_multiuniverse(sets,
        universe_p=universe_p), desired_output)

  def test_same_value_different_universe3(self):
    sets = { 0: { 0: set([1,2]), 1: set([2]) },
             1: { 0: set([1,2,3]) } }
    universe_p = { 0: 1.0, 1: 1.0 }
    desired_output = set([0,1])
    self.assertEqual(sc.approx_multiuniverse(sets,
        universe_p=universe_p), desired_output)

  def test_tuple_universe_id(self):
    sets = { 0: { (0,0): set([1,2]), (1,0): set([2]) },
             1: { (0,0): set([1,2,3]) } }
    universe_p = { (0,0): 1.0, (1,0): 1.0 }
    desired_output = set([0,1])
    self.assertEqual(sc.approx_multiuniverse(sets,
        universe_p=universe_p), desired_output)

  """Verifies that a computed set cover achieves the desired
  coverage of each universe.
  """
  def verify_partial_cover(self, sets, universe_p,
      output):
    universes = defaultdict(set)
    for sets_by_universe in sets.values():
      for universe_id, s in sets_by_universe.iteritems():
        universes[universe_id].update(s)
    for universe_id, universe in universes.iteritems():
      covered_elements = set()
      for set_id in output:
        if universe_id in sets[set_id]:
          covered_elements.update(sets[set_id][universe_id])
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

  def test_random(self):
    output_set = self.run_random(False)
    output_array = self.run_random(True)
    self.assertEqual(output_set, output_array)

  """Generates random instances of set cover, computes the solution,
  and verifies that the solution achieves the desired coverage. Also
  verifies that, on average, the solution achieves a reasonable
  reduction in the sum of weights of chosen sets versus choosing
  all sets.
  """
  def run_random(self, use_arrays):
    np.random.seed(1)
    weight_fracs = []
    outputs = []
    for n in xrange(20):
      # Generate the universes
      num_universes = np.random.randint(1, 10)
      universes = {}
      for universe_id in xrange(num_universes):
        universe_size = np.random.randint(100, 500)
        els = set(np.random.randint(0, 5000, size=universe_size))
        universes[universe_id] = els
      # Generate the sets
      num_sets = np.random.randint(500, 1000)
      sets = defaultdict(dict)
      sets_union = defaultdict(set)
      for set_id in xrange(num_sets):
        for universe_id in xrange(num_universes):
          set_size_from_universe = np.random.randint(0, 25)
          if set_size_from_universe > 0:
            els = set(np.random.choice(list(universes[universe_id]),
                      size=set_size_from_universe, replace=False))
            sets[set_id][universe_id] = els
            sets_union[universe_id].update(els)
      # Remove from all universes any elements that don't show
      # up in a set in order to ensure that we correctly verify
      # partial coverage
      for universe_id, universe in universes.iteritems():
        universe.intersection_update(sets_union[universe_id])
      # Generate random set costs and random coverage fractions
      costs = { set_id: 1.0+10.0*np.random.random() \
                for set_id in xrange(num_sets) }
      universe_p = { universe_id: np.random.random() \
                      for universe_id in xrange(num_universes) }
      # Compute the set cover
      output = sc.approx_multiuniverse(sets, costs, universe_p,
          use_arrays)
      self.verify_partial_cover(sets, universe_p, output)
      weight_fracs += [self.weight_frac(costs, output)]
      outputs += [output]
    # There's no guarantee that the average weight_frac should be
    # small, but in the average case it should be so test it anyway
    # (e.g., test that it's less than 0.01)
    self.assertLess(np.median(weight_fracs), 0.01)

  def tearDown(self):
    # Re-enable logging
    logging.disable(logging.NOTSET)

