"""Tests for set_cover module.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

import unittest

from hybseldesign.utils import set_cover as sc


"""Tests the set cover approximation implementation.
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

  def test_partial_weighted6(self):
    input = { 0: set([1,2]),
              1: set([2,3,4,5]),
              2: set([3]),
              3: set([4]),
              4: set([5]) }
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

