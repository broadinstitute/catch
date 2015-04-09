"""Tests for interval module.
"""

import unittest

from hybseldesign.utils import interval

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class TestMergeOverlapping(unittest.TestCase):

  """Tests the merge_overlapping function.
  """

  def compare(self, input, desired_output):
    self.assertEqual(interval.merge_overlapping(input),
      desired_output)

  def test_single(self):
    self.compare([(1,2)], [(1,2)])

  def test_two_nonoverlapping(self):
    self.compare([(1,2), (4,5)], [(1,2), (4,5)])

  def test_two_overlapping(self):
    self.compare([(1,3), (2,5)], [(1,5)])

  def test_two_touching(self):
    self.compare([(1,3), (3,5)], [(1,5)])

  def test_three(self):
    self.compare([(1,5), (3,7), (9,12)], [(1,7), (9,12)])

  def test_within(self):
    self.compare([(1,5), (3,4), (6,8)], [(1,5), (6,8)])

  def test_loner_in_middle(self):
    self.compare([(1,5), (3,5), (5,8), (10,12), (15,18), (17,20)],
                  [(1,8), (10,12), (15,20)])

class TestSchedule(unittest.TestCase):

  """Tests the schedule function.
  """

  def compare(self, input, desired_output):
    self.assertEqual(interval.schedule(input),
      desired_output)

  def test_single(self):
    self.compare([((1,2), 1)], [1])

  def test_two_nonoverlapping(self):
    self.compare([((1,2), 1), ((4,5), 2)], [1,2])
    self.compare([((4,5), 2), ((1,2), 1)], [1,2])

  def test_two_overlapping(self):
    self.compare([((1,3), 1), ((2,5), 2)], [1])
    self.compare([((2,5), 2), ((1,3), 1)], [1])

  def test_three(self):
    self.compare([((1,5), 1), ((3,7), 2), ((9,12), 3)], [1,3])

  def test_within(self):
    self.compare([((1,5), 1), ((3,4), 2), ((6,8), 3)], [2,3])

  def test_loner_in_middle(self):
    self.compare([((1,5), 1), ((3,4), 2), ((5,8), 3), ((10,12), 4),
                  ((15,18),5), ((17,20), 6)], [2,3,4,5])

