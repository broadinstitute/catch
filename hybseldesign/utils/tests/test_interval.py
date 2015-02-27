"""Tests for interval module.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

import unittest

from hybseldesign.utils import interval


"""Tests the merge_overlapping function.
"""
class TestMergeOverlapping(unittest.TestCase):

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

