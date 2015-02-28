"""Tests for set_cover module.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

import unittest

from hybseldesign.utils import set_cover as sc


"""Tests the set cover approximation implementation.
"""
class TestSetCoverApprox(unittest.TestCase):

  def test_basic(self):
    input = { 0: set([1,2]),
              1: set([1,2,4]),
              2: set([2,4]),
              3: set([4,5]),
              4: set([3]) }
    desired_output = set([1,3,4])
    self.assertEqual(sc.approx(input), desired_output)

