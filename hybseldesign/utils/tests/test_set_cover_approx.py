"""Tests for set_cover_approx module.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

import unittest

from hybseldesign.utils import set_cover_approx as sca


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
    self.assertEqual(sca.set_cover_approx(input), desired_output)

