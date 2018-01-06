"""Tests for interpolate_count module.
"""

import logging
import unittest

from hybseldesign.pool import interpolate_count as ic

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class TestInterpolateCount(unittest.TestCase):
    """Tests functions for interpolating probe count.
    """

    def setUp(self):
        # Disable logging
        logging.disable(logging.INFO)

    def test_round(self):
        self.assertEqual(ic._round_up(16.2, 5), 20)
        self.assertEqual(ic._round_up(19.2, 5), 20)
        self.assertEqual(ic._round_down(16.2, 5), 15)
        self.assertEqual(ic._round_down(19.2, 5), 15)

    def test_interp_probe_count_standard(self):
        pc = {'d1': {(1, 0): 6000, (2, 10): 5000, (2, 20): 4500,
                (4, 10): 4000, (4, 20): 3000, (4, 30): 2500}}
        interp_fn = ic._make_interp_probe_count_for_dataset_standard_fn(pc)

        # Test known exact count
        self.assertEqual(interp_fn('d1', [4, 20]), 3000)

        # Test interpolation
        num_probes = interp_fn('d1', [3, 15])
        self.assertGreater(num_probes, 3000)
        self.assertLess(num_probes, 5000)

    def test_interp_probe_count_nd(self):
        pc = {'d1': {(1, 0): 6000, (2, 10): 5000, (2, 20): 4500,
                (4, 10): 4000, (4, 20): 3000, (4, 30): 2500}}
        interp_fn = ic._make_interp_probe_count_for_dataset_nd_fn(pc)

        # Test known exact count
        self.assertEqual(interp_fn('d1', [4, 20]), 3000)

        # Test interpolation
        num_probes = interp_fn('d1', [3, 15])
        self.assertGreater(num_probes, 3000)
        self.assertLess(num_probes, 5000)


    def tearDown(self):
        # Re-enable logging
        logging.disable(logging.NOTSET)
