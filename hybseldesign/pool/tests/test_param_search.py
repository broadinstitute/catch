"""Tests for param_search module.
"""

import logging
import os
import unittest

import numpy as np

from hybseldesign.pool import param_search
from hybseldesign.utils import pool_probes_io

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class TestHelperFunctions(unittest.TestCase):
    """Tests helper functions for searching for optimal parameters.
    """

    def setUp(self):
        # Disable logging
        logging.disable(logging.INFO)

        # Set a seed for numpy's random generator so the initial
        # guess is always the same
        np.random.seed(1)

    def test_bounds_and_initial_guess_standard(self):
        pc = {'d1': {(1, 10): 6000, (2, 10): 5000, (2, 20): 4500,
                (4, 10): 4000, (4, 20): 3000, (5, 20): 2500}}
        bounds = param_search._make_param_bounds_standard(pc)

        x0 = param_search._make_initial_guess(pc, bounds, 2)

        # The bounding box is [(2, 10), (2, 20), (4, 10), (4, 20)]
        # so the bounds should ensure the guess is within this
        self.assertGreaterEqual(x0[0], 2)
        self.assertLessEqual(x0[0], 4)
        self.assertGreaterEqual(x0[1], 10)
        self.assertLessEqual(x0[1], 20)

    def test_bounds_and_initial_guess_nd(self):
        pc = {'d1': {(1, 0): 6000, (2, 10): 5000, (2, 20): 4500,
                (4, 10): 4000, (4, 20): 3000, (4, 30): 2500},
              'd2': {(6, 40): 10, (7, 60): 20}}

        x0 = param_search._make_initial_guess(pc, None, 2)

        # Since bounds=None, the guess should be one of the given
        # points
        self.assertIn(tuple(x0[0:2]), pc['d1'].keys())
        self.assertIn(tuple(x0[2:4]), pc['d2'].keys())

    def test_round_params(self):
        pc = {'d1': {(1, 0): 6000, (1, 10): 5500, (1, 20): 5400, (2, 10): 5000,
                (2, 20): 4500, (4, 10): 4000, (4, 20): 3000, (4, 30): 2500},
              'd2': {(2, 10): 10000, (3, 10): 1100, (4, 10): 1000,
                (2, 20): 9000, (3, 20): 900, (4, 20): 10}}

        # Test rounding cover_extension to 10
        params = [2.5, 12, 4, 15]
        rounded = param_search._round_params(params, pc, 4560,
            mismatches_round=1, cover_extension_round=10)
        self.assertEqual(rounded, [2, 20, 4, 20])

        # Test rounding cover_extension to 1
        params = [2.5, 12.3, 4, 14.2]
        rounded = param_search._round_params(params, pc, 5500,
            mismatches_round=1, cover_extension_round=1)
        for i in range(len(rounded)):
            self.assertEqual(rounded[i], int(rounded[i]))

    def tearDown(self):
        # Re-enable logging
        logging.disable(logging.NOTSET)


class TestSearchFunctions(unittest.TestCase):
    """Tests search functions for searching for optimal parameters.
    """

    def setUp(self):
        # Disable logging
        logging.disable(logging.INFO)

        # Set a seed for numpy's random generator so the initial
        # guess is always the same
        np.random.seed(1)

        # Read probe counts from actual data
        fn_vall = os.path.join(os.path.dirname(__file__), 'input',
                               'num-probes.V-All.201606.tsv')
        table_vall = pool_probes_io.read_table_of_probe_counts(fn_vall)
        self.param_names_vall, self.probe_counts_vall = table_vall

        fn_vwafr = os.path.join(os.path.dirname(__file__), 'input',
                               'num-probes.V-WAfr.201506.tsv')
        table_vwafr = pool_probes_io.read_table_of_probe_counts(fn_vwafr)
        self.param_names_vwafr, self.probe_counts_vwafr = table_vwafr

    def _search_vwafr_typical_counts(self, search_fn):
        self.assertEqual(self.param_names_vwafr, ('mismatches', 'cover_extension'))
        for max_total_count in [90000, 200000]:
            ss = search_fn(max_total_count)
            opt_params, opt_params_count, opt_param_loss = ss

            # For both counts, the total should be at most the constraint
            # but not too must less
            self.assertLessEqual(opt_params_count, max_total_count)
            self.assertGreater(opt_params_count, 0.9*max_total_count)

            # Parameter values for a relatively conserved dataset, like
            # ebola_zaire-with-2014, should not be too high
            ebov_mismatches, ebov_cover_extension = opt_params['ebola_zaire-with-2014']
            self.assertLessEqual(ebov_mismatches, 3)
            self.assertLessEqual(ebov_cover_extension, 20)

            # Parameter values for a relatively diverse dataset, like
            # hiv1_without_ltr, should not be too low
            hiv1_mismatches, hiv1_cover_extension = opt_params['hiv1_without_ltr']
            self.assertGreaterEqual(hiv1_mismatches, 2)
            self.assertGreaterEqual(hiv1_cover_extension, 20)

    def test_standard_search_vwafr_typical_counts(self):
        """Integration test with the V-WAfr probe set data."""
        def search_fn(max_total_count):
            return param_search.standard_search(self.probe_counts_vwafr,
                max_total_count)
        self._search_vwafr_typical_counts(search_fn)

    def test_standard_search_vwafr_high_count(self):
        """Integration test with the V-WAfr probe set data."""
        # At 1,000,000 probes, it should be able to use values of 0
        # (most conservative) for both parameters
        ss = param_search.standard_search(self.probe_counts_vwafr,
            1000000)
        opt_params, opt_params_count, opt_params_loss = ss

        self.assertLess(opt_params_count, 1000000)

        for dataset, param_vals in opt_params.items():
            mismatches, cover_extension = param_vals
            self.assertEqual(mismatches, 0)
            self.assertEqual(cover_extension, 0)

    def test_higher_dimensional_search_vwafr_typical_counts(self):
        """Integration test with the V-WAfr probe set data."""
        def search_fn(max_total_count):
            return param_search.higher_dimensional_search(
                self.param_names_vwafr, self.probe_counts_vwafr,
                max_total_count)
        self._search_vwafr_typical_counts(search_fn)

    def test_higher_dimensional_search_vwafr_with_third_param(self):
        """Integration test with the V-WAfr probe set data.

        This adds a third parameter to the data, which dominates the
        probe set count. The counts are set so that they are minimized
        when the third parameter is equal to 20. The loss function
        will still try to push this parameter down, but if it has
        a big enough impact on the barrier when it deviates from 25
        (i.e., by exceeding the max probe count) then the minimizer
        should pick a value of roughly 20.
        """
        param_names = ('mismatches', 'cover_extension', 'p3')

        pc = {}
        for dataset in self.probe_counts_vwafr.keys():
            pc[dataset] = {}
            for param_vals in self.probe_counts_vwafr[dataset]:
                count = self.probe_counts_vwafr[dataset][param_vals]
                for k in [0, 10, 20, 30, 40]:
                    # Make a new count that is small when k=25
                    new_count = count + 100000*(k/20.0 - 1)**2
                    param_vals_with_k = tuple(list(param_vals) + [k])
                    pc[dataset][param_vals_with_k] = new_count

        ss = param_search.higher_dimensional_search(param_names, pc, 150000)
        opt_params, opt_params_count, opt_params_loss = ss
        print(opt_params, opt_params_count)

        self.assertLess(opt_params_count, 150000)

        # Verify that the optimal value of the third parameter is
        # around 20 (namely, 18-22)
        for dataset, param_vals in opt_params.items():
            mismatches, cover_extension, p3 = param_vals
            self.assertTrue(18 <= p3 <= 22)

    def tearDown(self):
        # Re-enable logging
        logging.disable(logging.NOTSET)
