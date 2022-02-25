"""Tests for param_search module.
"""

import logging
import os
import unittest

import numpy as np

from catch.pool import param_search
from catch.utils import pool_probes_io

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
              'd2': {(2, 10): 10000, (3, 0): 2000, (3, 10): 1100,
                (4, 10): 1000, (2, 20): 9000, (3, 20): 900, (4, 20): 10}}
        loss_coeffs = (1.0, 1.0/100.0)
        weights = {'d1': 1.0, 'd2': 1.0}

        # Test rounding cover_extension to 10
        params = [2.5, 12, 4, 15]
        rounded = param_search._round_params(params, pc, 4560,
            loss_coeffs, weights, mismatches_round=1, cover_extension_round=10)
        self.assertEqual(rounded, [2, 20, 4, 20])

        # Test rounding cover_extension to 1
        params = [2.5, 12.3, 4, 14.2]
        rounded = param_search._round_params(params, pc, 5500,
            loss_coeffs, weights, mismatches_round=1, cover_extension_round=1)
        for i in range(len(rounded)):
            self.assertEqual(rounded[i], int(rounded[i]))

    def tearDown(self):
        # Re-enable logging
        logging.disable(logging.NOTSET)


class TestSearchFunctions(unittest.TestCase):
    """Tests search functions for searching for optimal parameters.
    """

    def setUp(self):
        # Disable logging (up to the WARNING level)
        logging.disable(logging.WARNING)

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
            # ebola_zaire-with-2014, should not be too high; check that
            # both of them are low
            ebov_mismatches, ebov_cover_extension = opt_params['ebola_zaire-with-2014']
            self.assertLessEqual(ebov_mismatches, 3)
            self.assertLessEqual(ebov_cover_extension, 20)

            # Parameter values for a relatively diverse dataset, like
            # hiv1_without_ltr, should not be too low; check that
            # at least one of them is high
            hiv1_mismatches, hiv1_cover_extension = opt_params['hiv1_without_ltr']
            self.assertTrue(hiv1_mismatches > 3 or hiv1_cover_extension > 20)

    def _search_vwafr_too_small_counts(self, search_fn):
        self.assertEqual(self.param_names_vwafr, ('mismatches', 'cover_extension'))
        for max_total_count in [1, 1000, 10000]:
            with self.assertRaises(param_search.CannotSatisfyProbeCountConstraintError):
                search_fn(max_total_count)

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

    def test_standard_search_vwafr_too_small_counts(self):
        """Integration test with the V-WAfr probe set data."""
        def search_fn(max_total_count):
            return param_search.standard_search(self.probe_counts_vwafr,
                max_total_count)
        self._search_vwafr_too_small_counts(search_fn)

    def test_higher_dimensional_search_vwafr_typical_counts(self):
        """Integration test with the V-WAfr probe set data."""
        def search_fn(max_total_count):
            return param_search.higher_dimensional_search(
                self.param_names_vwafr, self.probe_counts_vwafr,
                max_total_count, loss_coeffs=(1.0, 1.0/100.0))
        self._search_vwafr_typical_counts(search_fn)

    def test_higher_dimensional_search_vwafr_too_small_counts(self):
        """Integration test with the V-WAfr probe set data."""
        def search_fn(max_total_count):
            return param_search.higher_dimensional_search(
                self.param_names_vwafr, self.probe_counts_vwafr,
                max_total_count, loss_coeffs=(1.0, 1.0/100.0))
        self._search_vwafr_too_small_counts(search_fn)

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

        ss = param_search.higher_dimensional_search(param_names, pc, 150000,
            loss_coeffs=(1.0, 1.0, 1.0))
        opt_params, opt_params_count, opt_params_loss = ss

        self.assertLess(opt_params_count, 150000)

        # Verify that the optimal value of the third parameter is
        # around 20 (namely, 10-30)
        for dataset, param_vals in opt_params.items():
            mismatches, cover_extension, p3 = param_vals
            self.assertTrue(10 <= p3 <= 30)

    def test_standard_search_vwafr_with_coefficients(self):
        """Integration test with the V-WAfr probe set data.

        This sets the coefficients in the loss function such that
        mismatches has little impact on loss and cover_extension
        dominates. This should drive the mismatches parameter high
        because that will be yield with smaller probe counts.
        """
        # Note that, by default, loss_coeffs is (1.0, 0.01)
        loss_coeffs = (0.01, 1.0)

        # Test the standard search
        ss = param_search.standard_search(self.probe_counts_vwafr,
            50000, loss_coeffs=loss_coeffs)
        opt_params, opt_params_count, opt_params_loss = ss
        self.assertLess(opt_params_count, 50000)
        for dataset, param_vals in opt_params.items():
            mismatches, cover_extension = param_vals
            # Mismatches should be high
            self.assertGreater(mismatches, 5)

    def test_standard_search_vwafr_with_dataset_weights(self):
        """Integration test with the V-Wafr probe set data.

        This sets the dataset weights in the loss function such that
        all but two datasets have little impact on the loss. The
        other two datasets should dominate the loss and they should
        have small parameter values.
        """
        # Set two relatively diverse datasets (which would generally
        # have high parameter values) to have relatively high weights
        dataset_weights = {d: 1.0 for d in self.probe_counts_vwafr.keys()}
        dataset_weights['hiv1_without_ltr'] = 1000.0
        dataset_weights['hepatitis_c'] = 1000.0

        # Test the standard search
        ss = param_search.standard_search(self.probe_counts_vwafr,
            420000, dataset_weights=dataset_weights)
        opt_params, opt_params_count, opt_params_loss = ss
        self.assertLess(opt_params_count, 420000)

        # Check that both parameter values are small for these diverse
        # datasets
        for d in ['hiv1_without_ltr', 'hepatitis_c']:
            mismatches, cover_extension = opt_params[d]
            self.assertLessEqual(mismatches, 1)
            self.assertLessEqual(cover_extension, 10)

        # Check that at least one other dataset has larger parameter
        # values for each of the 2 parameters
        mismatches_num_large, cover_extension_num_large = 0, 0
        for d in (dataset_weights.keys() -
                set(['hiv1_without_ltr', 'hepatitis_c'])):
            mismatches, cover_extension = opt_params[d]
            if mismatches > 1:
                mismatches_num_large += 1
            if cover_extension > 10:
                cover_extension_num_large += 1
        self.assertGreater(mismatches_num_large, 0)
        self.assertGreater(cover_extension_num_large, 0)

    def tearDown(self):
        # Re-enable logging
        logging.disable(logging.NOTSET)
