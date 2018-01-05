"""Tests for param_search module.
"""

import logging
import os
import unittest

import numpy as np

from hybseldesign.pool import param_search
from hybseldesign.utils import pool_probes_io

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class TestStandardSearch(unittest.TestCase):
    """Tests searching for optimal mismatches and cover_extension parameters.
    """

    def setUp(self):
        # Disable logging (including calls of WARNING severity)
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

    def test_round(self):
        self.assertEqual(param_search._round_up(16.2, 5), 20)
        self.assertEqual(param_search._round_up(19.2, 5), 20)
        self.assertEqual(param_search._round_down(16.2, 5), 15)
        self.assertEqual(param_search._round_down(19.2, 5), 15)

    def test_interp_probe_count(self):
        pc = {'d1': {(1, 0): 6000, (2, 10): 5000, (2, 20): 4500,
                (4, 10): 4000, (4, 20): 3000, (4, 30): 2500}}
        interp_fn = param_search._make_interp_probe_count_for_dataset_fn(pc)

        # Test known exact count
        self.assertEqual(interp_fn('d1', 4, 20), 3000)

        # Test interpolation
        num_probes = interp_fn('d1', 3, 15)
        self.assertGreater(num_probes, 3000)
        self.assertLess(num_probes, 5000)

    def test_round_params(self):
        pc = {'d1': {(1, 0): 6000, (1, 10): 5500, (1, 20): 5400, (2, 10): 5000,
                (2, 20): 4500, (4, 10): 4000, (4, 20): 3000, (4, 30): 2500},
              'd2': {(4, 10): 1000, (3, 20): 900, (4, 20): 10}}
        params = [2.5, 12, 4, 15]
        rounded = param_search._round_params(params, pc, 4560)

        self.assertEqual(rounded, [2, 20, 4, 20])

    def test_standard_search_vwafr_typical_counts(self):
        """Integration test with the V-WAfr probe set data."""
        self.assertEqual(self.param_names_vwafr, ('mismatches', 'cover_extension'))
        for max_total_count in [90000, 200000]:
            ss = param_search.standard_search(self.probe_counts_vwafr,
                max_total_count)
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

    def tearDown(self):
        # Re-enable logging
        logging.disable(logging.NOTSET)
