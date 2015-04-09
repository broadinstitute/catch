"""Tests for duplicate_filter module.
"""

import unittest

from hybseldesign.filter import duplicate_filter
from hybseldesign import probe

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class TestDuplicateFilter(unittest.TestCase):
    """Tests the duplicate filter output on contrived input.
    """

    def test_basic(self):
        input = ['ATCGTCGCGG', 'ATCGTAGCGG', 'ATCGTCACGG', 'ATCGTAGCGG',
                 'ATTGTCGCGG', 'ATCGTCGCGG']
        desired_output = ['ATCGTCGCGG', 'ATCGTAGCGG', 'ATCGTCACGG',
                          'ATTGTCGCGG']
        input_probes = [probe.Probe.from_str(s) for s in input]
        desired_output_probes = [probe.Probe.from_str(s)
                                 for s in desired_output]
        f = duplicate_filter.DuplicateFilter()
        f.filter(input_probes)
        self.assertItemsEqual(f.input_probes, input_probes)
        # Order should be preserved, so use assertEqual rather than
        # assertItemsEqual
        self.assertEqual(f.output_probes, desired_output_probes)
