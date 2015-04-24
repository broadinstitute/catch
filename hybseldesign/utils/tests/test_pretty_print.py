"""Tests for pretty_print module.
"""

import unittest

from hybseldesign.utils import pretty_print as pp

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class TestTable(unittest.TestCase):
    """Tests the table function.
    """

    def test_basic(self):
        data = [["abcdef", "ghijklmno", "pqr"],
                ["row1", "1.1", "2"],
                ["row23", "3", "xyzxyz"]]
        col_justify = ["left", "center", "right"]
        expected = ("abcdef ghijklmno    pqr\n"
                    "------ --------- ------\n"
                    "row1      1.1         2\n"
                    "row23      3     xyzxyz\n")
        s = pp.table(data, col_justify, header_underline=True)
        self.assertEqual(s, expected)

    def test_with_heights(self):
        data = [["abcdef", "ghijklmno\n123", "pqr"],
                ["row1", "1.1", "2"],
                ["row23", "3", "xyzxyz\n456"]]
        col_justify = ["left", "center", "right"]
        expected = ("abcdef ghijklmno    pqr\n"
                    "          123          \n"
                    "------ --------- ------\n"
                    "row1      1.1         2\n"
                    "row23      3     xyzxyz\n"
                    "                    456\n")
        s = pp.table(data, col_justify, header_underline=True)
        self.assertEqual(s, expected)
