"""Tests for dynamic_load module.
"""

import os
import unittest

from catch.utils import dynamic_load

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class TestLoadModule(unittest.TestCase):
    """Tests loading a module.
    """

    def setUp(self):
        # Get the path to a test module relative to the path of this file
        self.module_path = os.path.join(os.path.dirname(__file__),
            'input/test_module.py')

    def test_load_module(self):
        # Load the module and check that it has the variable 'testvar'
        # equal to 'testval'
        module = dynamic_load.load_module_from_path(self.module_path)
        self.assertEqual(module.testvar, 'testval')


class TestLoadFunction(unittest.TestCase):
    """Tests loading a function from a module.
    """

    def setUp(self):
        # Get the path to a test module relative to the path of this file
        self.module_path = os.path.join(os.path.dirname(__file__),
            'input/test_module.py')

    def test_no_such_function(self):
        # Test loading a function that does not exist; should raise an
        # exception
        self.assertRaises(Exception, dynamic_load.load_function_from_path,
                self.module_path, 'no_such_fn')

    def test_load_function(self):
        # Load the module and execute the function 'sum', checking its
        # output
        sum_fn = dynamic_load.load_function_from_path(self.module_path, 'sum')

        self.assertEqual(sum_fn(1, 2), 3)
