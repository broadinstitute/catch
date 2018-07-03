"""Functions for dynamically loading modules and functions.
"""

import importlib
import os

__author__ = 'Hayden Metsky <hayden@mit.edu>'


def load_module_from_path(path):
    """Load Python module in the given path.

    Args:
        path: path to .py file

    Returns:
        Python module (before returning, this also executes
        the module)
    """
    path = os.path.abspath(path)

    # Use the filename (without extension) as the module name
    _, filename = os.path.split(path)
    module_name, _ = os.path.splitext(filename)

    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)

    # Execute the module
    spec.loader.exec_module(module)

    return module


def load_function_from_path(path, fn_name):
    """Load Python function in a module at the given path.

    Args:
        path: path to .py file
        fn_name: name of function in the module

    Returns:
        Python function

    Raises:
        Exception if the module at path does not contain a function
        with name fn_name
    """
    module = load_module_from_path(path)

    if not hasattr(module, fn_name):
        raise Exception(("Module at %s does not contain function %s" %
            (path, fn_name)))

    return getattr(module, fn_name)
