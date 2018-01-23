"""Utilities involving logging.
"""

import logging

__author__ = 'Hayden Metsky <hayden@mit.edu>'


def configure_logging(level=logging.WARNING):
    """Configure the global logging package.

    Args:
        level: logging level below which logging messages are ignored
    """
    format = '%(asctime)s - %(name)s [%(levelname)s] %(message)s'
    logging.basicConfig(format=format, level=level)
