"""Utilities involving logging.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

import logging


def configure_logging(level=logging.WARNING):
  format = '%(asctime)s - %(name)s [%(levelname)s] %(message)s'
  logging.basicConfig(format=format, level=level)

