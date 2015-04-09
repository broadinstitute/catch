"""Setup script for the hybseldesign distribution using Distutils.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

from setuptools import setup
from setuptools import find_packages

import hybseldesign

setup(name='hybseldesign',
      version=hybseldesign.__version__,
      description='Tools to design probes for use in hybrid selection',
      author='Hayden Metsky',
      author_email='hayden@mit.edu',
      packages=find_packages(),
      requires=['numpy'],
      )
