"""Setup script for the hybseldesign distribution using Distutils.
"""

from setuptools import find_packages
from setuptools import setup

import hybseldesign

__author__ = 'Hayden Metsky <hayden@mit.edu>'


setup(name='hybseldesign',
      version=hybseldesign.__version__,
      description='Tools to design probes for use in hybrid selection',
      author='Hayden Metsky',
      author_email='hayden@mit.edu',
      packages=find_packages(),
      requires=['numpy'],
      )
