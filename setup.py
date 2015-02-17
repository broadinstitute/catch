"""Setup script for the hybseldesign distribution using Distutils.
"""

# Author: Hayden Metsky <hayden@mit.edu>

from setuptools import setup
from setuptools import find_packages

setup(name='hybseldesign',
      version='1.0',
      description='Tools to design probes for use in hybrid selection',
      author='Hayden Metsky',
      author_email='hayden@mit.edu',
      packages=find_packages(),
      requires=['numpy'],
)
