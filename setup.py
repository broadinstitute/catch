"""Setup script for the catch distribution using setuptools.
"""

from setuptools import find_packages
from setuptools import setup

import catch

__author__ = 'Hayden Metsky <hayden@mit.edu>'

setup(name='catch',
      version=catch.__version__,
      description='A package for designing compact and comprehensive probe sets.',
      author='Hayden Metsky',
      author_email='hayden@mit.edu',
      packages=find_packages(),
      install_requires=['numpy>=1.15.2', 'scipy>=1.2.0'],
      scripts=[
          'bin/analyze_probe_coverage.py',
          'bin/design.py',
          'bin/design_naively.py',
          'bin/pool.py',
      ])
