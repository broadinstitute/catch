"""Setup script for the catch distribution using Distutils.
"""

from setuptools import find_packages
from setuptools import setup

import catch

__author__ = 'Hayden Metsky <hayden@mit.edu>'

setup(name='catch',
      version=catch.__version__,
      description='Tools to design compact, comprehensive probe sets',
      author='Hayden Metsky',
      author_email='hayden@mit.edu',
      packages=find_packages(),
      install_requires=['numpy>=1.9.0', 'scipy>=1.0.0'],
      scripts=[
          'bin/analyze_probe_coverage.py',
          'bin/make_probes.py',
          'bin/make_probes_naively.py',
          'bin/pool_probes.py',
      ])
