#!/usr/bin/env python

from setuptools import setup
import sys

if sys.version_info.major != 3:
    raise RuntimeError('STREAM-atac requires Python 3')


setup(name='stream-atac',
      version="0.1.0",
      description='Single-cell Trajectories Reconstruction, Exploration And Mapping of single-cell data. Preprocessing for single cell atac-seq data',
      url='https://github.com/pinellolab/stream-atac',
      author='Huidong Chen',
      author_email='huidong.chen AT mgh DOT harvard DOT edu',
      license='Affero',
      packages=['stream-atac'],
      package_dir={'stream-atac': 'stream-atac'},
      install_requires=[''],)

