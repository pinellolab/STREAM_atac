#!/usr/bin/env python

from setuptools import setup
import sys

if sys.version_info.major != 3:
    raise RuntimeError('STREAM-atac requires Python 3')


setup(name='stream_atac',
      version="0.3.5",
      description='STREAM preprocessing for single cell atac-seq data',
      url='https://github.com/pinellolab/stream-atac',
      author='Huidong Chen',
      author_email='huidong.chen AT mgh DOT harvard DOT edu',
      license='Affero',
      packages=['stream_atac'],
      package_dir={'stream_atac': 'stream_atac'},
      package_data={'stream_atac': ['run_preprocess.R']},
      include_package_data=True,
      install_requires=[''],
      entry_points = {'console_scripts': ['stream_atac=stream_atac.command_line:main']})

