# -*- coding: utf-8 -*-

from distutils.core import setup

# Dynamically calculate the version based on src.VERSION.
version = __import__('src').get_version()

setup(name='stepic_bioinformatics',
      version=version,
      package_dir={'': 'src'},
      packages=['stepic_bioinformatics'],
      install_requires=['biopython']
      )

