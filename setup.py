# -*- coding: utf-8 -*-
# setup.py
from setuptools import setup, find_packages

setup(
      name="sfftk-plus",
      version="0.1.0.dev0",
      packages=find_packages(),
      author="Paul K. Korir, PhD",
      author_email="pkorir@ebi.ac.uk, paul.korir@gmail.com",
      description="Backend toolkit for working with EMDB-SFF with respect to an OMERO server",
      license="Apache License",
      keywords="EMDB-SFF, SFF, segmentation",
      install_requires=['sfftk', 'psycopg2==2.6', ],
      entry_points={
          'console_scripts': [
              'sffp = sfftkplus.sffplus:main',
              ]
          },
)
