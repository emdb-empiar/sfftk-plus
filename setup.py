# -*- coding: utf-8 -*-
# setup.py
from setuptools import setup, find_packages

from sfftkplus import SFFTKPLUS_VERSION

with open('README.md') as f:
    long_description = f.read()

setup(
    name="sfftk-plus",
    version=SFFTKPLUS_VERSION,
    packages=find_packages(),
    author="Paul K. Korir, PhD",
    author_email="pkorir@ebi.ac.uk, paul.korir@gmail.com",
    description="Extended toolkit for working with EMDB-SFF files",
    long_description=long_description,
    license="Apache License",
    keywords="EMDB-SFF, SFF, segmentation",
    install_requires=['sfftk', 'vtk<9.0', 'psycopg2-binary', ],
    entry_points={
        'console_scripts': [
            'sff = sfftkplus.sffplus:main',
        ]
    },
    package_data={
        'sfftkplus': ['sffp.conf', 'schema/roi.xsd', ],
    }
)
