#!/usr/bin/env python

import os
from setuptools import setup, find_packages

requirements = [
"numpy",
"futures",
"astropy",
"montage-wrapper",
"memory-profiler",
]

PACKAGE_NAME = 'mosaic-queen'
__version__ = '1.0.1'

setup(name = PACKAGE_NAME,
    version = __version__,
    description = "A package with mosaicking commands from montage",
    author = "Sarah White and CaraCal pipeline tean",
    author_email = "sarahwhite.astro@gmail.com",
    url = "https://github.com/caracal-pipeline/MosaicQueen",
    packages = find_packages(),
    install_requires = requirements,
    include_package_data = True,
    license="GNU GPL v2",
    scripts = ["bin/" + i for i in os.listdir("bin")],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: Astronomy"
    ]
     )
