#!/usr/bin/env python

import os
from setuptools import setup

requirements = [
"numpy",
"futures",
"astropy",
"montage-wrapper"
]

PACKAGE_NAME = 'MosaicSteward'
__version__ = '0.0.3'

setup(name = PACKAGE_NAME,
    version = __version__,
    description = "A package with mosaicking commands from montage",
    author = "Sarah White",
    author_email = "sarahwhite.astro@gmail.com",
    url = "https://github.com/svw26",
    packages=["MosaicSteward"],
    install_requires = requirements,
    include_package_data = True,
    license=["GNU GPL v2"],
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
