#!/usr/bin/env python

import os
from setuptools import setup, find_packages

requirements = [
"numpy",
"astropy",
"memory-profiler",
"matplotlib",
"scipy"
]

build_root = os.path.dirname(__file__)

def readme():
    """Get readme content for package long description"""
    with open(os.path.join(build_root, 'README.rst')) as f:
        return f.read()

PACKAGE_NAME = 'mosaic-queen'
__version__ = '1.1.1'

setup(name = PACKAGE_NAME,
    version = __version__,
    description = "A package with mosaicking commands from montage",
    long_description = readme(),
    long_description_content_type="text/x-rst",
    author = "Sarah White and CaraCal pipeline tean",
    author_email = "sarahwhite.astro@gmail.com",
    url = "https://github.com/caracal-pipeline/MosaicQueen",
    packages = find_packages(),
    install_requires = requirements,
    include_package_data = True,
    scripts = ["bin/" + i for i in os.listdir("bin")],
    classifiers = [])
