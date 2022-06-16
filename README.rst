=============
MosaicQueen
=============

|Pypi Version|

A package with mosaicking commands from montage, that also does primary-beam correction for both 2D and 3D images. Originally written to be used by stimela (https://github.com/SpheMakh/Stimela) and, ultimately, CARACal (https://github.com/caracal-pipeline/caracal) for automated radio-data reduction, but it also works independently. 

==============
Installation
==============

First install montage via

.. code-block:: bash
  
    sudo apt install montage

Package is available on PyPI (https://pypi.org/project/mosaic-queen) and so it is installable via

.. code-block:: bash
  
    pip install mosaic-queen


Work in progress: convolution of input images to a common synthesised beam before mosaicking, and employing parallel processing to speed this up...

.. |Pypi Version| image:: https://img.shields.io/pypi/v/mosaic-queen.svg
                  :target: https://pypi.org/project/mosaic-queen/
                  :alt:
