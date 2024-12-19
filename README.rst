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

=======
Usage
=======

::

  usage: mosaic-queen [-h] -i INPUT -t TARGET_IMAGES [TARGET_IMAGES ...] -o
                    OUTPUT -n NAME [-j NUM_WORKERS] [-a] [-r] [-f]
                    [-bc BEAM_CUTOFF] [-mc MOSAIC_CUTOFF] [-u]
                    [-s {mad,rms,fit}] [-g GUESS_STD] [-ra RA] [-dec DEC]
                    [-v VELOCITY] [-dra DRA] [-ddec DDEC] [-dv DV]

  Run make_mosaic over the targets

  optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        The directory that contains the (2D or 3D) images and
                        beams.
  -t TARGET_IMAGES [TARGET_IMAGES ...], --target-images TARGET_IMAGES [TARGET_IMAGES ...]
                        The filenames of each target/pointing image to be
                        mosaicked. A suffix of 'image.fits' is expected, and
                        this is replaced by 'pb.fits' in order to locate the
                        corresponding beams (which are also required as
                        input).
  -o OUTPUT, --output OUTPUT
                        The directory for all output files.
  -n NAME, --name NAME  The prefix to be used for output files.
  -j NUM_WORKERS, --num-workers NUM_WORKERS
                        Number of worker threads. Default=0 means all
                        available threads.
  -a, --associated-mosaics
                        Also make mosaics of the associated 'model' and
                        'residual' .fits files.
  -r, --regrid          Use montage for regridding the images and beams.Also
                        regrid the 'model' and 'residual' files, if '--
                        associated-mosaics' is enabled.
  -f, --force-regrid    If the user wants newly-regridded files, this '--
                        force-regrid' argument should be enabled.(If '--
                        regrid' is enabled instead, the package will first
                        check whether regridded files already exist.If they
                        are found, regridding will not proceed because this is
                        a time-consuming step.)
  -bc BEAM_CUTOFF, --beam-cutoff BEAM_CUTOFF
                        The cutoff in the primary beam to use.E.g. The default
                        of 0.1 means going down to the 10 percent level for
                        each pointing. Set to zero for no primary beam cutoff.
  -mc MOSAIC_CUTOFF, --mosaic-cutoff MOSAIC_CUTOFF
                        Sensitivity cutoff in the final mosaic. Pixels with a
                        noise level > minimum mosaic noise / cutoff are
                        blanked in all final products. E.g. The default of 0.2
                        means blanking in the mosaic all pixels with a noise
                        level > 5x the minimum mosaic noise level. Set to zero
                        for no cutoff (but some cutoff may still result from
                        the -bc setting).
  -u, --unity-weights   Build the mosaic using weight=1 instead of
                        weight=1/noise**2 for the input images.
  -s {mad,rms,fit}, --statistic {mad,rms,fit}
                        State 'mad' (median absolute deviation), 'rms' (root
                        mean square) or 'fit' (Gaussian fit) as the statistic
                        to be used for estimating the noise level in the input
                        images. This will be derived using the negative pixel-
                        values. The noise levels set the weights=1/noise**2
                        used when mosaicking. Not used if the '-u' option is
                        enabled. Default is mad.
  -g GUESS_STD, --guess-std GUESS_STD
                        An initial guess of the noise level in the input
                        images, if user has set '--statistic' to 'fit'.(This
                        is to aid a Gaussian fit to the negative pixel-
                        values.) The default of 0.02 assumes that the pixel
                        values are in units of Jy/beam, so a std of ~ 20
                        mJy/beam).
  -ra RA                Central RA (in degrees) of the output mosaic
                        image/cube, if the user does not want to image the
                        entire FoV covered by the input images/cubes.
  -dec DEC              Central Dec (in degrees) of the output mosaic
                        image/cube, if the user does not want to image the
                        entire FoV covered by the input images/cubes.
  -v VELOCITY, --velocity VELOCITY
                        Central velocity/frequency of the output mosaic cube
                        (in the appropriate units of the input cubes) if the
                        user does not want to image the entire
                        velocity/frequency range covered by the input cubes.
  -dra DRA              RA range of the output mosaic image/cube (in degrees),
                        if the user does not want to image the entire FoV
                        covered by the input images/cubes.
  -ddec DDEC            Dec range of the output mosaic image/cube (in
                        degrees), if the user does not want to image the
                        entire FoV covered by the input images/cubes.
  -dv DV                Velocity/frequency range of the output mosaic cube (in
                        the unit used by the input images), if the user does
                        not want to image the entire velocity/frequency range
                        covered by the input cubes.

