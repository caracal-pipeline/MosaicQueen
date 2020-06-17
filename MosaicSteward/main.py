#!/usr/bin/env python

# ------------------------------------------------------------------------------------------------------
# Author of package: Sarah White (sarahwhite.astro@gmail.com) and Sphe Makhathini (sphemakh@gmail.com)
# Based on a mosaicking script by Paolo Serra (paolo80serra@gmail.com) 
# and a convolving script by Landman Bester (landman.bester@gmail.com)
# ------------------------------------------------------------------------------------------------------

from MosaicSteward import make_mosaic, image_convolver
from argparse import ArgumentParser
import MosaicSteward
import os
import sys
import multiprocessing

log = MosaicSteward.log

# So that error handling is compatible with Python 2 as well as Python 3
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

def main(argv):

    parser = ArgumentParser(description="Run make_mosaic over the targets")

    parser.add_argument("-i", "--input",
                        help="The directory that contains the (2D or 3D) images "
                             "and beams.")
    parser.add_argument("-m", "--mosaic-type",
                        help="State 'continuum' or 'spectral' as the type of "
                             "mosaic to be made.")
    parser.add_argument("-d", "--domontage", action="store_true",
                        help="Use montage for regridding the images and beams.")
    parser.add_argument("-u", "--uniform-resolution", action="store_true",
                        help="Convolve all images and beams to a uniform "
                             "resolution before mosaicking.")
    parser.add_argument("-pm", "--psf-mode", default="auto",
                        help="State 'auto' or 'override' for determining the "
                             "uniform resolution (psf) to be used (if enabled). \n"
                             "The default is 'auto', meaning that the largest psf "
                             "across the input images will be found and used.")
    parser.add_argument("-pp", "--psf-pars", default=None, nargs='+', type=float,
                        help="Beam parameters, specified as emaj emin pa, for the psf "
                             "to be used for enforcing uniform resolution, if wishing to " 
                             "override the 'auto' setting.")
    parser.add_argument("-cp", "--circ-psf", action="store_true",
                        help="Passing this flag will convolve with a circularised "
                             "beam instead of an elliptical one.")
    parser.add_argument("-ncpu", '--ncpu', default=0, type=int,
                        help="Number of threads to use for convolution. \n"
                             "Default of zero means use all threads")
    parser.add_argument("-c", "--cutoff", type=float, default=0.1,
                        help="The cutoff in the primary beam to use (assuming a "
                             "Gaussian at the moment). \n"
                             "E.g. The default of 0.1 means going down to the "
                             "10 percent level for each pointing.")
    parser.add_argument("-n", "--name", default="mymosaic",
                        help="The prefix to be used for output files.")
    parser.add_argument("-t", "--target-images", action="append",
                        help="The filenames of each target/pointing image to be "
                             "mosaicked. A suffix of 'image.fits' is expected, and "
                             "this is replaced by 'pb.fits' in order to locate the "
                             "corresponding beams (which are also required as input).")
    parser.add_argument("-o", "--output",
                        help="The directory for all output files.")

    args = parser.parse_args(argv)
    input_dir = args.input
    mosaic_type = args.mosaic_type
    psf_mode = args.psf_mode
    cutoff = args.cutoff
    outname = args.name
    output_dir = args.output

    if args.target_images: 
        log.info('Target images = {}'.format(" ".join(args.target_images)))
        images = args.target_images
    else:
        log.error(
            "Must specify the (2D or 3D) images to be mosaicked, each prefixed by '-t '.")
        raise LookupError("Must specify the (2D or 3D) images to be mosaicked, each prefixed by '-t '.")

    # Throw an error if the user provides only one image
    if len(images) < 2:
        log.error('At least two images must be specified for mosaicking')
        raise ValueError('At least two images must be specified for mosaicking')

    beams = [tt.replace('image.fits', 'pb.fits') for tt in images]
    imagesR = [tt.replace('image.fits', 'imageR.fits') for tt in images]
    beamsR = [tt.replace('image.fits', 'pbR.fits') for tt in images]

    for tt in images:
        try:
            open(input_dir+'/'+tt)
        except FileNotFoundError:
            log.error('File {0:s} does not exist'.format(input_dir+'/'+tt))
            raise FileNotFoundError('File {0:s} does not exist'.format(input_dir+'/'+tt))

    for bb in beams:
        try:
            open(input_dir+'/'+bb)
        except FileNotFoundError:
            log.error('File {0:s} does not exist'.format(input_dir+'/'+bb)) 
            raise FileNotFoundError('File {0:s} does not exist'.format(input_dir+'/'+bb))

    log.info('All images and beams found on disc')

    # Stage where uniform-resolution is 'applied' (if enabled)

    # Multiprocessing to speed up convolution
    if not args.ncpu:
        args.ncpu = multiprocessing.cpu_count()
    log.info("Using {0:i} threads".format(args.ncpu))

    if args.uniform_resolution:
     
        if psf_mode = 'auto':

            psf_to_use = make_mosaic.find_largest_BMAJ(input_dir, images, mosaic_type, 'images')
            psf_to_use_arcsec = psf_to_use*3600.0   # Since BMAJ is (or should be) in units of deg 
            log.info(
                    "With psf-mode set to 'auto', the input images will be convolved so that they "
                    "have a uniform resolution of {0:f} arcsec".format(psf_to_use_arcsec))

            ### To simplify things for the moment, have the 'auto' setting being to convolve with  cicularised beam
            beampars = tuple([psf_to_use, psf_to_use, 0.0])  ### CHECK I'VE SET THIS UP RIGHT

        else:
            
            if args.psf_pars is None:
                log.error("If wishing to override the 'auto' setting, "
                          "the user must specify the psf parameters to be used for convolution.")
                raise TypeError("If wishing to override the 'auto' setting, "
                          "the user must specify the psf parameters to be used for convolution.")
            else:
                beampars = tuple(args.psf_pars)

            psf_to_use_arcsec = beampars[0]  # User is asked to pass this value in units of arcsec
            psf_to_use = psf_to_use_arcsec/3600.0  # Need to pass this to convolve_image in units of deg
            log.info(
                    "With psf-mode set to 'override', the input images will be convolved so that "
                    "they have a uniform resolution of {0:f} arcsec".format(psf_to_use_arcsec))

            if args.circ_psf:  ### 'auto' setting is a circularised beam, so no warning needed
                log.info("WARNING: Enabling circularised beam. User must set circ-psf to 'False' if they "
                         "want their 'emin' and 'pa' values set through psf-pars to be used.")
                beampars[1] = beampars[0]  # If BPA is varying a lot over the input images, then best to set emin to emaj
                beampars[2] = 0.0  # pa of psf set to zero

        log.info("Psf paramters to be used: emaj = {0:.3f}, emin = {0:.3f}, PA = {0:.3f}".format(beampars[0], beampars[1], beampars[2])) ### CHECK FORMATTING

        for image in images:
            image_convolver.convolve_image(input_dir, image, beampars, args.ncpu)
    
        ### Need to ensure that the primary beams and *convolved* images are going to be passed to make_mosaic 

    else:

        log.info(
                "Will use the 'native' synthesised beams of the input images, with no convolution "
                "to a single resolution before mosaicking. If uniform resolution across the input "
                "images is desired, before mosaicking, please enable 'uniform-resolution' and re-run "
                "this worker (with consideration of the related settings).")
    
        
    # Ready for re-gridding
    if args.domontage:
        make_mosaic.use_montage_for_regridding(
            input_dir, output_dir, mosaic_type, images, beams, imagesR, beamsR, outname)
    else:
        log.info(
                "Will use mosaic header {0:s}.hdr and regridded images and beams available on disc. "
                "WARNING: We assume that the user is happy with the resolution used for these existing, "
                "regridded images. If not, please re-run this worker after enabling 'uniform-resolution' "
                "and 'domontage' (in order to redo the regridding).".format(outname))

    make_mosaic.check_for_regridded_files(output_dir, imagesR, beamsR)


    # Now to mosaic
    make_mosaic.make_mosaic_using_beam_info(input_dir, output_dir, mosaic_type, outname, imagesR, beamsR, cutoff, images)

    # Move the log file to the output directory
    os.system('mv log-make_mosaic.txt '+output_dir+'/')

    return 0
