#!/usr/bin/env python

# ------------------------------------------------------------------------------------------------------
# Author of package: Sarah White (sarahwhite.astro@gmail.com) and Sphe Makhathini (sphemakh@gmail.com)
# Based on a mosaicking script by Paolo Serra (paolo80serra@gmail.com)
# ------------------------------------------------------------------------------------------------------

from montage_mosaic import make_mosaic
from argparse import ArgumentParser
import montage_mosaic
import os
import sys

log = montage_mosaic.log


def main(argv):

    parser = ArgumentParser(description="Run make_mosaic over the targets")

    parser.add_argument("-i", "--input",
                        help="The directory that contains the (2D or 3D) images and beams.")
    parser.add_argument("-m", "--mosaic-type",
                        help="State 'continuum' or 'spectral' as the type of mosaic to be made.")
    parser.add_argument("-d", "--domontage", action="store_true",
                        help="Use montage for regridding the images and beams.")
    parser.add_argument("-c", "--cutoff", type=float, default=0.1,
                         help="The cutoff in the primary beam to use (assuming a Gaussian at the moment)."
                              "E.g. The default of 0.1 means going down to the 10 percent level for each pointing.")
    parser.add_argument("-n", "--name", default="mymosaic",
                        help="The prefix to be used for output files.")
    parser.add_argument("-t", "--target-images", action="append",
                         help="The filenames of each target/pointing image to be mosaicked. A suffix of 'image.fits' is expected, and this is replaced by 'pb.fits' in order to locate the corresponding beams (which are also required as input).")
    parser.add_argument("-o", "--output",
                         help="The directory for all output files.")

    args = parser.parse_args(argv)
    input_dir = args.input
    mosaic_type = args.mosaic_type
    cutoff = args.cutoff
    outname = args.name
    output_dir = args.output

    if args.target_images: 
        log.info('Target images = {}'.format(" ".join(args.target_images)))
        images = args.target_images
    else:
        log.error(
            "Must specify the (2D or 3D) images to be mosaicked, each prefixed by '-t '.")
        sys.exit()

    # Throw an error if the user provides only one image
    if len(images) < 2:
        log.error('At least two images must be specified for mosaicking')
        sys.exit()

    beams = [tt.replace('image.fits', 'pb.fits') for tt in images]
    imagesR = [tt.replace('image.fits', 'imageR.fits') for tt in images]
    beamsR = [tt.replace('image.fits', 'pbR.fits') for tt in images]

    for tt in images:
        if not os.path.exists(input_dir+'/'+tt):
            log.error('File {0:s} does not exist'.format(input_dir+'/'+tt))
            sys.exit()

    for bb in beams:
        if not os.path.exists(input_dir+'/'+bb):
            log.error('File {0:s} does not exist'.format(input_dir+'/'+bb))
            sys.exit()

    log.info('All images and beams found on disc')


    if args.domontage:
        make_mosaic.use_montage_for_regridding(
            input_dir, output_dir, mosaic_type, images, beams, imagesR, beamsR, outname)
    else:
        log.info(
            'Will use mosaic header {0:s}.hdr and regridded images and beams available on disc'.format(outname))

    make_mosaic.check_for_regridded_files(output_dir, imagesR, beamsR)

    make_mosaic.make_mosaic_using_beam_info(input_dir, output_dir, mosaic_type, outname, imagesR, beamsR, cutoff, images)

    # Move the log file to the output directory
    os.system('mv log-make_mosaic.txt '+output_dir+'/')

    return 0
