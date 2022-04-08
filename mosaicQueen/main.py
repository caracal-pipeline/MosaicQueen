#!/usr/bin/env python

# ------------------------------------------------------------------------------------------------------
# Author of package: Sarah White (sarahwhite.astro@gmail.com) and Sphe Makhathini (sphemakh@gmail.com)
# Based on a mosaicking script by Paolo Serra (paolo80serra@gmail.com)
# ------------------------------------------------------------------------------------------------------

from mosaicQueen import make_mosaic
from argparse import ArgumentParser
import mosaicQueen
import os
import sys
import glob

log = mosaicQueen.log

# So that error handling is compatible with Python 2 as well as Python 3
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

def check_for_files(input_dir, images):

    dont_exist = False 
    for tt in images:
        try:
            open(input_dir+'/'+tt)
        except FileNotFoundError:
            log.error('File {0:s} does not exist'.format(input_dir+'/'+tt))
            dont_exist = True 

    return dont_exist


def main(argv):

    parser = ArgumentParser(description="Run make_mosaic over the targets")

    parser.add_argument("-i", "--input",
                        help="The directory that contains the (2D or 3D) images and beams.")
    parser.add_argument("-m", "--mosaic-type", choices= ["spectral", "continuum"], required = True,
                        help="State 'continuum' or 'spectral' as the type of mosaic to be made.")
    parser.add_argument("-a", "--associated-mosaics", action="store_true",
                        help="Also make mosaics of the associated 'model' and 'residual' .fits files.")
    parser.add_argument("-r", "--regrid", action="store_true",
                        help="Use montage for regridding the images and beams."
                              "Also regrid the 'model' and 'residual' files, if '--associated-mosaics' is enabled.")
    parser.add_argument("-f", "--force-regrid", action="store_true",
                        help="If the user wants newly-regridded files, this '--force-regrid' argument should be enabled."
                              "(If '--regrid' is enabled instead, the package will first check whether regridded files already exist." 
                              "If they are found, regridding will not proceed because this is a time-consuming step.)")
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
    os.makedirs(output_dir, exist_ok=True)

    if args.target_images: 
        if len(args.target_images) == 1:
            images = glob.glob(os.path.join(args.input, args.target_images[0]))
            images = [os.path.basename(item) for item in images]
        else:
            images = args.target_images
        log.info('Target images = {}'.format(" ".join(args.target_images)))
    else:
        log.error(
            "Must specify the (2D or 3D) images to be mosaicked, each prefixed by '-t '.")
        raise LookupError("Must specify the (2D or 3D) images to be mosaicked, each prefixed by '-t '.")

    # # Throw an error if the user provides only one image
    # if len(images) < 2:
    #    log.error('At least two images must be specified for mosaicking')
    #    raise ValueError('At least two images must be specified for mosaicking')

    # 'R' to signify the regridded versions of the different .fits files
    imagesR = [tt.replace('image.fits', 'imageR.fits') for tt in images]
    beams = [tt.replace('image.fits', 'pb.fits') for tt in images]
    beamsR = [tt.replace('image.fits', 'pbR.fits') for tt in images]
    if args.associated_mosaics:
        log.info('Will generate mosaics made from the input images, their associated models, and their residuals')
        models = [tt.replace('image.fits', 'model.fits') for tt in images]
        modelsR = [tt.replace('image.fits', 'modelR.fits') for tt in images]
        residuals = [tt.replace('image.fits', 'residual.fits') for tt in images]
        residualsR = [tt.replace('image.fits', 'residualR.fits') for tt in images]
    else:
        log.info("Will generate a mosaic from the input images. If you would also like mosaics to be made from the associated models and residuals, please re-run with the '--associated-mosaics' argument enabled.")

    log.info('Checking for images and beams')
    make_mosaic.final_check_for_files(input_dir, images, beams)  # This function raises an error and exits if files are not found

    if args.associated_mosaics:  # Want to be sure that all of the ingredients are in place before doing any mosaicking
        log.info('Checking for models and residuals')
        make_mosaic.final_check_for_files(input_dir, models, residuals)  # Function raises an error and should exit if files are not found

    if args.force_regrid:
        log.info('You have asked for all regridded files to be created by this run, even if they are already on disk') 
        make_mosaic.use_montage_for_regridding(
            input_dir, output_dir, mosaic_type, 'image', images, imagesR, beams, beamsR, outname)
        make_mosaic.use_montage_for_regridding(
            input_dir, output_dir, mosaic_type, 'pb', images, imagesR, beams, beamsR, outname)
    elif args.regrid:
        log.info('Checking for regridded images and beams')
        imagesR_dont_exist = check_for_files(output_dir, imagesR)
        if imagesR_dont_exist:
            log.info('Regridded images are not all in place, so using montage to create them')
            make_mosaic.use_montage_for_regridding(
                input_dir, output_dir, mosaic_type, 'image', images, imagesR, beams, beamsR, outname)
        else:
            log.info('Regridded images are all in place')
        beamsR_dont_exist = check_for_files(output_dir, beamsR)
        if beamsR_dont_exist:  
            log.info('Regridded beams are not all in place, so using montage to create them')
            make_mosaic.use_montage_for_regridding(
                input_dir, output_dir, mosaic_type, 'pb', images, imagesR, beams, beamsR, outname)
        else:
            log.info('Regridded beams are all in place')
    else:
        log.info(
            'Will use mosaic header {0:s}.hdr and regridded images and beams available on disk'.format(outname))
        make_mosaic.final_check_for_files(output_dir, imagesR, beamsR)  # This function raises an error and exits if files are not found 

    make_mosaic.make_mosaic_using_beam_info(input_dir, output_dir, mosaic_type, 'image', outname, imagesR, beamsR, cutoff, images)


    if args.associated_mosaics:  # Code is more readable by keeping these mosaics separate

        if args.force_regrid:
            log.info('You have asked for all regridded files to be created by this run, even if they are already on disk')
            make_mosaic.use_montage_for_regridding(
                input_dir, output_dir, mosaic_type, 'model', models, modelsR, beams, beamsR, outname)
            make_mosaic.use_montage_for_regridding(
                input_dir, output_dir, mosaic_type, 'residual', residuals, residualsR, beams, beamsR, outname)
        elif args.regrid:
            log.info('Checking for regridded models and residuals')
            modelsR_dont_exist = check_for_files(output_dir, modelsR)
            if modelsR_dont_exist:
                log.info('Regridded models are not all in place, so using montage to create them')
                make_mosaic.use_montage_for_regridding(
                    input_dir, output_dir, mosaic_type, 'model', models, modelsR, beams, beamsR, outname)
            else:
                log.info('Regridded models are all in place')
            residualsR_dont_exist = check_for_files(output_dir, residualsR)
            if residualsR_dont_exist:
                log.info('Regridded residuals are not all in place, so using montage to create them')
                make_mosaic.use_montage_for_regridding(
                    input_dir, output_dir, mosaic_type, 'residual', residuals, residualsR, beams, beamsR, outname)
            else:
                log.info('Regridded residuals are all in place')
        else:
            log.info(
                'Will use regridded models and residuals available on disk'.format(outname))
            make_mosaic.final_check_for_files(output_dir, modelsR, residualsR)

        make_mosaic.make_mosaic_using_beam_info(input_dir, output_dir, mosaic_type, 'model', outname, modelsR, beamsR, cutoff, models)
        make_mosaic.make_mosaic_using_beam_info(input_dir, output_dir, mosaic_type, 'residual', outname, residualsR, beamsR, cutoff, residuals)

    # Move the log file to the output directory
    os.system('mv log-make_mosaic.txt '+output_dir+'/')

    return 0
