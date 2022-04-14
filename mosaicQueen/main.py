#!/usr/bin/env python

# ------------------------------------------------------------------------------------------------------
# Authors of package: Sarah White (sarahwhite.astro@gmail.com)
#                     Sphe Makhathini (sphemakh@gmail.com)
#                     Paolo Serra (paolo.serra@inaf.it)
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

def check_for_files(directory, fits_files, type_of_fits_file, regrid_boolean):
    # Remember that if regridded files were produced by MosaicKing, we expect them to be in the output directory
    # Check location for beams made by MosaicKing

    dont_exist = False
    exit_if_not_found_list = [ 'images', 'beams', 'models', 'residuals' ]
    for ff in fits_files:
        try:
            open(directory+'/'+ff)
        except FileNotFoundError:
            dont_exist = True
            if type_of_fits_file in exit_if_not_found_list:
                log.error('File {0:s} does not exist'.format(directory+'/'+ff))
                raise FileNotFoundError
            else:  # Intended for the regridded files
                if regrid_boolean:
                    log.warning('File {0:s} does not exist'.format(directory+'/'+ff))
                else:
                    log.error('File {0:s} does not exist'.format(directory+'/'+ff))
                    raise FileNotFoundError
    if dont_exist == False:  # i.e they do exist(!)
        log.info('All {0:s} found on disk'.format(type_of_fits_file))

    return dont_exist


def main(argv):

    parser = ArgumentParser(description="Run make_mosaic over the targets")

    parser.add_argument("-i", "--input",
                        help="The directory that contains the (2D or 3D) images and beams.")
    parser.add_argument("-t", "--target-images", nargs='+', required = True,
                        help="The filenames of each target/pointing image to be mosaicked. A suffix of 'image.fits' is expected, "
                             "and this is replaced by 'pb.fits' in order to locate the corresponding beams (which are also required as input).")
    parser.add_argument("-o", "--output",
                        help="The directory for all output files.")
    parser.add_argument("-n", "--name", default="mymosaic",
                        help="The prefix to be used for output files.")
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
    parser.add_argument("-bc", "--beam-cutoff", type=float, default=0.1,
                        help="The cutoff in the primary beam to use."
                              "E.g. The default of 0.1 means going down to the 10 percent level for each pointing.")
    parser.add_argument("-mc", "--mosaic-cutoff", type=float, default=0.2,
                        help="Sensitivity cutoff in the final mosaic. Pixels with a noise level > minimum mosaic noise / cutoff are blanked in all final products. "
                              "E.g. The default of 0.2 means blanking in the mosaic all pixels with a noise level > 5x the minimum mosaic noise level.")
    parser.add_argument("-u", "--unity-weights", action="store_true",
                        help="Build the mosaic using weight=1 instead of weight=1/noise**2 for the input images.")
    parser.add_argument("-s", "--statistic", choices= ["mad", "rms", "fit"], required = False, default = "mad",
                        help="State 'mad' (median absolute deviation), 'rms' (root mean square) or 'fit' (Gaussian fit) as the statistic to be "
                             "used for estimating the noise level in the input images. This will be derived using the negative pixel-values. "
                             "The noise levels set the weights=1/noise**2 used when mosaicking. Not used if the '-u' option is enabled. Default is mad.")
    parser.add_argument("-g", "--guess-std", type=float, default=0.02,
                        help="An initial guess of the noise level in the input images, if user has set '--statistic' to 'fit'."
                             "(This is to aid a Gaussian fit to the negative pixel-values.) The default of 0.02 assumes that "
                             "the pixel values are in units of Jy/beam, so a std of ~ 20 mJy/beam).")

    args = parser.parse_args(argv)
    input_dir = args.input
    mosaic_type = args.mosaic_type
    beam_cutoff = args.beam_cutoff
    mosaic_cutoff = args.mosaic_cutoff
    statistic = args.statistic
    unity_weights = args.unity_weights
    sigma_guess = args.guess_std
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
    check_for_files(input_dir, images, 'images', args.regrid)  # This raises an error and exits if files are not found
    check_for_files(input_dir, beams, 'beams', args.regrid)  # This raises an error and exits if files are not found

    if args.associated_mosaics:  # Want to be sure that all of the ingredients are in place before doing any mosaicking
        log.info('Checking for models and residuals')
        check_for_files(input_dir, models, 'models', args.regrid)  # This raises an error and exits if files are not found
        check_for_files(input_dir, residuals, 'residuals', args.regrid)  # This raises an error and exits if files are not found

    # Check here BITPIX of input images
    bitpix = make_mosaic.find_lowest_precision(input_dir, images)

    if args.force_regrid:
        log.info('You have asked for all regridded files to be created by this run, even if they are already on disk')
        make_mosaic.use_montage_for_regridding(
            input_dir, output_dir, mosaic_type, 'image', images, imagesR, beams, beamsR, outname, bitpix)
        make_mosaic.use_montage_for_regridding(
            input_dir, output_dir, mosaic_type, 'pb', images, imagesR, beams, beamsR, outname, bitpix)
    elif args.regrid:
        log.info('Checking for regridded images and beams')
        imagesR_dont_exist = check_for_files(output_dir, imagesR, 'regridded images', args.regrid)
        if imagesR_dont_exist:
            log.info('Regridded images are not all in place, so using montage to create them')
            make_mosaic.use_montage_for_regridding(
                input_dir, output_dir, mosaic_type, 'image', images, imagesR, beams, beamsR, outname, bitpix)
        else:
            log.info('Regridded images are all in place')
        beamsR_dont_exist = check_for_files(output_dir, beamsR, 'regridded beams', args.regrid)
        if beamsR_dont_exist:
            log.info('Regridded beams are not all in place, so using montage to create them')
            make_mosaic.use_montage_for_regridding(
                input_dir, output_dir, mosaic_type, 'pb', images, imagesR, beams, beamsR, outname, bitpix)
        #else:  # redundant
        #    log.info('Regridded beams are all in place')
    else:
        log.info('User specified neither --force-regrid nor --regrid in the mosaic-queen command')
        log.info(
            'Will use mosaic header {0:s}.hdr and regridded images and beams available on disk'.format(outname))
        make_mosaic.final_check_for_files(output_dir, imagesR, beamsR)  # This function raises an error and exits if files are not found


    log.info("Mosaicking 'image' files")
    noises = make_mosaic.make_mosaic_using_beam_info(input_dir, output_dir, mosaic_type, 'image', outname, imagesR, beamsR, beam_cutoff, unity_weights, statistic, sigma_guess, images, mosaic_cutoff, bitpix)


    if args.associated_mosaics:  # Code is more readable by keeping these mosaics separate

        log.info("Mosaicking 'model' and 'residual' files")
        if args.force_regrid:
            log.info('You have asked for all regridded files to be created by this run, even if they are already on disk')
            make_mosaic.use_montage_for_regridding(
                input_dir, output_dir, mosaic_type, 'model', models, modelsR, beams, beamsR, outname, bitpix)
            make_mosaic.use_montage_for_regridding(
                input_dir, output_dir, mosaic_type, 'residual', residuals, residualsR, beams, beamsR, outname, bitpix)
        elif args.regrid:
            log.info('Checking for regridded models and residuals')
            modelsR_dont_exist = check_for_files(output_dir, modelsR, 'regridded models', args.regrid)
            if modelsR_dont_exist:
                log.info('Regridded models are not all in place, so using montage to create them')
                make_mosaic.use_montage_for_regridding(
                    input_dir, output_dir, mosaic_type, 'model', models, modelsR, beams, beamsR, outname, bitpix)
            #else:  # redundant
            #    log.info('Regridded models are all in place')
            residualsR_dont_exist = check_for_files(output_dir, residualsR, 'regridded residuals', args.regrid)
            if residualsR_dont_exist:
                log.info('Regridded residuals are not all in place, so using montage to create them')
                make_mosaic.use_montage_for_regridding(
                    input_dir, output_dir, mosaic_type, 'residual', residuals, residualsR, beams, beamsR, outname, bitpix)
            #else:  # redundant
            #    log.info('Regridded residuals are all in place')
        else:
            log.info('User specified neither --force-regrid nor --regrid in the mosaic-queen command')
            log.info(
                'Will use regridded models and residuals available on disk'.format(outname))
            check_for_files(output_dir, modelsR, 'regridded models', args.regrid)
            check_for_files(output_dir, residualsR, 'regridded residuals', args.regrid)

        make_mosaic.make_mosaic_using_beam_info(input_dir, output_dir, mosaic_type, 'model'   , outname, modelsR   , beamsR, beam_cutoff, unity_weights, statistic, sigma_guess, models   , mosaic_cutoff, bitpix, all_noise_estimates=noises)
        make_mosaic.make_mosaic_using_beam_info(input_dir, output_dir, mosaic_type, 'residual', outname, residualsR, beamsR, beam_cutoff, unity_weights, statistic, sigma_guess, residuals, mosaic_cutoff, bitpix, all_noise_estimates=noises)

    # Move the log file to the output directory
    os.system('mv log-make_mosaic.txt '+output_dir+'/')

    return 0
