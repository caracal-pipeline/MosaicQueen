#!/usr/bin/env python

# ------------------------------------------------------------------------------------------------------
# Authors of package: Sarah White (sarahwhite.astro@gmail.com)
#                     Sphe Makhathini (sphemakh@gmail.com)
#                     Paolo Serra (paolo.serra@inaf.it)
# ------------------------------------------------------------------------------------------------------

from mosaicqueen import make_mosaic
from argparse import ArgumentParser
import mosaicqueen
import os
import sys
import glob
import psutil
import tracemalloc

log = mosaicqueen.log

def check_for_files(directory, fits_files, type_of_fits_file, regrid_boolean):
    # Remember that if regridded files were produced by MosaicKing, we expect them to be in the output directory
    # Check location for beams made by MosaicKing

    dont_exist = False
    exit_if_not_found_list = [ 'images', 'beams', 'masks', 'models', 'residuals' ]
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


def main():

    parser = ArgumentParser(description="Run make_mosaic over the targets")

    parser.add_argument("-i", "--input", required = True,
                        help="The directory that contains the (2D or 3D) images and beams.")
    parser.add_argument("-t", "--target-images", nargs='+', required=True,
                        help="Space-separated list of images to be mosaicked. The file names must end with 'image.fits', and this is "
                             "automatically replaced by 'pb.fits' to locate the corresponding beams (which are also required as input).")
    parser.add_argument("-o", "--output", required=True,
                        help="The directory for all output files.")
    parser.add_argument("-n", "--name", required = True,
                        help="The prefix to be used for output files.")
    parser.add_argument("-j", "--num-workers", type=int, default=0,
                        help="Number of worker threads. Default=0 means all available threads.")
    parser.add_argument("-a", "--associated-mosaics", nargs='+', default="", choices=["mask", "model", "residual"],
                        help="Also make mosaics of the associated .fits files, selecting any of 'mask', 'model' and 'residual'. Users "
                             "should give their choice as a space-separated list. The defult is to not make any associated mosaics.")
    parser.add_argument("-r", "--regrid", action="store_true",
                        help="Use montage for regridding the images and beams."
                              "Also regrid the 'mask', 'model' and/or 'residual' files, if '--associated-mosaics' is enabled.")
    parser.add_argument("-f", "--force-regrid", action="store_true",
                        help="If the user wants newly-regridded files, this '--force-regrid' argument should be enabled."
                              "(If '--regrid' is enabled instead, the package will first check whether regridded files already exist."
                              "If they are found, regridding will not proceed because this is a time-consuming step.)")
    parser.add_argument("-bc", "--beam-cutoff", type=float, default=0.1,
                        help="The cutoff in the primary beam to use."
                              "E.g. The default of 0.1 means going down to the 10 percent level for each pointing."
                              " Set to zero for no primary beam cutoff.")
    parser.add_argument("-mc", "--mosaic-cutoff", type=float, default=0.2,
                        help="Sensitivity cutoff in the final mosaic. Pixels with a noise level > minimum mosaic noise / cutoff are blanked in all final products. "
                              "E.g. The default of 0.2 means blanking in the mosaic all pixels with a noise level > 5x the minimum mosaic noise level."
                              " Set to zero for no cutoff (but some cutoff may still result from the -bc setting).")
    parser.add_argument("-u", "--unity-weights", action="store_true",
                        help="Build the mosaic using weight=1 instead of weight=1/noise**2 for the input images.")
    parser.add_argument("-s", "--statistic", choices=["mad", "rms", "fit", "exptime"], required=False, default="mad",
                        help="State 'mad' (median absolute deviation), 'rms' (root mean square), 'fit' (Gaussian fit) or 'exptime' (exposure "
                             "time) as the statistic to be used for estimating the noise level in the input images. If selected, the exposure time is "
                             "taken from the header key 'EXPTIME' and the noise is set to 1./sqrt(exptime). In all other cases the noise is derived using the "
                             "negative pixel-values. The noise levels set the weights=1/noise**2 used when mosaicking. Not used if the '-u' option is enabled. "
                             "The default is 'mad'.")
    parser.add_argument("-g", "--guess-std", type=float, default=0.02,
                        help="An initial guess of the noise level in the input images, if user has set '--statistic' to 'fit'."
                             "(This is to aid a Gaussian fit to the negative pixel-values.) The default of 0.02 assumes that "
                             "the pixel values are in units of Jy/beam, so a std of ~ 20 mJy/beam).")
    parser.add_argument("-ra", type=float,
                        help="Central RA (in degrees) of the output mosaic image/cube, if the user does not want to image "
                        "the entire FoV covered by the input images/cubes.")
    parser.add_argument("-dec", type=float,
                        help="Central Dec (in degrees) of the output mosaic image/cube, if the user does not want to image "
                        "the entire FoV covered by the input images/cubes.")
    parser.add_argument("-v", "--velocity", type=float,
                        help="Central velocity/frequency of the output mosaic cube (in the appropriate units of the input cubes) "
                        "if the user does not want to image the entire velocity/frequency range covered by the input cubes.")
    parser.add_argument("-dra", type=float,
                        help="RA range of the output mosaic image/cube (in degrees), if the user does not want to image "
                        "the entire FoV covered by the input images/cubes.")
    parser.add_argument("-ddec", type=float,
                        help="Dec range of the output mosaic image/cube (in degrees), if the user does not want to image "
                        "the entire FoV covered by the input images/cubes.")
    parser.add_argument("-dv", type=float,
                        help="Velocity/frequency range of the output mosaic cube (in the unit used by the input images), if the user does not want to image "
                        "the entire velocity/frequency range covered by the input cubes.")

    args = parser.parse_args(sys.argv[1:])
    input_dir = args.input
    num_workers = args.num_workers
    beam_cutoff = args.beam_cutoff
    mosaic_cutoff = args.mosaic_cutoff
    statistic = args.statistic
    unity_weights = args.unity_weights
    sigma_guess = args.guess_std
    outname = args.name
    output_dir = args.output

    os.makedirs(output_dir, exist_ok=True)
    if not num_workers:
        num_workers = psutil.cpu_count()

    subimage_dict = {
        'CRVAL1': args.ra,
        'CRVAL2': args.dec,
        'CRVAL3': args.velocity,
        'dRA':    args.dra,
        'dDec':   args.ddec,
        'dv':     args.dv,
    }

    if args.target_images:
        if len(args.target_images) == 1:
            images = glob.glob(os.path.join(args.input, args.target_images[0]))
            images = [os.path.basename(item) for item in images]
        else:
            images = args.target_images
        log.info('Target images:')
        for im in images:
            log.info('    {}'.format(im))
    else:
        log.error(
            "Must specify the (2D or 3D) images to be mosaicked.")
        raise LookupError("Must specify the (2D or 3D) images to be mosaicked.")

    # check NAXIS and NAXISi of input FITS files to decide whether to follow the 2D or 3D mosaic workflow
    naxis, nlongaxis = make_mosaic.find_naxis(input_dir, images)
    if nlongaxis == 2:
        mosaic_type = 'continuum'
    elif nlongaxis == 3:
        mosaic_type = 'spectral'
    else:
        log.error('The number of non-trivial axes in the Input FITS files is {} but must be 2 or 3 for MosaicQueen to work. Aborting.'.format(nlongaxis))

    # remove from input image list those images that do not fall within the requested RA,Dec region
    if subimage_dict['CRVAL1'] or subimage_dict['CRVAL2']:
        images = make_mosaic.filter_images_list(images, subimage_dict, input_dir, mosaic_type)
        if not images:
            log.error(
                "None of the input images overlap with the RA,Dec region requested for mosaicking.")
            raise LookupError("None of the input images overlap with the RA,Dec region requested for mosaicking.")
        log.info('Target images overlapping with the RA,Dec region requested for mosaicking:')
        for im in images:
            log.info('    {}'.format(im))

    if (mosaic_type == 'spectral') and (args.velocity and not args.dv):
        log.error(
            "Must define a cube size along the z-axis ('-dv'), if a mosaic central coordinate along this axis was requested ('-v').")
        raise LookupError("Must define a cube size along the z-axis ('-dv'), if a mosaic central coordinate along this axis was requested ('-v').")

    imagesR = [tt.replace('image.fits', 'imageR.fits') for tt in images]
    beams = [tt.replace('image.fits', 'pb.fits') for tt in images]
    beamsR = [tt.replace('image.fits', 'pbR.fits') for tt in images]

    log.info("Will generate a mosaic from the input images.")
    if not args.associated_mosaics:
        log.info("If you would also like mosaics to be made from the associated masks, models and/or residuals, please re-run giving your choice at the '-a' argument.")
        masks, masksR, models, modelsR, residuals, residualsR = [], [], [], [], [], []

    else:
        if 'mask' in args.associated_mosaics:
            log.info('Will also generate a mosaic from the masks associated with the input images')
            masks = [tt.replace('image.fits', 'mask.fits') for tt in images]
            if subimage_dict['CRVAL3'] and mosaic_type == 'spectral' and (args.force_regrid or args.regrid):
                masksR = [tt.replace('image.fits', 'z_cut_maskR.fits') for tt in images]
            else:
                masksR = [tt.replace('image.fits', 'maskR.fits') for tt in images]
        else:
            masks, masksR = [], []

        if 'model' in args.associated_mosaics:
            log.info('Will also generate a mosaic from the models associated with the input images')
            models = [tt.replace('image.fits', 'model.fits') for tt in images]
            if subimage_dict['CRVAL3'] and mosaic_type == 'spectral' and (args.force_regrid or args.regrid):
                modelsR = [tt.replace('image.fits', 'z_cut_modelR.fits') for tt in images]
            else:
                modelsR = [tt.replace('image.fits', 'modelR.fits') for tt in images]
        else:
            models, modelsR = [], []

        if 'residual' in args.associated_mosaics:
            log.info('Will also generate a mosaic from the residuals associated with the input images')
            residuals = [tt.replace('image.fits', 'residual.fits') for tt in images]
            if subimage_dict['CRVAL3'] and mosaic_type == 'spectral' and (args.force_regrid or args.regrid):
                residualsR = [tt.replace('image.fits', 'z_cut_residualR.fits') for tt in images]
            else:
                residualsR = [tt.replace('image.fits', 'residualR.fits') for tt in images]
        else:
            residuals, residualsR = [], []

    log.info('Checking for images and beams')
    check_for_files(input_dir, images, 'images', args.regrid)  # This raises an error and exits if files are not found
    check_for_files(input_dir, beams, 'beams', args.regrid)     # This raises an error and exits if files are not found
    if masks:  # Want to be sure that all of the ingredients are in place before doing any mosaicking
        log.info('Checking for masks')
        check_for_files(input_dir, masks, 'masks', args.regrid)              # This raises an error and exits if files are not found
    if models:  # Want to be sure that all of the ingredients are in place before doing any mosaicking
        log.info('Checking for models')
        check_for_files(input_dir, models, 'models', args.regrid)           # This raises an error and exits if files are not found
    if residuals:  # Want to be sure that all of the ingredients are in place before doing any mosaicking
        log.info('Checking for residuals')
        check_for_files(input_dir, residuals, 'residuals', args.regrid)  # This raises an error and exits if files are not found

    # Check BITPIX of input images
    bitpix = make_mosaic.find_lowest_precision(input_dir, images)

    tmp_inputs = []

    if args.force_regrid:
        log.info('You have asked for all regridded files to be created by this run, even if they are already on disk')
        if subimage_dict['CRVAL3'] and mosaic_type == 'spectral' and (args.force_regrid or args.regrid):
            log.info('Creating spectral slabs for mosaicking ...')
            images, imagesR = make_mosaic.create_spectral_slab(images, input_dir, outname, 'image', subimage_dict)
            beams, beamsR = make_mosaic.create_spectral_slab(beams, input_dir, outname, 'pb', subimage_dict)
            tmp_inputs += images + beams
            print(tmp_inputs)
            log.info('Files we plan to delete later so far are: {}'.format(tmp_inputs))
        make_mosaic.use_montage_for_regridding(
            input_dir, output_dir, mosaic_type, 'image', images, imagesR, beams, beamsR, outname, bitpix, naxis, subimage_dict, num_workers)
        make_mosaic.use_montage_for_regridding(
            input_dir, output_dir, mosaic_type, 'pb', images, imagesR, beams, beamsR, outname, bitpix, naxis, subimage_dict, num_workers)
    elif args.regrid:
        log.info('Checking for regridded images and beams')
        imagesR_dont_exist = check_for_files(output_dir, imagesR, 'regridded images', args.regrid)
        if imagesR_dont_exist:
            log.info('Regridded images are not all in place, so using montage to create them')
            if subimage_dict['CRVAL3'] and mosaic_type == 'spectral' and (args.force_regrid or args.regrid):
                log.info('Creating spectral slabs for mosaicking ...')
                images, imagesR = make_mosaic.create_spectral_slab(images, input_dir, outname, 'image', subimage_dict)
                tmp_inputs += images
            make_mosaic.use_montage_for_regridding(
                input_dir, output_dir, mosaic_type, 'image', images, imagesR, beams, beamsR, outname, bitpix, naxis, subimage_dict, num_workers)
        else:
            log.info('Regridded images are all in place')
        beamsR_dont_exist = check_for_files(output_dir, beamsR, 'regridded beams', args.regrid)
        if beamsR_dont_exist:
            if subimage_dict['CRVAL3'] and mosaic_type == 'spectral' and (args.force_regrid or args.regrid):
                log.info('Creating spectral slabs for mosaicking ...')
                beams, beamsR = make_mosaic.create_spectral_slab(beams, input_dir, outname, 'pb', subimage_dict)
                tmp_inputs += beams
            log.info('Regridded beams are not all in place, so using montage to create them')
            make_mosaic.use_montage_for_regridding(
                input_dir, output_dir, mosaic_type, 'pb', images, imagesR, beams, beamsR, outname, bitpix, naxis, subimage_dict, num_workers)
        #else:  # redundant
        #    log.info('Regridded beams are all in place')
    else:
        log.info('User specified neither --force-regrid nor --regrid in the mosaicqueen command')
        log.info(
            'Will use mosaic header {0:s}.hdr and regridded images and beams available on disk'.format(outname))
        make_mosaic.final_check_for_files(output_dir, imagesR, beamsR)  # This function raises an error and exits if files are not found


    log.info("Mosaicking 'image' files")

    noises = make_mosaic.make_mosaic_using_beam_info(input_dir, output_dir, mosaic_type, 'image', outname, imagesR, beamsR, beam_cutoff, unity_weights, statistic, sigma_guess, images, mosaic_cutoff, bitpix, naxis)


    if 'mask' in args.associated_mosaics:  # Code is more readable by keeping these mosaics separate
        log.info("Mosaicking 'mask' files")
        if args.force_regrid:
            log.info('You have asked for all regridded files to be created by this run, even if they are already on disk')
            if subimage_dict['CRVAL3'] and mosaic_type == 'spectral' and (args.force_regrid or args.regrid):
                log.info('Creating spectral slabs of mask cubes for mosaicking ...')
                masks, masksR = make_mosaic.create_spectral_slab(masks, input_dir, outname, 'mask', subimage_dict)
                tmp_inputs += masks
            make_mosaic.use_montage_for_regridding(
                input_dir, output_dir, mosaic_type, 'mask', masks, masksR, beams, beamsR, outname, bitpix, naxis, subimage_dict, num_workers)
        elif args.regrid:
            log.info('Checking for regridded masks')
            masksR_dont_exist = check_for_files(output_dir, masksR, 'regridded masks', args.regrid)
            if masksR_dont_exist:
                log.info('Regridded masks are not all in place, so using montage to create them')
                if subimage_dict['CRVAL3'] and mosaic_type == 'spectral' and (args.force_regrid or args.regrid):
                    log.info('Creating spectral slabs for mosaicking ...')
                    masks, masksR = make_mosaic.create_spectral_slab(masks, input_dir, outname, 'mask', subimage_dict)
                    tmp_inputs += masks
                make_mosaic.use_montage_for_regridding(
                    input_dir, output_dir, mosaic_type, 'mask', masks, masksR, beams, beamsR, outname, bitpix, naxis, subimage_dict, num_workers)
        else:
            log.info('User specified neither --force-regrid nor --regrid in the mosaicqueen command')
            log.info('Will use regridded masks, available on disk'.format(outname))
            check_for_files(output_dir, masksR, 'regridded masks', args.regrid)
        make_mosaic.make_mosaic_using_beam_info(input_dir, output_dir, mosaic_type, 'mask'    , outname, masksR    , beamsR, beam_cutoff, unity_weights, statistic, sigma_guess, models   , mosaic_cutoff, bitpix, naxis, all_noise_estimates=noises)

    if 'model' in args.associated_mosaics:  # Code is more readable by keeping these mosaics separate
        log.info("Mosaicking 'model' files")
        if args.force_regrid:
            log.info('You have asked for all regridded files to be created by this run, even if they are already on disk')
            if subimage_dict['CRVAL3'] and mosaic_type == 'spectral' and (args.force_regrid or args.regrid):
                log.info('Creating spectral slabs of model cubes for mosaicking ...')
                models, modelsR = make_mosaic.create_spectral_slab(models, input_dir, outname, 'model', subimage_dict)
                tmp_inputs += models
            make_mosaic.use_montage_for_regridding(
                input_dir, output_dir, mosaic_type, 'model', models, modelsR, beams, beamsR, outname, bitpix, naxis, subimage_dict, num_workers)
        elif args.regrid:
            log.info('Checking for regridded models')
            modelsR_dont_exist = check_for_files(output_dir, modelsR, 'regridded models', args.regrid)
            if modelsR_dont_exist:
                log.info('Regridded models are not all in place, so using montage to create them')
                if subimage_dict['CRVAL3'] and mosaic_type == 'spectral' and (args.force_regrid or args.regrid):
                    log.info('Creating spectral slabs for mosaicking ...')
                    models, modelsR = make_mosaic.create_spectral_slab(models, input_dir, outname, 'model', subimage_dict)
                    tmp_inputs += models
                make_mosaic.use_montage_for_regridding(
                    input_dir, output_dir, mosaic_type, 'model', models, modelsR, beams, beamsR, outname, bitpix, naxis, subimage_dict, num_workers)
        else:
            log.info('User specified neither --force-regrid nor --regrid in the mosaicqueen command')
            log.info('Will use regridded models available on disk'.format(outname))
            check_for_files(output_dir, modelsR, 'regridded models', args.regrid)
        make_mosaic.make_mosaic_using_beam_info(input_dir, output_dir, mosaic_type, 'model'   , outname, modelsR   , beamsR, beam_cutoff, unity_weights, statistic, sigma_guess, models   , mosaic_cutoff, bitpix, naxis, all_noise_estimates=noises)

    if 'residual' in args.associated_mosaics:  # Code is more readable by keeping these mosaics separate
        log.info("Mosaicking 'residual' files")
        if args.force_regrid:
            log.info('You have asked for all regridded files to be created by this run, even if they are already on disk')
            if subimage_dict['CRVAL3'] and mosaic_type == 'spectral' and (args.force_regrid or args.regrid):
                log.info('Creating spectral slabs of residual cubes for mosaicking ...')
                residuals, residualsR = make_mosaic.create_spectral_slab(residuals, input_dir, outname, 'residual', subimage_dict)
                tmp_inputs += residuals
            make_mosaic.use_montage_for_regridding(
                input_dir, output_dir, mosaic_type, 'residual', residuals, residualsR, beams, beamsR, outname, bitpix, naxis, subimage_dict, num_workers)
        elif args.regrid:
            log.info('Checking for regridded residuals')
            residualsR_dont_exist = check_for_files(output_dir, residualsR, 'regridded residuals', args.regrid)
            if residualsR_dont_exist:
                log.info('Regridded residuals are not all in place, so using montage to create them')
                if subimage_dict['CRVAL3'] and mosaic_type == 'spectral' and (args.force_regrid or args.regrid):
                    log.info('Creating spectral slabs for mosaicking ...')
                    residuals, residualsR = make_mosaic.create_spectral_slab(residuals, input_dir, outname, 'residual', subimage_dict)
                    tmp_inputs += residuals
                make_mosaic.use_montage_for_regridding(
                    input_dir, output_dir, mosaic_type, 'residual', residuals, residualsR, beams, beamsR, outname, bitpix, naxis, subimage_dict, num_workers)
        else:
            log.info('User specified neither --force-regrid nor --regrid in the mosaicqueen command')
            log.info('Will use regridded residuals available on disk'.format(outname))
            check_for_files(output_dir, residualsR, 'regridded residuals', args.regrid)
        make_mosaic.make_mosaic_using_beam_info(input_dir, output_dir, mosaic_type, 'residual', outname, residualsR, beamsR, beam_cutoff, unity_weights, statistic, sigma_guess, residuals, mosaic_cutoff, bitpix, naxis, all_noise_estimates=noises)

    # Move the log file to the output directory
    os.system('mv log-make_mosaic.txt '+output_dir+'/')

    # Remove temporary spectral slab cubes if there are any
    log.info('Files we plan to delete are: {}'.format(tmp_inputs))
    if tmp_inputs:
        log.info('Removing temporary spectral slab cubes from the input folder ({}) ...'.format(input_dir))
        for z_cut in tmp_inputs:
            to_del = os.path.join(input_dir,z_cut)
            log.info('Trying to delete {}'.format(to_del))
            os.remove(to_del)
    log.info('Done')

    return 0
