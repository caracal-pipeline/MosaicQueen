#!/usr/bin/env python

# ------------------------------------------------------------------------------------------------------
# Author of package: Sarah White (sarahwhite.astro@gmail.com) and Sphe Makhathini (sphemakh@gmail.com)
# Based on a mosaicking script by Paolo Serra (paolo.serra@inaf.it)
# ------------------------------------------------------------------------------------------------------

from astropy.io import fits
import subprocess
from sys import argv
import sys
import os
import numpy as np
import mosaicQueen
import argparse
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from memory_profiler import profile

log = mosaicQueen.log

# So that error handling is compatible with Python 2 as well as Python 3
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

# -------------------- Edited functions from the original script ------------------------ #

def create_montage_list(inputfiles, outputfile):
    list_file = open(outputfile, 'w')
    # Need to know length of the longest filename so that none of them get truncated
    max_filename_length = len( max(inputfiles, key=len) )
    start_of_line = '|'
    start_and_whitespace = start_of_line.ljust(max_filename_length) # This will add trailing whitespace of the required length
    list_file.write( start_and_whitespace+'fname|\n' )
    list_file.write( start_and_whitespace+' char|\n' )
    for ff in inputfiles:
        list_file.write(' {0:s}\n'.format(ff))
    list_file.close()


def Run(command, verb1=1, verb2=0, getout=0):
    if verb1:
        log.info('      '+command)
    result = subprocess.check_output(command.split())
    if verb2:
        for jj in result:
            log.info(jj)
    if getout:
        return result


def make_mosaic_header(mosaic_type, t_head):
    astro_t_head = fits.Header()
    for hh in t_head.keys():
        astro_t_head[hh] = t_head[hh]
    if mosaic_type == 'spectral':
        for hh in 'crpix1,crval1,cdelt1,crpix2,crval2,cdelt2,crota2,crpix3,crval3,cdelt3'.split(','):
            astro_t_head[hh] = float(astro_t_head[hh])
        for hh in 'naxis,naxis1,naxis2,naxis3,equinox'.split(','):
            astro_t_head[hh] = int(astro_t_head[hh])
    if mosaic_type == 'continuum':
        for hh in 'crpix1,crval1,cdelt1,crpix2,crval2,cdelt2,crota2'.split(','):
            astro_t_head[hh] = float(astro_t_head[hh])
        for hh in 'naxis,naxis1,naxis2,equinox'.split(','):
            astro_t_head[hh] = int(astro_t_head[hh])
    for hh in 'ctype1,ctype2'.split(','):
        astro_t_head[hh] = astro_t_head[hh].replace("'", "")
    return astro_t_head


# ---------------------------- New/re-structured functions --------------------------------- #

def use_montage_for_regridding(input_dir, output_dir, mosaic_type, image_type, images, imagesR, beams, beamsR, outname, bitpix):
                               # image_type should be 'image', 'pb', 'model', or 'residual'

    dtype = f"float{bitpix}"

    # Which montage program is used for regridding depends on whether the image is 2D or 3D
    if mosaic_type == 'spectral':
        montage_projection = 'mProjectCube'
        montage_add = 'mAddCube'
    elif mosaic_type == 'continuum':
        montage_projection = 'mProject'
        montage_add = 'mAdd'

    if image_type != 'pb':  # i.e. creating a header for 'image', 'model', or 'residual'
        log.info('Running montage tasks to create mosaic header ...')
        # Create an image list
        create_montage_list(images, '{0:s}/{1:s}_{2:s}_fields'.format(output_dir,outname,image_type))
        #print(sys.stdout)
        Run('mImgtbl -t {0:s}/{1:s}_{2:s}_fields {3:s} {0:s}/{1:s}_{2:s}_fields.tbl'.format(output_dir,outname,image_type,input_dir))
        # Create mosaic header
        Run('mMakeHdr {0:s}/{1:s}_{2:s}_fields.tbl {0:s}/{1:s}_{2:s}.hdr'.format(output_dir,outname,image_type))
        log.info('Running montage tasks to regrid files ...')
        # Reproject the input images
        for cc in images:
            Run(montage_projection + ' {0:s}/{1:s} {2:s}/{3:s} {2:s}/{4:s}_{5:s}.hdr'.format(
                input_dir, cc, output_dir, cc.replace(image_type+'.fits', image_type+'R.fits'), outname, image_type))
            # CONVERT FROM 64-bit TO 32-bit HERE
            with fits.open('{0:s}/{1:s}'.format(output_dir, cc.replace(image_type+'.fits', image_type+'R.fits'))) as Rfits:
                head = Rfits[0].header
                if head['bitpix'] != -bitpix:
                    log.info('      Convert from {}-bit to {}-bit ...'.format(np.abs(head['bitpix']), bitpix))
                    head['bitpix'] = -bitpix
                    fits.writeto('{0:s}/{1:s}'.format(output_dir, cc.replace(image_type+'.fits', image_type+'R.fits')), Rfits[0].data.astype(dtype), header=head, overwrite=True)
        # Create a reprojected-image metadata file
        create_montage_list(imagesR, '{0:s}/{1:s}_{2:s}_fields_regrid'.format(output_dir,outname, image_type))
        Run('mImgtbl -d -t {0:s}/{1:s}_{2:s}_fields_regrid {0:s} {0:s}/{1:s}_{2:s}_fields_regrid.tbl'.format(output_dir,outname,image_type)) # '-d' flag added to aid de-bugging
        # Co-add the reprojected images
        #Run(montage_add + ' -p . {0:s}_{1:s}_fields_regrid.tbl {0:s}_{1:s}.hdr {0:s}.fits'.format(outname,image_type))

    else:  # i.e. for image_type == 'pb', to maintain the familiar _beams_regrid filenames
        log.info('Running montage tasks to regrid beams ...')
        # Reproject the input beams
        for bb in beams:
            # Assuming that an image_type == 'image' run of this function was called beforehand
            Run(montage_projection + ' {0:s}/{1:s} {2:s}/{3:s} {2:s}/{4:s}_image.hdr'.format(
                input_dir, bb, output_dir, bb.replace('pb.fits', 'pbR.fits'), outname))
            # CONVERT FROM 64-bit TO 32-bit HERE
            with fits.open('{0:s}/{1:s}'.format(output_dir, bb.replace('pb.fits', 'pbR.fits'))) as Rfits:
                head = Rfits[0].header
                if head['bitpix'] != -bitpix:
                    log.info('      Convert from {}-bit to {}-bit ...'.format(np.abs(head['bitpix']), bitpix))
                    head['bitpix'] = -bitpix
                    fits.writeto('{0:s}/{1:s}'.format(output_dir, bb.replace('pb.fits', 'pbR.fits')), Rfits[0].data.astype(dtype), header=head, overwrite=True)
        # Create a reprojected-beams metadata file
        create_montage_list(beamsR, '{0:s}/{1:s}_beams_regrid'.format(output_dir,outname))
        Run('mImgtbl -t {0:s}/{1:s}_beams_regrid {0:s} {0:s}/{1:s}_beams_regrid.tbl'.format(output_dir,outname))
        # Co-add the reprojected beams
        #Run(montage_add + ' -p . {0:s}_beams_regrid.tbl {0:s}.hdr {0:s}_pb.fits'.format(outname))

    return 0


def final_check_for_files(directory, imagesR, beamsR):
    # As the regridded files were produced by montage_mosaic, we expect them to be in the output directory

    for cc in imagesR:
        try:
            open(directory+'/'+cc)
        except FileNotFoundError:
            log.error('File {0:s} does not exist'.format(directory+'/'+cc))
            raise FileNotFoundError

    for bb in beamsR:
        try:
            open(directory+'/'+bb)
        except FileNotFoundError:
            log.error('File {0:s} does not exist'.format(directory+'/'+bb))
            raise FileNotFoundError

    log.info('All files found on disk')

    return 0


def gauss(x, *p):  # Define model function to be used to fit to the data in estimate_noise, below
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))


def estimate_noise(image_regrid_hdu, statistic, sigma_guess, check_Gaussian_filename):

    if statistic == 'mad':
        log.info('... using the median absolute deviation of all negative pixels (assuming median = 0) ...')
        image_tmp = image_regrid_hdu[0].data
        image_noise_estimate = 1.4826 * np.median(np.abs(image_tmp[(image_tmp < 0) * (~np.isnan(image_tmp))]))
        log.info('... noise estimate = {0:.3e} Jy/beam'.format(image_noise_estimate)) # Assumed units

    if statistic == 'rms':
        log.info('... using the rms of all negative pixels ...')
        image_tmp = image_regrid_hdu[0].data
        image_noise_estimate = np.sqrt(np.nanmean(image_tmp[image_tmp < 0]**2))
        log.info('... noise estimate = {0:.3e} Jy/beam'.format(image_noise_estimate)) # Assumed units

    elif statistic == 'fit':

        log.info('... using a Gaussian fit to the negative values ...')

        image_tmp = np.nan_to_num(image_regrid_hdu[0].data)
        mask = image_tmp < 0.0
        negative_values = image_tmp[mask]
        positive_values = -1.0*negative_values # Flipping to get the other side of the Gaussian
        values = np.append( negative_values, positive_values )

        n, bin_edges, patches = plt.hist(values, bins=100, density=True, facecolor='lightblue')
        bin_centres = (bin_edges[:-1] + bin_edges[1:])/2

        # p0 is the initial guess for the fitting coefficients of the Gaussian (A, mu and sigma above)
        p0 = [np.max(n), 0., sigma_guess]
        #log.info('sigma_guess = ' + str(sigma_guess))

        coeff, var_matrix = curve_fit(gauss, bin_centres, n, p0=p0)

        # Get the fitted curve
        n_fit = gauss(bin_centres, *coeff)

        # Plot to check
        plt.plot(bin_centres, n, label='Image data', color='blue')
        plt.plot(bin_centres, n_fit, label='Fitted line', color='red')
        plt.legend(loc='upper right')
        plt.xlabel('pixel value')
        plt.ylabel('number')
        plt.savefig(check_Gaussian_filename, dpi=72, bbox_inches='tight')
        plt.close() # Close figure, so that lines don't remain for subsequent calls

        # Get the fitting parameters, i.e. the mean and standard deviation:
        log.info('... noise estimate = {0:.3e} Jy/beam'.format(coeff[2])) # Assumed units
        log.info('    (see ' + check_Gaussian_filename + ')')
        image_noise_estimate = coeff[2]

    return image_noise_estimate   # This returns a single value


def update_norm(norm, slc, beam_regrid_hdu, cutoff, noise):
    """
        update normalization array
    """

    tmp = np.nan_to_num(beam_regrid_hdu[0].data)
    mask = np.nan_to_num(beam_regrid_hdu[0].data)>cutoff
    tmp = tmp*tmp / noise**2
    # See update_mos below for an explanation of the masking
    norm[slc] += tmp*mask


def update_mos(mos, slc, image_regrid_hdu, beam_regrid_hdu, cutoff, noise, finite):
    """
        update mosaic array
    """

    weighted_image_tmp = image_regrid_hdu[0].data
    beam_tmp = beam_regrid_hdu[0].data
    mask = beam_tmp > cutoff
    # Set finite = True only for pixels 1) above the beam cutoff and 2) !=NaN in both cube and beam
    finite[slc] += mask * ~np.isnan(beam_tmp) * ~np.isnan(weighted_image_tmp)
    # Pixels not meeting the above requirements 1) and 2) for this cube and beam are temporarily
    # Set to zero before being added to mos and norm array
    mos[slc] +=  np.nan_to_num(weighted_image_tmp) * np.nan_to_num(beam_tmp) * mask / noise**2


def find_lowest_precision(input_dir, images):
    """
        find lowest floating-point precision of the input, and return so that it can be used to set the precision of the output
    """
    bitpix_list = []
    for image in images:
        with fits.open(os.path.join(input_dir,image)) as hdul:
            bitpix_list.append( abs(int(hdul[0]._bitpix)) )
    bitpix = min(bitpix_list)

    return bitpix


#@profile
def make_mosaic_using_beam_info(input_dir, output_dir, mosaic_type, image_type, outname, imagesR, beamsR, cutoff, uwei, statistic, sigma_guess, images, mos_cutoff, bitpix, all_noise_estimates=[]):


    log.info("Creating a mosaic from '{0:s}' files ...".format(image_type))

    # Details for the mosaic header
    moshead = [jj.strip().replace(' ', '').split('=')
            for jj in open('{0:s}/{1:s}_{2:s}.hdr'.format(output_dir,outname,image_type)).readlines()]
    if ['END'] in moshead:
        del(moshead[moshead.index(['END'])])
    moshead = {k: v for (k, v) in moshead}   # Creating a dictionary, where 'k' stands for 'keyword' and 'v' stands for 'value'

    dtype = f"float{bitpix}"
    # delete BITPIX from montage (always 64-bit) so that precision is from input FITS files
    del moshead['BITPIX']

    # Initialise zero-valued mosaic and normalisation arrays
    if mosaic_type == 'spectral':
        shape = ((int(moshead['NAXIS3']), int(moshead['NAXIS2']), int(moshead['NAXIS1'])))
    if mosaic_type == 'continuum':
        shape = (int(moshead['NAXIS2']), int(moshead['NAXIS1']))

    mos_array = np.zeros(shape, dtype=dtype)
    norm_array = np.zeros(shape, dtype=dtype)
    finite_array = np.zeros(shape, dtype='bool')

    # Gathering noise estimates for each of the input images
    if uwei:
        log.info("Will use weight=1 instead of weight=1/noise**2 for all images")
        all_noise_estimates = [1. for image in imagesR]
    else:
        log.info("Will use 1/noise**2 weights")
        if image_type == 'image':
            for image in imagesR:
                log.info('Estimating the noise level of {0:s} ...'.format(image))
                image_regrid_hdu = fits.open(output_dir+'/'+image, mmap=True)  # i.e. open a specific re-gridded image
                check_Gaussian_filename = output_dir + '/' + image.replace('.fits', '_check_Gaussian_fit.png')
                image_noise_estimate = estimate_noise(image_regrid_hdu, statistic, sigma_guess, check_Gaussian_filename)
                all_noise_estimates.append(image_noise_estimate)
        log.info("Summary of noise levels estimated from the 'image' files:")
        for ee in all_noise_estimates:
            log.info('    {0:.3e} Jy/beam'.format(ee)) # Assumed units

    # The mosaicking part: iterate over input regridded arrays and add to mosaic and normalisation arrays at each step
    weighting_index = 0

    log.info('Adding single fields to mosaic using primary beam cutoff = {0:.3f} (set by --beam-cutoff)'.format(cutoff))
    for ii, bb, ss in zip(imagesR, beamsR, all_noise_estimates):
        log.info('Adding {0:s} to the mosaic ...'.format(ii))
        image_regrid_hdu = fits.open(output_dir+'/'+ii, mmap=True)  # i.e. open a specific re-gridded image
        head = image_regrid_hdu[0].header
        beam_regrid_hdu = fits.open(output_dir+'/'+bb, mmap=True)  # i.e. open a specific re-gridded beam

        # Calculate the location of this regridded input array within the mosaic array
        y1 = int(float(moshead['CRPIX2']) - head['CRPIX2'])
        y2 = int(float(moshead['CRPIX2']) - head['CRPIX2'] + head['NAXIS2'])
        x1 = int(float(moshead['CRPIX1']) - head['CRPIX1'])
        x2 = int(float(moshead['CRPIX1']) - head['CRPIX1'] + head['NAXIS1'])
        if mosaic_type == 'spectral':
            slc = slice(None), slice(y1,y2), slice(x1,x2)
        elif mosaic_type == 'continuum':
            slc = slice(y1,y2), slice(x1,x2)

        # Add the regridded input array with appropriate weights
        update_norm(norm_array, slc, beam_regrid_hdu, cutoff, ss)
        update_mos(mos_array, slc, image_regrid_hdu, beam_regrid_hdu , cutoff, ss, finite_array)
        weighting_index = weighting_index + 1
        image_regrid_hdu.close()
        beam_regrid_hdu.close()

    # Blank pixels below the sensitivity cutoff = mos_cutoff in the final mosaic
    # Since sigma = 1/sqrt(norm), then pixels should be blanked if:
    #   sigma > sigma_min / mos_cutoff  ->  norm < norm_max * mos_cutoff**2
    log.info('Blanking mosaic pixels with noise level > minimum mosaic noise level / {0:.3f} (set by --mosaic_cutoff)'.format(mos_cutoff))
    finite_array[norm_array < np.nanmax(norm_array) * mos_cutoff**2] = False

    # Pixels whose value is False in the finite array are blanked in the final mos and norm arrays
    mos_array[~finite_array] = np.nan
    norm_array[~finite_array] = np.nan

    # Fixing mosaic header (add missing keys that exist in the original input cubes but not yet in the mosaic cube)
    # Header keys which should be identical in all input cubes (if not, we take the first one)
    single_keys = [
                   'bunit',
                   'specsys',
                   'specsys3',
                   'restfreq',
                   'restfrq',
                   'velref',
                   'epoch',
                   'equinox',
                   'altrval',
                   'altrpix',
                   'telescop',
                   'observer',
                   'obsgeo-x',
                   'obsgeo-y',
                   'obsgeo-z',
                   'cellscal',
                   ]

    # header keys which can vary from one input cube to another (we take the median)
    multi_keys  = 'bmaj bmin bpa'.split()
    all_keys    = single_keys + multi_keys
    for kk in all_keys:
        if kk not in moshead:
            keyvals = []
            for ii in images:
                with fits.open('{0}/{1}'.format(input_dir,ii)) as org_hdu:
                    org_head = org_hdu[0].header
                if kk in org_head:
                    keyvals.append(org_head[kk])
            if len(keyvals):
                if kk in single_keys:
                    if np.unique(np.array(keyvals)).shape[0] > 1:
                        log.info('Inconsistent {0} values in input cubes. Will use {0} of the first input cube.'.format(kk))
                    moskey = keyvals[0]
                elif kk in multi_keys:
                    moskey = np.median(np.array(keyvals))
                log.info('Setting {0} = {1} in mosaic FITS header'.format(kk, moskey))
                moshead[kk] = moskey

    # Tidying up and writing
    moshead = make_mosaic_header(mosaic_type, moshead)
    if mosaic_type == 'spectral':
        f = fits.open(input_dir+'/'+images[0])
        zhead = f[0].header
        moshead['ctype3'] = zhead['ctype3']
        f.close()
    fits.writeto('{0:s}/{1:s}_{2:s}.fits'.format(output_dir,outname,image_type), mos_array /
                 norm_array, overwrite=True, header=moshead)

    # Creating the accompanying 'noise' and 'weights' mosaics
    if image_type == 'image':  # Only want one copy of the _noise and _weights mosaics to be produced
        fits.writeto('{0:s}/{1:s}_noise.fits'.format(output_dir,outname), 1. /
                     np.sqrt(norm_array), overwrite=True, header=moshead)
        fits.writeto('{0:s}/{1:s}_weights.fits'.format(output_dir,outname),
                     np.sqrt(norm_array), overwrite=True, header=moshead)
        log.info('The following mosaic FITS were written to disk: {0:s}_{1:s}.fits {0:s}_noise.fits {0:s}_weights.fits'.format(outname,image_type))
    else:  # i.e. when making a mosaic of the 'model' or 'residual' .fits files
        log.info('The following mosaic FITS was written to disk: {0:s}_{1:s}.fits'.format(outname,image_type))

    log.info("Mosaicking of '{}' files completed".format(image_type))

    return all_noise_estimates
