#!/usr/bin/env python

# ------------------------------------------------------------------------------------------------------
# Author of package: Sarah White (sarahwhite.astro@gmail.com) and Sphe Makhathini (sphemakh@gmail.com)
# Based on a mosaicking script by Paolo Serra (paolo80serra@gmail.com)
# ------------------------------------------------------------------------------------------------------

from astropy.io import fits
import subprocess
from sys import argv
import sys
import os
import numpy as np
import montage_mosaic
import argparse

log = montage_mosaic.log



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
        for hh in 'crpix1,crval1,cdelt1,crpix2,crval2,cdelt2,crpix3,crval3,cdelt3'.split(','):
            astro_t_head[hh] = float(astro_t_head[hh])
        for hh in 'naxis,naxis1,naxis2,naxis3,equinox'.split(','):
            astro_t_head[hh] = int(astro_t_head[hh])
    if mosaic_type == 'continuum':
        for hh in 'crpix1,crval1,cdelt1,crpix2,crval2,cdelt2'.split(','):
            astro_t_head[hh] = float(astro_t_head[hh])
        for hh in 'naxis,naxis1,naxis2,equinox'.split(','):
            astro_t_head[hh] = int(astro_t_head[hh])
    for hh in 'ctype1,ctype2'.split(','):
        astro_t_head[hh] = astro_t_head[hh].replace("'", "")
    return astro_t_head


# ---------------------------- New/re-structured functions --------------------------------- #

def use_montage_for_regridding(input_dir, output_dir, mosaic_type, images, beams, imagesR, beamsR, outname):
    log.info('Running montage tasks to create mosaic header ...')
    # Create an image list
    create_montage_list(images, '{0:s}/{1:s}_fields'.format(output_dir,outname))
    #print(sys.stdout)
    Run('mImgtbl -t {0:s}/{1:s}_fields {2:s} {0:s}/{1:s}_fields.tbl'.format(output_dir,outname,input_dir))
    # Create mosaic header
    Run('mMakeHdr {0:s}/{1:s}_fields.tbl {0:s}/{1:s}.hdr'.format(output_dir,outname))
    
    # Which montage program is used for regridding depends on whether the image is 2D or 3D
    if mosaic_type == 'spectral':
        montage_projection = 'mProjectCube'
        montage_add = 'mAddCube'
    if mosaic_type == 'continuum':
        montage_projection = 'mProject'
        montage_add = 'mAdd'

    log.info('Running montage tasks to regrid the input images ...')
    # Reproject the input images
    for cc in images:
        Run(montage_projection + ' {0:s}/{1:s} {2:s}/{3:s} {2:s}/{4:s}.hdr'.format(
            input_dir, cc, output_dir, cc.replace('image.fits', 'imageR.fits'), outname))
    # Create a reprojected-image metadata file
    create_montage_list(imagesR, '{0:s}/{1:s}_fields_regrid'.format(output_dir,outname))
    Run('mImgtbl -d -t {0:s}/{1:s}_fields_regrid {0:s} {0:s}/{1:s}_fields_regrid.tbl'.format(output_dir,outname)) # '-d' flag added to aid de-bugging
    # Co-add the reprojected images
    #Run(montage_add + ' -p . {0:s}_fields_regrid.tbl {0:s}.hdr {0:s}.fits'.format(outname))
    
    log.info('Running montage tasks to regrid the input beams ...')
    # Reproject the input beams
    for bb in beams:
        Run(montage_projection + ' {0:s}/{1:s} {2:s}/{3:s} {2:s}/{4:s}.hdr'.format(
            input_dir, bb, output_dir, bb.replace('pb.fits', 'pbR.fits'), outname))
    # Create a reprojected-beams metadata file
    create_montage_list(beamsR, '{0:s}/{1:s}_beams_regrid'.format(output_dir,outname))
    Run('mImgtbl -t {0:s}/{1:s}_beams_regrid {0:s} {0:s}/{1:s}_beams_regrid.tbl'.format(output_dir,outname))
    # Co-add the reprojected beams
    #Run(montage_add + ' -p . {0:s}_beams_regrid.tbl {0:s}.hdr {0:s}pb.fits'.format(outname))

    return 0


def check_for_regridded_files(output_dir, imagesR, beamsR):
    # As the regridded files were produced by montage_mosaic, we expect them to be in the output directory

    for cc in imagesR:
        if not os.path.exists(output_dir+'/'+cc):
            log.error('File {0:s} does not exist'.format(output_dir+'/'+cc))
            sys.exit()

    for bb in beamsR:
        if not os.path.exists(output_dir+'/'+bb):
            log.error('File {0:s} does not exist'.format(output_dir+'/'+bb))
            sys.exit()

    log.info('All regridded images and beams found on disc')

    return 0


def make_mosaic_using_beam_info(input_dir, output_dir, mosaic_type, outname, imagesR, beamsR, cutoff, images):

    log.info('Mosaicking ...')
    moshead = [jj.strip().replace(' ', '').split('=')
            for jj in open('{0:s}/{1:s}.hdr'.format(output_dir,outname)).readlines()]
    if ['END'] in moshead:
        del(moshead[moshead.index(['END'])])
    moshead = {k: v for (k, v) in moshead}   # Creating a dictionary, where 'k' stands for 'keyword' and 'v' stands for 'value'
    if mosaic_type == 'spectral':
        mos_array = np.zeros((int(moshead['NAXIS3']), int(
            moshead['NAXIS2']), int(moshead['NAXIS1'])), dtype='float32')
        norm_array = np.zeros((int(moshead['NAXIS3']), int(
            moshead['NAXIS2']), int(moshead['NAXIS1'])), dtype='float32')
    if mosaic_type == 'continuum':
        mos_array = np.zeros((int(moshead['NAXIS2']), int(moshead['NAXIS1'])), dtype='float32')
        norm_array = np.zeros((int(moshead['NAXIS2']), int(moshead['NAXIS1'])), dtype='float32')
    for ii, bb in zip(imagesR, beamsR):
        log.info('Adding {0:s} to the mosaic ...'.format(ii))
        f = fits.open(output_dir+'/'+ii)  # i.e. open a specific re-gridded image
        head = f[0].header
        g = fits.open(output_dir+'/'+bb)  # i.e. open a specific re-gridded beam
        y1 = int(float(moshead['CRPIX2']) - head['CRPIX2'])
        y2 = int(float(moshead['CRPIX2']) - head['CRPIX2'] + head['NAXIS2'])
        x1 = int(float(moshead['CRPIX1']) - head['CRPIX1'])
        x2 = int(float(moshead['CRPIX1']) - head['CRPIX1'] + head['NAXIS1'])
        # mosaicking with no PB correction
        # norm_array[:,y1:y2,x1:x2]+=~np.isnan(f[0].data)
        # mos_array[:,y1:y2,x1:x2]+=np.nan_to_num(f[0].data)
        # mosaicking with PB correction
        if mosaic_type == 'spectral':
            norm_array[:, y1:y2, x1:x2] += (np.nan_to_num(g[0].data)
                                      * (np.nan_to_num(g[0].data) > cutoff))**2
            mos_array[:, y1:y2, x1:x2] += np.nan_to_num(f[0].data)*np.nan_to_num(
                g[0].data)*(np.nan_to_num(g[0].data) > cutoff)
        if mosaic_type == 'continuum':
            norm_array[y1:y2, x1:x2] += (np.nan_to_num(g[0].data)
                                      * (np.nan_to_num(g[0].data) > cutoff))**2
            mos_array[y1:y2, x1:x2] += np.nan_to_num(f[0].data)*np.nan_to_num(
                g[0].data)*(np.nan_to_num(g[0].data) > cutoff)
        f.close()
        g.close()

    norm_array[norm_array == 0] = np.nan
    moshead = make_mosaic_header(mosaic_type, moshead)
    if mosaic_type == 'spectral':
        f = fits.open(input_dir+'/'+images[0]) 
        zhead = f[0].header
        moshead['ctype3'] = zhead['ctype3']
    fits.writeto('{0:s}/{1:s}.fits'.format(output_dir,outname), mos_array /
                 norm_array, overwrite=True, header=moshead)
    fits.writeto('{0:s}/{1:s}_noise.fits'.format(output_dir,outname), 1. /
                 np.sqrt(norm_array), overwrite=True, header=moshead)
    fits.writeto('{0:s}/{1:s}_weights.fits'.format(output_dir,outname),
                 np.sqrt(norm_array), overwrite=True, header=moshead)

    log.info('The following mosaic FITS were written to disc: {0:s}.fits {0:s}_noise.fits {0:s}_weights.fits'.format(outname))
    log.info('Mosaicking completed')

    return 
