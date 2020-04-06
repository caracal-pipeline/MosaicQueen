#!/usr/bin/env python
# -*- coding: utf-8 -*-
# flake8: noqa

# ------------------------------------------------------------------------------------------------------
# Author of package: Sarah White (sarahwhite.astro@gmail.com) and Sphe Makhathini (sphemakh@gmail.com)
# This script based on a convolving script by Landman Bester (landman.bester@gmail.com)
# ------------------------------------------------------------------------------------------------------

import numpy as np
from astropy.io import fits
from africanus.model.spi.examples.utils import load_fits_contiguous, get_fits_freq_space_info, Gaussian2D  ### SO JUST NEED TO ADD AFRICANUS TO REQUIREMENTS?
from pypocketfft import r2c, c2r   ### ADD PYPOCKETFFT TO REQUIREMENTS IN SETUP.PY TOO?
iFs = np.fft.ifftshift
Fs = np.fft.fftshift

log = MosaicSteward.log

# So that error handling is compatible with Python 2 as well as Python 3
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

# -------------------- Turning main.py from the original script into a function ------------------------ #

def convolve_image(input_dir, imagename, beampars, ncpu):
 
    # read resolution of fits file
    hdr = fits.getheader(input_dir+'/'+imagename)
    l_coord, m_coord, freqs, _, freq_axis = get_fits_freq_space_info(hdr)
    nchan = freqs.size
    #psf_pars = {}     ### psf parameters determined in main.py but this may be handy
    #for i in range(1,nchan+1):
    #    key = 'BMAJ' + str(i)
    #    if key in hdr.keys():
    #        emaj = hdr[key]
    #        emin = hdr['BMIN' + str(i)]
    #        pa = hdr['BPA' + str(i)]
    #        psf_pars[i] = (emaj, emin, pa)

    # update header
    for i in range(1, nchan+1):
        hdr['BMAJ' + str(i)] = beampars[0]
        hdr['BMIN' + str(i)] = beampars[1]
        hdr['BPA' + str(i)] = beampars[2]

    # coodinate grid
    xx, yy = np.meshgrid(l_coord, m_coord, indexing='ij') 

    # get padding
    npix_l, npix_m = xx.shape
    pfrac = 0.2/2
    npad_l = int(pfrac*npix_l)
    npad_m = int(pfrac*npix_m)
    
    # get fast FFT sizes and update padding
    from scipy.fftpack import next_fast_len
    nfft = next_fast_len(npix_l + 2*npad_l)
    npad_ll = (nfft - npix_l)//2
    npad_lr = nfft - npix_l - npad_ll
    nfft = next_fast_len(npix_m + 2*npad_m)
    npad_ml = (nfft - npix_m)//2
    npad_mr = nfft - npix_m - npad_ml
    padding = ((0, 0), (npad_ll, npad_lr), (npad_ml, npad_mr))
    unpad_l = slice(npad_ll, -npad_lr)
    unpad_m = slice(npad_ml, -npad_mr)
    
    ax = (1, 2)  # axes over which to perform fft
    lastsize = npix_m + np.sum(padding[-1])

    # kernel of desired resolution
    gausskern = Gaussian2D(xx, yy, beampars)
    gausskern = np.pad(gausskern[None], padding, mode='constant')
    gausskernhat = r2c(iFs(gausskern, axes=ax), axes=ax, forward=True, nthreads=ncpu, inorm=0)

    # FT of image
    image = load_fits_contiguous(input_dir+'/'+imagename)
    image = np.pad(image, padding, mode='constant')
    imhat = r2c(iFs(image, axes=ax), axes=ax, forward=True, nthreads=ncpu, inorm=0)

    # convolve to desired resolution
    if len(psf_pars) == 0:
        imhat *= gausskernhat
    else:
        for i in range(nchan):
            thiskern = Gaussian2D(xx, yy, psf_pars[i+1])
            thiskern = np.pad(thiskern[None], padding, mode='constant')
            thiskernhat = r2c(iFs(thiskern, axes=ax), axes=ax, forward=True, nthreads=ncpu, inorm=0)

            convkernhat = np.where(np.abs(thiskernhat)>0.0, gausskernhat/thiskernhat, 0.0)

            imhat[i] *= convkernhat[0]

    image = Fs(c2r(imhat, axes=ax, forward=False, lastsize=lastsize, inorm=2, nthreads=ncpu), axes=ax)[:, unpad_l, unpad_m]

    ### BEAM CORRECTION DONE AT MOSAICKING STEP, SO REMOVED FROM HERE

    # Save next to input image if no outfile is provided
    tmp = image[::-1]  ### CHECK THIS
    idx = tmp.find('.')
    outfile = args.model[0:-idx]

    hdu = fits.PrimaryHDU(header=hdr)
    # save it
    if freq_axis == 3:
        hdu.data = np.transpose(image, axes=(0, 2, 1))[None, :, :, ::-1].astype(np.float32)
    elif freq_axis == 4:
        hdu.data = np.transpose(image, axes=(0, 2, 1))[:, None, :, ::-1].astype(np.float32)
    name = input_dir + '/' + outfile + '.convolved.fits'
    hdu.writeto(name, overwrite=True)
    log.info("Wrote convolved image: {0:s}".format(name))

    return
    
    ### STUFF ABOUT MULTIPROCESSING MOVED TO MAIN.PY
    
