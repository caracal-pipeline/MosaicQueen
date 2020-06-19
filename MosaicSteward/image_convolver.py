#!/usr/bin/env python
# -*- coding: utf-8 -*-
# flake8: noqa

# ------------------------------------------------------------------------------------------------------
# Author of package: Sarah White (sarahwhite.astro@gmail.com) and Sphe Makhathini (sphemakh@gmail.com)
# This script is based on a convolving script by Landman Bester (landman.bester@gmail.com)
# ------------------------------------------------------------------------------------------------------

import MosaicSteward
import numpy as np
from astropy.io import fits
from ducc0.fft import r2c, c2r  ### replaced pypocketfft with ducc.fft
iFs = np.fft.ifftshift
Fs = np.fft.fftshift

log = MosaicSteward.log

# So that error handling is compatible with Python 2 as well as Python 3
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError


# -------------------- Functions from codex-africanus (model.spi.examples.utils) ----------------------- #

def load_fits_contiguous(name):

    arr = fits.getdata(name)
    if arr.ndim == 4: 
        # figure out which axes to keep from header
        hdr = fits.getheader(name)
        if hdr["CTYPE4"].lower() == 'stokes':
            arr = arr[0]
        else:
            arr = arr[:, 0]
    # reverse last spatial axis then transpose (f -> c contiguous)
    arr = np.transpose(arr[:, :, ::-1], axes=(0, 2, 1))
    
    return np.ascontiguousarray(arr, dtype=np.float64)


def get_fits_freq_space_info(hdr):
    
    if hdr['CUNIT1'].lower() != "deg":
        log.error("Image coordinates must be in degrees")
        raise ValueError("Image coordinates must be in degrees")
    npix_l = hdr['NAXIS1']
    refpix_l = hdr['CRPIX1']
    delta_l = hdr['CDELT1']
    l_coord = np.arange(1 - refpix_l, 1 + npix_l - refpix_l)*delta_l

    if hdr['CUNIT2'].lower() != "deg":
        log.error("Image coordinates must be in degrees")
        raise ValueError("Image coordinates must be in degrees")
    npix_m = hdr['NAXIS2']
    refpix_m = hdr['CRPIX2']
    delta_m = hdr['CDELT2']
    m_coord = np.arange(1 - refpix_m, 1 + npix_m - refpix_m)*delta_m

    # get frequencies
    if hdr["CTYPE4"].lower() == 'freq':
        freq_axis = 4
        nband = hdr['NAXIS4']
        refpix_nu = hdr['CRPIX4']
        delta_nu = hdr['CDELT4']  # assumes units are Hz
        ref_freq = hdr['CRVAL4']
        ncorr = hdr['NAXIS3']
    elif hdr["CTYPE3"].lower() == 'freq':
        freq_axis = 3
        nband = hdr['NAXIS3']
        refpix_nu = hdr['CRPIX3']
        delta_nu = hdr['CDELT3']  # assumes units are Hz
        ref_freq = hdr['CRVAL3']
        ncorr = hdr['NAXIS4']
    else:
        log.error("Freq axis must be 3rd or 4th")
        raise ValueError("Freq axis must be 3rd or 4th")

    if ncorr > 1:
        log.error("Only Stokes-I cubes supported")
        raise ValueError("Only Stokes-I cubes supported")

    freqs = ref_freq + np.arange(1 - refpix_nu,
                                 1 + nband - refpix_nu) * delta_nu

    return l_coord, m_coord, freqs, ref_freq, freq_axis


def Gaussian2D(xin, yin, GaussPar=(1., 1., 0.)):
    
    S0, S1, PA = GaussPar
    Smaj = np.maximum(S0, S1)
    Smin = np.minimum(S0, S1)
    A = np.array([[1. / Smin ** 2, 0],
                  [0, 1. / Smaj ** 2]])

    c, s, t = np.cos, np.sin, np.deg2rad(-PA)
    R = np.array([[c(t), -s(t)],
                  [s(t), c(t)]])
    A = np.dot(np.dot(R.T, A), R)
    sOut = xin.shape
    # only compute the result out to 5 * emaj
    extent = (5 * Smaj)**2
    xflat = xin.squeeze()
    yflat = yin.squeeze()
    ind = np.argwhere(xflat**2 + yflat**2 <= extent).squeeze()
    idx = ind[:, 0]
    idy = ind[:, 1]
    x = np.array([xflat[idx, idy].ravel(), yflat[idx, idy].ravel()])
    R = np.einsum('nb,bc,cn->n', x.T, A, x)
    # need to adjust for the fact that GaussPar corresponds to FWHM
    fwhm_conv = 2*np.sqrt(2*np.log(2))
    tmp = np.exp(-fwhm_conv*R)
    gausskern = np.zeros(xflat.shape, dtype=np.float64)
    gausskern[idx, idy] = tmp

    return np.ascontiguousarray(gausskern.reshape(sOut),
                                dtype=np.float64)


# ---------------------------- A new function for uniform-resolution mode ------------------------------- #

def find_largest_BMAJ(input_dir, images, mosaic_type, data_type): 

    # data_type is to allow this to be run over 'images' and 'regridded-images' and still give useful log messages

    # To save going through the headers of the images again, keep a record of the beam keywords
    # Will need to be 2D arrays, where the second axis is given by the number of images used as input

    n_images = len(images)

    if mosaic_type == 'continuum':

        log.info('Checking the synthesised-beam information for continuum {0:s}'.format(data_type))
        BMAJ_array = np.empty([n_images, 1]) # Looks strange but to be in parallel with arrays for spectral mode
        BMIN_array = np.empty([n_images, 1])
        BPA_array = np.empty([n_images, 1])

        index_image = 0
        for ii in images:
            
            f = fits.open(input_dir+'/'+ii)
            head = f[0].header
            BMAJ_array[index_image][0] = float(head['BMAJ'])
            BMIN_array[index_image][0] = float(head['BMIN'])
            BPA_array[index_image][0] = float(head['BPA'])
            index_image += 1

        #print('BMAJ_array = ', BMAJ_array)  ### FOR CHECKING 
        #print('BMIN_array = ', BMIN_array)
        #print('BPA_array = ', BPA_array)

    else:

        log.info('Checking the synthesised-beam information for spectral {0:s}'.data_type)

        f = fits.open(input_dir+'/'+images[0])  # i.e. open the first input image
        head = f[0].header
        n_channels = int(head['NAXIS'])
        log.info('Based on the first input, we expect each of the input {0:s} to have {1:i} channels'.format(data_type,n_channels))

        BMAJ_array = np.empty([n_images, n_channels]) 
        BMIN_array = np.empty([n_images, n_channels])
        BPA_array = np.empty([n_images, n_channels])

        index_image = 0
        for ii in images:
            
            f = fits.open(input_dir+'/'+ii)
            head = f[0].header



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

    ## kernel of desired resolution  ### ORIGINAL LOCATION
    #gausskern = Gaussian2D(xx, yy, beampars)
    #gausskern = np.pad(gausskern[None], padding, mode='constant')
    #gausskernhat = r2c(iFs(gausskern, axes=ax), axes=ax, forward=True, nthreads=ncpu, inorm=0)

    # FT of image
    image = load_fits_contiguous(input_dir+'/'+imagename)
    image = np.pad(image, padding, mode='constant')
    imhat = r2c(iFs(image, axes=ax), axes=ax, forward=True, nthreads=ncpu, inorm=0)

    # convolve to desired resolution
    if len(psf_pars) == 0:
        imhat *= gausskernhat
    else:
        for i in range(nchan):

            # kernel of desired resolution   ### WHY MOVED HERE? GAUSSKERNHAT IS NOW DEFINED TOO LATE
            gausskern = Gaussian2D(xx, yy, beampars)
            gausskern = np.pad(gausskern[None], padding, mode='constant')
            gausskernhat = r2c(iFs(gausskern, axes=ax), axes=ax, forward=True, nthreads=ncpu, inorm=0)

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
    
