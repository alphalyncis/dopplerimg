import numpy as np
from astropy.table import Table
import modelfitting as mf
import sys
from scipy import signal
import pickle
import glob
from astropy.io import fits
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import os
homedir = os.path.expanduser('~')

def fit(modelpath):
    """
    Fit avged observations, but each obs with own wcoef.
    """
    ###################################
    #  Open data
    ###################################

    filelist = sorted(glob.glob(f'{homedir}/uoedrive/data/IGRINS/SDCK*_1f.spec.fits'))

    fluxes = []
    wls = []

    for filename in filelist:
        wlname = filename.split('_1f')[0]+'.wave.fits'

        flux = fits.getdata(filename)
        wl = fits.getdata(wlname)

        # trim first and last 100 columns
        flux = flux[:, 100:1948]
        wl = wl[:, 100:1948]

        #print(wl.shape)
        
        fluxes.append(flux)
        wls.append(wl)

    ###############################
    ## Collapse data over the whole observation into one mean, smoothed spectrum
    ###############################

    dims = np.array(fluxes).shape
    fluxes = np.array(fluxes)
    wls = np.array(wls)

    # remove first order with all NaNs when making various arrays
    obs0 = np.median(fluxes[:, 1:, :], axis=0)*dims[0]  # make mean spectrum
    eobs0 = np.median(fluxes[:, 1:, :], axis=0)*np.sqrt(dims[0])  # make mean noise spectrum (assuming just photon noise)
    fobs0 = np.vstack([signal.medfilt(obs0[jj], 3) for jj in range(dims[1]-1)])  # smooth spectrum
    eobs0 /= np.nanmedian(fobs0, 1).reshape(dims[1]-1,1)
    fobs0 /= np.nanmedian(fobs0, 1).reshape(dims[1]-1,1)
    wfobs0 = 1./eobs0**2   # make noise spectrum
    #ind = np.where(np.isinf(wfobs0)) # set infinite values in noise spectrum to NaNs
    #wfobs0[ind] = np.nan
    wind = np.isinf(wfobs0)  # remove points with infinite values
    wfobs0[wind] = 0.

    # fix wavelength array to have same shape
    wls = wls[:, 1:24, :] 

    ##########################
    # open model
    ##########################

    model = Table.read(modelpath, format='fits')
    modelname = modelpath.split("/")[-1]
    model['wl'] = model['Wavelength']
    model['flux'] = model['Flux']

    # set # of observations and # of orders to process
    nobs = wls.shape[0]
    norders = 20  # skip last few non-fitting orders for now
    npix = dims[2]

    NPW = 4  # number of terms for polynomial fit to the wavelengths
    pix = np.arange(npix, dtype=float)/npix
    chipfits = []
    chipmods = np.zeros((nobs, norders, npix), dtype=float)
    chiplams = np.zeros((nobs, norders, npix), dtype=float)
    chipmodnobroad = np.zeros((nobs, norders, npix), dtype=float)
    chipguesses = np.zeros((nobs, norders, npix), dtype=float)
    chisqarr = np.zeros((nobs, norders), dtype=float)

    # set up tables for best fit vsini, limb darkening values, and rv
    orderval=[]
    obsval=[]
    vsini = []
    limbdark = []
    rv = []
    chisq = []

    for jj in (np.arange(norders)):
        lolim = wls[:, jj, :].min() - 0.003
        hilim = wls[:, jj, :].max() + 0.003
        tind = (model['wl']>lolim) * (model['wl'] < hilim)
        lam_template = model['wl'][tind]
        template = model['flux'][tind]
        template /= np.median(template)

        # fit continuum coefficients / flux scaling
        ccoef = [-0.1, 1.2/np.median(template)]
        NPC = len(ccoef)
        ind90 = np.sort(fobs0[jj])[int(0.9*npix)]  
        ccoef = np.polyfit(pix[fobs0[jj]>ind90], fobs0[jj][fobs0[jj]>ind90], NPC-1)
        
        print(jj)
        for obs in np.arange(nobs):
            wcoef = np.polyfit(pix, wls[obs, jj, :], NPW-1)

            guess = np.concatenate(([21, 0.3, 9e-5], wcoef, ccoef))
            #mygmod, mygw = fit_atmo.modelspec_tel_template(guess, lam_template, template, lam_atmo, atmo, NPW, NPC, npix, retlam=True)

            thingy = mf.modelspec_template(guess, lam_template, template, NPW, NPC, npix)
            chipguesses[obs,jj] = thingy
            #sys.exit()
            fitargs = (mf.modelspec_template, lam_template, template, NPW, NPC, npix, fobs0[jj], wfobs0[jj])

            fit = mf.fmin(mf.errfunc, guess, args=fitargs, full_output=True, disp=True, maxiter=10000, maxfun=100000)
            print(fit[0][0:4])
            mymod, myw = mf.modelspec_template(fit[0], *fitargs[1:-2], retlam=True)
            #mycor = mf.modelspec_tel_template(fit[0], lam_template, np.ones(template.size), *fitargs[3:-2], retlam=False)
            chipfits.append(fit)
            chipmods[obs,jj] = mymod
            chiplams[obs,jj] = myw

            # save best parameters
            orderval.append(jj)
            obsval.append(obs)
            vsini.append(fit[0][0])
            limbdark.append(fit[0][1])
            rv.append(fit[0][2])
            chisq.append(fit[1])
            
            # make non-broadened model
            fit[0][0:2] = 0
            mymodnobroad, mywnobroad = mf.modelspec_template(fit[0], *fitargs[1:-2], retlam=True)
            chipmodnobroad[obs,jj] = mymodnobroad

            chisqarr[obs,jj] = fit[1]


    # make table of best parameters
    results = Table()
    results['order'] = orderval
    results['obs'] = obsval
    results['chisq'] = chisq
    results['vsini'] = vsini
    results['limbdark'] = limbdark
    results['rv'] = rv

    results.write(f'result/IGRINS_W1049B_fitting_results_{modelname[:12]}.txt', format='ascii')

    #fits.writeto('IGRINS_W1049B_chipmodnobroad.fits', chipmodnobroad)
    fits.writeto(f'result/IGRINS_W1049B_chipmods_{modelname[:12]}.fits', chipmods)
    fits.writeto(f'result/IGRINS_W1049B_chiplams_{modelname[:12]}.fits', chiplams)

def fit_nonstacked(modelpath):
    """
    Fit non-stacked observations, but each with specific wcoef.
    """
    ###################################
    #  Open data
    ###################################

    filelist = sorted(glob.glob(f'{homedir}/uoedrive/data/IGRINS/SDCK*_1f.spec.fits'))

    fluxes = []
    wls = []

    for filename in filelist:
        wlname = filename.split('_1f')[0]+'.wave.fits'

        flux = fits.getdata(filename)
        wl = fits.getdata(wlname)

        # trim first and last 100 columns
        flux = flux[:, 100:1948]
        wl = wl[:, 100:1948]

        #print(wl.shape)
        
        fluxes.append(flux)
        wls.append(wl)

    ###############################
    ## Collapse data over the whole observation into one mean, smoothed spectrum
    ###############################

    dims = np.array(fluxes).shape
    fluxes = np.array(fluxes)
    wls = np.array(wls)

    # remove first order with all NaNs when making various arrays
    obs0 = fluxes[:, 1:, :]*dims[0]
    eobs0 = fluxes[:, 1:, :]*np.sqrt(dims[0])  # make noise spectrum (assuming just photon noise)
    fobs0 = np.array([[signal.medfilt(obs0[obs,jj], 3) for jj in range(dims[1]-1)] for obs in range(dims[0])])
    fobs0 /= np.array([np.nanmedian(fobs0[obs], axis=1) for obs in range(dims[0])]).reshape(dims[0], dims[1]-1, 1)
    eobs0 /= np.array([np.nanmedian(fobs0[obs], axis=1) for obs in range(dims[0])]).reshape(dims[0], dims[1]-1, 1)
    wfobs0 = 1./eobs0**2   # make noise spectrum
    wind = np.isinf(wfobs0)  # remove points with infinite values
    wfobs0[wind] = 0.

    # fix wavelength array to have same shape
    wls = wls[:, 1:24, :] 

    ##########################
    # open model
    ##########################

    model = Table.read(modelpath, format='fits')
    modelname = modelpath.split("/")[-1]
    model['wl'] = model['Wavelength']
    model['flux'] = model['Flux']

    # set # of observations and # of orders to process
    nobs = wls.shape[0]
    norders = 20  # skip last few non-fitting orders for now
    npix = dims[2] 

    NPW = 4  # number of terms for polynomial fit to the wavelengths
    pix = np.arange(npix, dtype=float)/npix
    chipfits = []
    chipmods = np.zeros((nobs, norders, npix), dtype=float)
    chiplams = np.zeros((nobs, norders, npix), dtype=float)
    chipmodnobroad = np.zeros((nobs, norders, npix), dtype=float)
    chipguesses = np.zeros((nobs, norders, npix), dtype=float)

    # loop over each chip and fit the model to that chip
    chisqarr = np.zeros((nobs, norders), dtype=float)

    # set up tables for best fit vsini, limb darkening values, and rv
    orderval=[]
    obsval=[]
    vsini = []
    limbdark = []
    rv = []
    chisq = []

    for jj in (np.arange(norders)):
        lolim = wls[:, jj, :].min() - 0.003
        hilim = wls[:, jj, :].max() + 0.003
        tind = (model['wl']>lolim) * (model['wl'] < hilim)
        lam_template = model['wl'][tind]
        template = model['flux'][tind]
        template /= np.median(template)
        
        print(jj)
        for obs in np.arange(nobs):
            wcoef = np.polyfit(pix, wls[obs, jj, :], NPW-1)
            # fit continuum coefficients / flux scaling
            ccoef = [-0.1, 1.2/np.median(template)]
            NPC = len(ccoef)
            ind90 = np.sort(fobs0[obs, jj])[int(0.9*npix)]  
            ccoef = np.polyfit(pix[fobs0[obs,jj]>ind90], fobs0[obs,jj][fobs0[obs,jj]>ind90], NPC-1)

            guess = np.concatenate(([21, 0.3, 9e-5], wcoef, ccoef))
            #mygmod, mygw = fit_atmo.modelspec_tel_template(guess, lam_template, template, lam_atmo, atmo, NPW, NPC, npix, retlam=True)

            thingy = mf.modelspec_template(guess, lam_template, template, NPW, NPC, npix)
            chipguesses[obs,jj] = thingy
            #sys.exit()
            fitargs = (mf.modelspec_template, lam_template, template, NPW, NPC, npix, fobs0[obs,jj], wfobs0[obs,jj])

            fit = mf.fmin(mf.errfunc, guess, args=fitargs, full_output=True, disp=True, maxiter=10000, maxfun=100000)
            print(fit[0][0:4])
            mymod, myw = mf.modelspec_template(fit[0], *fitargs[1:-2], retlam=True)
            #mycor = mf.modelspec_tel_template(fit[0], lam_template, np.ones(template.size), *fitargs[3:-2], retlam=False)
            chipfits.append(fit)
            chipmods[obs,jj] = mymod
            chiplams[obs,jj] = myw

            # save best parameters
            orderval.append(jj)
            obsval.append(obs)
            vsini.append(fit[0][0])
            limbdark.append(fit[0][1])
            rv.append(fit[0][2])
            chisq.append(fit[1])

            chisqarr[obs,jj] = fit[1]

    # make table of best parameters
    results = Table()
    results['order'] = orderval
    results['obs'] = obsval
    results['chisq'] = chisq
    results['vsini'] = vsini
    results['limbdark'] = limbdark
    results['rv'] = rv

    results.write(f'result/IGRINS_W1049B_K_nonstack_fitting_results_{modelname[:12]}.txt', format='ascii')

    #fits.writeto('IGRINS_W1049B_chipmodnobroad.fits', chipmodnobroad)
    fits.writeto(f'result/IGRINS_W1049B_K_nonstack_chipmods_{modelname[:12]}.fits', chipmods)
    fits.writeto(f'result/IGRINS_W1049B_K_nonstack_chiplams_{modelname[:12]}.fits', chiplams)


if __name__ == "__main__":
    modellist = sorted(glob.glob(f'{homedir}/uoedrive/data/BTSettlTemps5.0/*.fits'))
    for model in modellist:
        print("running model fit", model)
        fit_nonstacked(modelpath=model)