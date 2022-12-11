from __future__ import print_function, division

import numpy as np
from astropy.table import Table
from astropy.io import fits
import modelfitting as mf
import sys
from scipy import signal
import pickle
import matplotlib.patches as patches
import matplotlib.pyplot as plt

# from PyAstronomy import pyasl

#open data

filename = '../data/fainterspectral-fits_6.pickle'
f = open(filename, 'rb')
ret = pickle.load(f, encoding="latin1")


obs0 = np.median(ret['obs0'], axis=0)*14.
eobs0 = np.median(ret['obs0'], axis=0)*np.sqrt(14.)
fobs0 = np.vstack([signal.medfilt(obs0[jj], 3) for jj in range(4)])
eobs0 /= np.median(fobs0, 1).reshape(4,1)
fobs0 /= np.median(fobs0, 1).reshape(4,1)
wfobs0 = 1./eobs0**2
wind = np.isinf(wfobs0)
wfobs0[wind] = 0.


# open model
"""
#CIFIST BT-Settl models
model = Table.read('/Users/bbiller/Documents/proposals/CRIRES_P108_ETC/models_1615994358/bt-settl/bt-settl-cifist/lte015.0-5.0-0.0a+0.0.BT-Settl.spec.7.dat.txt', format='ascii', names=('wl', 'flux'))

ind = np.where((model['wl'] > 22500.) & (model['wl'] < 24000.))
model = model[ind]

# correct back to vacuum wavelengths
model['wl'] = pyasl.airtovac2(model['wl'], mode="edlen53")

# correct to same flux units as original model
model['flux'] *= 10.
"""

"""
# original BT-Settl models
#model = Table.read('/Users/bbiller/Data/Doppler_imaging_code/btsettl_orig/models_1633965013/bt-settl/lte014-5.0-0.0.BT-Settl.7.dat.txt', format='ascii', names=('wl', 'flux'))

# AGSS 2009 BT-Settl models
#model = Table.read('/Users/bbiller/Data/Doppler_imaging_code/btsettl_AGSS2009/models_1633965947/bt-settl-agss/lte014-5.0-0.0.BT-Settl.7.dat.txt', format='ascii', names=('wl', 'flux'))
"""

# convert to um
#model['wl'] /= 10000.
#model['Wavelength'] = model['wl']/10000.
#model['Flux'] = model['flux']

# original model from Ian's daofind directory
#model = Table.read('/Users/bbiller/Data/Doppler_imaging_code/daofiles/lte015-5.0-0.0a+0.0.BT-Settl.spec.7.xz', format='fits')
model = Table.read('../data/BT-Settl_lte015-5.0-0.0+0.0_orig.fits', format='fits')

model['wl'] = model['Wavelength']
model['flux'] = model['Flux']

# open telluric spectrum
telfn = '/Users/bbiller/Data/Doppler_imaging_code/transdata_0,5-14_mic_hires.fits'
#transdata_2,28-2,37_mic_hires.fits'
atm0 = fits.getdata(telfn)

NPW = 4
npix = ret['wobs'].shape[1]
pix = np.arange(npix, dtype=float)/npix
chipfits = []
chipmods = np.zeros((4, 1024), dtype=float)
chiplams = np.zeros((4, 1024), dtype=float)
chipmodnobroad = np.zeros((4, 1024), dtype=float)
  
for jj in np.arange(4):
    lolim = ret['wobs'][jj].min() - 0.003
    hilim = ret['wobs'][jj].max() + 0.003
    tind = (model['wl']>lolim) * (model['wl'] < hilim)
    aind = (atm0[:, 0]>lolim) * (atm0[:, 0] < hilim)
    
    lam_template = model['wl'][tind]
    template = model['flux'][tind]
    template /= np.median(template)

    lam_atmo = atm0[aind, 0]
    atmo = atm0[aind, 1]
    
    wcoef = np.polyfit(pix, ret['wobs'][jj], NPW-1)
    ccoef = [-0.1, 1.2/np.median(template)]
    NPC = len(ccoef)
    ind90 = np.sort(fobs0[jj])[int(0.9*npix)]  
    ccoef = np.polyfit(pix[fobs0[jj]>ind90], fobs0[jj][fobs0[jj]>ind90], NPC-1)

    guess = np.concatenate(([21, 0.3, 9e-5, 42, 1.3], wcoef, ccoef))
    #mygmod, mygw = fit_atmo.modelspec_tel_template(guess, lam_template, template, lam_atmo, atmo, NPW, NPC, npix, retlam=True)

    #thingy = mf.modelspec_template(guess, lam_template, template, NPW, NPC, npix)
    #sys.exit()
#    fitargs = (mf.modelspec_template, lam_template, template, NPW, NPC, npix, fobs0[jj], wfobs0[jj])  # No telluric correction

    fitargs = (mf.modelspec_tel_template, lam_template, template, lam_atmo, atmo, NPW, NPC, npix, fobs0[jj], wfobs0[jj])  # with telluric correction
    
    fit = mf.fmin(mf.errfunc, guess, args=fitargs, full_output=True, disp=True, maxiter=10000, maxfun=10000)
    #mymod, myw = mf.modelspec_template(fit[0], *fitargs[1:-2], retlam=True)
    mymod, myw = mf.modelspec_tel_template(fit[0], *fitargs[1:-2], retlam=True)
    mycor = mf.modelspec_tel_template(fit[0], lam_template, np.ones(template.size), *fitargs[3:-2], retlam=False)
    print(fit[0])
    chipfits.append(fit)
    chipmods[jj] = mymod
    chiplams[jj] = myw

    # make non-broadened model
    
    """
    nonbroadfit = fit[0]
    nonbroadfit[0:2] = 0.

    fitnobroad = fit[0].copy()
    fitnobroad[0] = 0.
    mymodnobroad = mf.modelspec_tel_template(fitnobroad, *fitargs[1:-2])
    
    #mymodnobroad, mywnobroad = mf.modelspec_tel_template(fit, *fitargs[1:-2], retlam=True)
    chipmodnobroad[jj] = mymodnobroad
    """

    # make non-broadened model
    fit[0][0:2] = 0
    mymodnobroad, mywnobroad = mf.modelspec_tel_template(fit[0], *fitargs[1:-2], retlam=True)
    chipmodnobroad[jj] = mymodnobroad
    

#sys.exit()
    
#allfits.append(chipfits)
#allmods.append(chipmods)
#alllams.append(chiplams)

fig = plt.figure()
#fig = plt.figure()
#ax = fig.add_subplot(111)
    
for jj in np.arange(4):
    ax = fig.add_subplot(4,1,jj+1)

    ax.plot(ret['wobs'][jj, :], fobs0[jj, :], color='black', label='Data')

#    ax.plot(chiplams[jj, :], chipmods[jj, :], color='cyan')

#    ax.plot(chiplams[jj, :], chipmodnobroad[jj, :], color='magenta')

    ax.plot(ret['wobs'][jj, :], chipmods[jj, :], color='salmon', label='Best fit old model')

    #ax.plot(ret['wobs'][jj, :], chipmodnobroad[jj, :], color='chartreuse', label='Best fit non-broadened CIFIST model')

    #ax.plot(ret['wobs'][jj, :], ret['chipmodnobroad'][0, jj, :], color='blue', label='original non-broadened model fit')

    ax.plot(ret['wobs'][jj, :], ret['chipmods'][0, jj, :], color='purple', label='original broadened model fit')

    ax.set_xlabel('Wavelength ($\mu$m)')

ax.legend()
    
fig.savefig('chipmod_CRIRES_origmodel_nonbroadened.pdf')

    

    
  
