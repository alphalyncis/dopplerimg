import numpy as np
from astropy.table import Table
import modelfitting as mf
from scipy import signal
import pickle
import matplotlib.pyplot as plt

# open data

filename = '../data/fainterspectral-fits_6.pickle'
f = open(filename, 'rb')
ret = pickle.load(f, encoding="latin1")

obs0 = np.median(ret['obs0'], axis=0)*14.
eobs0 = np.median(ret['obs0'], axis=0)*np.sqrt(14.) #TODO: why is obs0 not sqrted
fobs0 = np.vstack([signal.medfilt(obs0[jj], 3) for jj in range(4)])
eobs0 = eobs0 / np.median(fobs0, 1).reshape(4,1)
fobs0 = fobs0 / np.median(fobs0, 1).reshape(4,1)
wfobs0 = 1./eobs0**2
wind = np.isinf(wfobs0)
wfobs0[wind] = 0.

# open model

model = Table.read('../data/BT-Settl_lte015-5.0-0.0+0.0_orig.fits', format='fits')
model['wl'] = model['Wavelength']
model['flux'] = model['Flux']

NPW = 4
npix = ret['wobs'].shape[1]
pix = np.arange(npix, dtype=float)/npix
chipfits = []
chipmods = np.zeros((4, 1024), dtype=float)
chiplams = np.zeros((4, 1024), dtype=float)
chipmodnobroad = np.zeros((4, 1024), dtype=float)
  
for jj in np.arange(4): # TODO: separate fit for each spectra, but physically some parameter should be same among spectra?
    lolim = ret['wobs'][jj].min() - 0.003  #TODO: where this value from
    hilim = ret['wobs'][jj].max() + 0.003
    tind = (model['wl']>lolim) * (model['wl'] < hilim)
    lam_template = model['wl'][tind]
    template = model['flux'][tind]
    template /= np.median(template)

    wcoef = np.polyfit(pix, ret['wobs'][jj], NPW-1)
    ccoef = [-0.1, 1.2/np.median(template)]
    NPC = len(ccoef)
    ind90 = np.sort(fobs0[jj])[int(0.9*npix)]  
    ccoef = np.polyfit(pix[fobs0[jj]>ind90], fobs0[jj][fobs0[jj]>ind90], NPC-1)

    guess = np.concatenate(([21, 0.3, 9e-5], wcoef, ccoef))
    #mygmod, mygw = fit_atmo.modelspec_tel_template(guess, lam_template, template, lam_atmo, atmo, NPW, NPC, npix, retlam=True)

    #thingy = mf.modelspec_template(guess, lam_template, template, NPW, NPC, npix)
    #sys.exit()
    fitargs = (mf.modelspec_template, lam_template, template, NPW, NPC, npix, fobs0[jj], wfobs0[jj])
    fit = mf.fmin(mf.errfunc, guess, args=fitargs, full_output=True, disp=True, maxiter=10000, maxfun=10000)
    mymod, myw = mf.modelspec_template(fit[0], *fitargs[1:-2], retlam=True)
    #mycor = mf.modelspec_tel_template(fit[0], lam_template, np.ones(template.size), *fitargs[3:-2], retlam=False)
    chipfits.append(fit)
    chipmods[jj] = mymod
    chiplams[jj] = myw

    # make non-broadened model
    fit[0][0:2] = 0
    mymodnobroad, mywnobroad = mf.modelspec_template(fit[0], *fitargs[1:-2], retlam=True)
    chipmodnobroad[jj] = mymodnobroad
    
#allfits.append(chipfits)
#allmods.append(chipmods)
#alllams.append(chiplams)

fig = plt.figure(figsize=(15,11))
    
for jj in np.arange(4):
    ax = fig.add_subplot(4,1,jj+1)
    ax.plot(ret['wobs'][jj, :], fobs0[jj, :], color='black', label="observation")
    ax.plot(ret['wobs'][jj, :], chipmods[jj, :], color='red', label="best-fit")
    ax.plot(ret['wobs'][jj, :], chipmodnobroad[jj, :], color='red', alpha=0.1, label="best-fit but no broaden")
    ax.plot(ret['wobs'][jj, :], ret['chipmods'][0, jj, :], color='green', alpha=0.4, linewidth=1.5, label="???")

plt.legend()

    
  
