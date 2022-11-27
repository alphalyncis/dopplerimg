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

#open data

filelist = sorted(glob.glob('/Users/bbiller/Data/Doppler_imaging_code/IGRINS_data/2M0036/retell/SDCK*.spec_a0v.fits'))

fluxes = []
wls = []

for filename in filelist:
    hdu = fits.open(filename)

    flux = hdu[0].data
    wl = hdu[1].data

    # trim first and last 100 columns
    flux = flux[:, 100:1948]
    wl = wl[:, 100:1948]

    print(wl.shape)
    
    fluxes.append(flux)
    wls.append(wl)

    
    
dims = np.array(fluxes).shape
fluxes = np.array(fluxes)
wls = np.array(wls)

obs0 = np.median(fluxes, axis=0)*dims[0]
eobs0 = np.median(fluxes, axis=0)*np.sqrt(dims[0])
fobs0 = np.vstack([signal.medfilt(obs0[jj], 3) for jj in range(dims[1])])
eobs0 /= np.nanmedian(fobs0, 1).reshape(dims[1],1)
fobs0 /= np.nanmedian(fobs0, 1).reshape(dims[1],1)
wfobs0 = 1./eobs0**2
wind = np.isinf(wfobs0)
wfobs0[wind] = 0.

# open model
model = Table.read('/Users/bbiller/Documents/proposals/CRIRES_P108_ETC/models_1615994358/bt-settl/bt-settl-cifist/lte014.0-5.0-0.0a+0.0.BT-Settl.spec.7.dat.txt', format='ascii', names=('wl', 'flux'))
model['wl'] /= 10000.

NPW = 4
npix = dims[2]
pix = np.arange(npix, dtype=float)/npix
chipfits = []
chipmods = np.zeros((dims[1], dims[2]), dtype=float)
chiplams = np.zeros((dims[1], dims[2]), dtype=float)
chipmodnobroad = np.zeros((dims[1], dims[2]), dtype=float)
  
for jj in (np.arange(dims[1]))[1:23]:
    print(jj)
    lolim = wls[:, jj, :].min() - 0.003
    hilim = wls[:, jj, :].max() + 0.003
    tind = (model['wl']>lolim) * (model['wl'] < hilim)
    lam_template = model['wl'][tind]
    template = model['flux'][tind]
    template /= np.median(template)

    wcoef = np.polyfit(pix, wls[0, jj, :], NPW-1)
    ccoef = [-0.1, 1.2/np.median(template)]
    NPC = len(ccoef)
    ind90 = np.sort(fobs0[jj])[int(0.9*npix)]  
    ccoef = np.polyfit(pix[fobs0[jj]>ind90], fobs0[jj][fobs0[jj]>ind90], NPC-1)

    guess = np.concatenate(([21, 0.3, 9e-5], wcoef, ccoef))
    #mygmod, mygw = fit_atmo.modelspec_tel_template(guess, lam_template, template, lam_atmo, atmo, NPW, NPC, npix, retlam=True)

    thingy = mf.modelspec_template(guess, lam_template, template, NPW, NPC, npix)
    #sys.exit()
    fitargs = (mf.modelspec_template, lam_template, template, NPW, NPC, npix, fobs0[jj], wfobs0[jj])
    fit = mf.fmin(mf.errfunc, guess, args=fitargs, full_output=True, disp=False)
    mymod, myw = mf.modelspec_template(fit[0], *fitargs[1:-2], retlam=True)
    #mycor = mf.modelspec_tel_template(fit[0], lam_template, np.ones(template.size), *fitargs[3:-2], retlam=False)
    chipfits.append(fit)
    chipmods[jj] = mymod
    chiplams[jj] = myw

    # make non-broadened model
    fit[0][0:2] = 0
    mymodnobroad, mywnobroad = mf.modelspec_template(fit[0], *fitargs[1:-2], retlam=True)
    chipmodnobroad[jj] = mymodnobroad



fig = plt.figure()
#fig = plt.figure()
#ax = fig.add_subplot(111)
    
#for jj in (np.arange(dims[1]))[1:]:
for jj in (np.arange(10) + 10):
#    ax = fig.add_subplot(dims[1],1,jj+1)
    ax = fig.add_subplot(10, 1, jj-9)
    
    ax.plot(wl[jj, :], fobs0[jj, :], color='black')

#    ax.plot(chiplams[jj, :], chipmods[jj, :], color='cyan')

#    ax.plot(chiplams[jj, :], chipmodnobroad[jj, :], color='magenta')

    ax.plot(wl[jj, :], chipmods[jj, :], color='salmon')

    ax.plot(wl[jj, :], chipmodnobroad[jj, :], color='chartreuse')

    
sys.exit()

    

fig = plt.figure()
ax1 = fig.add_subplot(111)

for jj in (np.arange(dims[1]))[1:23]:    
    ax1.plot(chiplams[jj, :], chipmodnobroad[jj, :])
    
#allfits.append(chipfits)
#allmods.append(chipmods)
#alllams.append(chiplams)
# load crires dataset
datdir='/Users/bbiller/Data/Doppler_imaging_code/'
filename = datdir+'fainterspectral-fits_6.pickle'
f = open(filename, 'rb')
ret = pickle.load(f, encoding="latin1")

# get min and max wavelengths for CRIRES chips
minvals = []
maxvals = []

for i in np.arange(4):
    minvals.append(np.min(ret['chiplams'][:, i, :]))
    maxvals.append(np.max(ret['chiplams'][:, i, :]))


for l1, l2 in zip(minvals, maxvals):
    ax1.add_patch(
    patches.Rectangle(
        (l1, 0),   # (x,y)
        (l2 - l1),          # width
        2.0,          # height
        color='yellow', alpha=0.5))#, label=lname )

    
ax1.plot(chiplams[1:, :], fobs0[1:, :], 'k.')
