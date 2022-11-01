import numpy as np
from astropy.table import Table, vstack
import modelfitting as mf
import sys
from scipy import signal
import pickle
import glob
from astropy.io import fits
import matplotlib.patches as patches
import matplotlib.pyplot as plt

from PyAstronomy import pyasl

###################################
#  Open data
###################################

filelist = sorted(glob.glob('/Users/bbiller/Data/Doppler_imaging_code/IGRINS_data/W1049B_Kband/20200211/SDCH*_1f.spec.fits'))

fluxes = []
wls = []

for filename in filelist:
    wlname = filename.split('_1f')[0]+'.wave.fits'

    flux = fits.getdata(filename)
    wl = fits.getdata(wlname)

    # trim first and last 100 columns
    flux = flux[:, 100:1948]
    wl = wl[:, 100:1948]

    print(wl.shape)
    
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

#model = Table.read('/Users/bbiller/Documents/proposals/CRIRES_P108_ETC/models_1615994358/bt-settl/bt-settl-cifist/lte014.0-5.0-0.0a+0.0.BT-Settl.spec.7.dat.txt', format='ascii', names=('wl', 'flux'))
#model['wl'] /= 10000.  # convert to microns

"""
# original model from Ian's daofind directory
#model = Table.read('/Users/bbiller/Data/Doppler_imaging_code/daofiles/lte015-5.0-0.0a+0.0.BT-Settl.spec.7.xz', format='fits')
model = Table.read('/Users/bbiller/Data/Doppler_imaging_code/BT-Settl_lte015-5.0-0.0+0.0_orig.fits', format='fits')
model['wl'] = model['Wavelength']
model['flux'] = model['Flux']
"""


model = Table.read('/Users/bbiller/Documents/proposals/CRIRES_P108_ETC/models_1615994358/bt-settl/bt-settl-cifist/lte015.0-5.0-0.0a+0.0.BT-Settl.spec.7.dat.txt', format='ascii', names=('wl', 'flux'))

ind = np.where((model['wl'] > 1.48*10000.) & (model['wl'] < 1.81*10000.))

ind = np.where((model['wl'] > 1.48*10000.) & (model['wl'] < 1.6*10000.))

ind = np.where((model['wl'] > 1.48*10000.) & (model['wl'] < 1.638*10000.))  # works up to 1.638 um
modelblue = model[ind]

ind = np.where((model['wl'] > 1.63834*10000.) & (model['wl'] < 1.81*10000.)) # works between 1.64 to 1.81 um

#ind = np.where((model['wl'] >= 1.6*10000.) & (model['wl'] < 1.7*10000.))

modelred = model[ind]

# convert to vacuum
modelred['wl'] = pyasl.airtovac2(modelred['wl'], mode="edlen53", maxiter=50, precision=1.e-14)
modelblue['wl'] = pyasl.airtovac2(modelblue['wl'], mode="edlen53", maxiter=50, precision=1.e-14)

model = vstack([modelblue, modelred])
model['wl'] /= 10000.  # convert to microns



NPW = 4  # number of terms for polynomial fit to the wavelengths
npix = dims[2]   # number of pixels across the entire spectrum 
pix = np.arange(npix, dtype=float)/npix
chipfits = []
chipmods = np.zeros((dims[1]-1, dims[2]), dtype=float)
chiplams = np.zeros((dims[1]-1, dims[2]), dtype=float)
chipmodnobroad = np.zeros((dims[1]-1, dims[2]), dtype=float)
chipguesses = np.zeros((dims[1]-1, dims[2]), dtype=float)

# loop over each chip and fit the model to that chip
chisqarr = np.zeros((dims[1]-1), dtype=float)

for jj in (np.arange(dims[1]-1))[0:23]:
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
    chipguesses[jj] = thingy
    #sys.exit()
    fitargs = (mf.modelspec_template, lam_template, template, NPW, NPC, npix, fobs0[jj], wfobs0[jj])

    fit = mf.fmin(mf.errfunc, guess, args=fitargs, full_output=True, disp=True, maxiter=10000, maxfun=100000)
    print(fit[0][0:4])
    mymod, myw = mf.modelspec_template(fit[0], *fitargs[1:-2], retlam=True)
    #mycor = mf.modelspec_tel_template(fit[0], lam_template, np.ones(template.size), *fitargs[3:-2], retlam=False)
    chipfits.append(fit)
    chipmods[jj] = mymod
    chiplams[jj] = myw

    # make non-broadened model
    fit[0][0:2] = 0
    mymodnobroad, mywnobroad = mf.modelspec_template(fit[0], *fitargs[1:-2], retlam=True)
    chipmodnobroad[jj] = mymodnobroad

    chisqarr[jj] = fit[1]
    # make non-broadened model

    #mymodnobroad, mywnobroad = mf.modelspec_tel_template((np.concatenate((np.array([0, 0]), fit[0][2:]))), *fitargs[1:-2], retlam=True)
    #chipmodnobroad[jj] = mymodnobroad

wl = wl[1:, :]

pltstart = np.arange(5)*5

fignames = ['part1', 'part2', 'part3', 'part4', 'part5']

for i in np.arange(5):
    fig = plt.figure(figsize=(8, 25))
    #ax = fig.add_subplot(111)

    for jj in np.arange(5):
        ax = fig.add_subplot(5,1,jj+1)
        ax.plot(wl[jj+pltstart[i], :], fobs0[jj+pltstart[i], :], color='black')
        ax.plot(wl[jj+pltstart[i], :], chipguesses[jj+pltstart[i], :], color='green')
        ax.plot(wl[jj+pltstart[i], :], chipmods[jj+pltstart[i], :], color='salmon')
        #ax.plot(wl[jj+pltstart[i], :], chipmodnobroad[jj+pltstart[i], :], color='magenta')
        #ax.text(np.min(wl[jj+pltstart[i], :]), 1, chisqarr[jj+pltstart[i]])
        ax.text(0.1, 0.1, chisqarr[jj+pltstart[i]], horizontalalignment='center',
                    verticalalignment='center', transform=ax.transAxes)

       
    ax.set_xlabel('Wavelength ($\mu$m)')

    fig.savefig('bestfits_'+fignames[i]+'_Hband_IGRINS_origmodnotelluric_chisq_panelled.png')


sys.exit()    

fig = plt.figure()
#fig = plt.figure()
#ax = fig.add_subplot(111)

nplots = 4
firstorder = 6
#for jj in (np.arange(dims[1]))[1:]:
for jj in (np.arange(nplots) + firstorder):
#    ax = fig.add_subplot(dims[1],1,jj+1)
    ax = fig.add_subplot(nplots, 1, jj-(firstorder-1))
    
    ax.plot(wl[jj, :], fobs0[jj, :], color='black')

    #ax.plot(chiplams[jj, :], chipmods[jj, :], color='cyan')

#    ax.plot(chiplams[jj, :], chipmodnobroad[jj, :], color='magenta')

    #ax.plot(wl[jj, :], chipmods[jj, :], color='salmon')

    ax.plot(wl[jj, :], chipmodnobroad[jj, :], color='chartreuse')
#    ax.text(np.min(wl[jj], :]), 1, chisqarr[jj])



fig = plt.figure()
    

ind = np.where((chisqarr < 400.) & (chisqarr > 0.) & (np.isfinite(chisqarr)))# & (chisqarr > 1000.))

ngood = len(ind[0])

for jj in (np.arange(ngood)):
    ax = fig.add_subplot(ngood,1,jj+1)

    ax.plot(wl[ind[0][jj], :], fobs0[ind[0][jj], :], color='black')

    #ax.plot(chiplams[jj, :], chipmods[jj, :], color='cyan')

#    ax.plot(chiplams[jj, :], chipmodnobroad[jj, :], color='magenta')

    ax.plot(wl[ind[0][jj], :], chipmods[ind[0][jj], :], color='salmon')
    

    ax.text(np.min(wl[ind[0][jj], :]), 0.8, chisqarr[ind][jj])
    
sys.exit()
    
    
sys.exit()
    
fig = plt.figure()
ax1 = fig.add_subplot(111)

for jj in (np.arange(dims[1]))[1:]:    
    ax1.plot(chiplams[jj, :], chipmodnobroad[jj, :])
    ax1.text(np.median(chiplams[jj, :]), np.max(chipmodnobroad[jj, :]), jj)
    
#allfits.append(chipfits)
#allmods.append(chipmods)
#alllams.append(chiplams)
# load crires dataset, overplot CRIRES chips
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

    

#linelist = Table.read('/Users/bbiller/Data/Doppler_imaging_code/daofiles/lte015.5-5.5-0.0a+0.0.BT-Settl.spec.7_dao_edited.clineslsd', format='ascii')

#for xc in linelist['col1']:
#    ax1.axvline(x=xc/10000.)

    
