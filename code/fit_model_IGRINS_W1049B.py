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


# original model from Ian's daofind directory
#model = Table.read('/Users/bbiller/Data/Doppler_imaging_code/daofiles/lte015-5.0-0.0a+0.0.BT-Settl.spec.7.xz', format='fits')
model = Table.read('/Users/bbiller/Data/Doppler_imaging_code/BT-Settl_lte015-5.0-0.0+0.0_orig.fits', format='fits')
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
        # make non-broadened model

        #mymodnobroad, mywnobroad = mf.modelspec_tel_template((np.concatenate((np.array([0, 0]), fit[0][2:]))), *fitargs[1:-2], retlam=True)
        #chipmodnobroad[obs,jj] = mymodnobroad

# make table of best parameters
results = Table()
results['order'] = orderval
results['obs'] = obsval
results['chisq'] = chisq
results['vsini'] = vsini
results['limbdark'] = limbdark
results['rv'] = rv

results.write('IGRINS_W1049B_fitting_results.txt', format='ascii')

fits.writeto('IGRINS_W1049B_chipmodnobroad.fits', chipmodnobroad)
fits.writeto('IGRINS_W1049B_chipmods.fits', chipmods)
fits.writeto('IGRINS_W1049B_chiplams.fits', chiplams)


sys.exit()
        
wl = wl[1:, :]

pltstart = np.arange(5)*5

fignames = ['part1', 'part2', 'part3', 'part4', 'part5']

for i in np.arange(5):
    fig = plt.figure(figsize=(8, 25))
    #ax = fig.add_subplot(111)

    for jj in np.arange(5):
        ax = fig.add_subplot(5,1,jj+1)
        ax.plot(wl[jj+pltstart[i], :], fobs0[jj+pltstart[i], :], color='black')
        #ax.plot(wl[jj+pltstart[i], :], chipguesses[jj+pltstart[i], :], color='green')
        #ax.plot(wl[jj+pltstart[i], :], chipmods[jj+pltstart[i], :], color='salmon')
        ax.plot(wl[jj+pltstart[i], :], chipmodnobroad[jj+pltstart[i], :], color='magenta')
        #ax.text(np.min(wl[jj+pltstart[i], :]), 1, chisqarr[jj+pltstart[i]])
        ax.text(0.1, 0.1, chisqarr[jj+pltstart[i]], horizontalalignment='center',
                    verticalalignment='center', transform=ax.transAxes)

       
    ax.set_xlabel('Wavelength ($\mu$m)')

    fig.savefig('bestfits_'+fignames[i]+'_IGRINS_origmodnotelluric_nonbroadonly_chisq_panelled.png')


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

for jj in (np.arange(6)+1):
    ax1.plot(chiplams[jj, :], fobs0[jj, :], color='black')
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

ax1.set_xlabel('Wavelength ($\mu$m)', fontsize=25)
ax1.set_ylabel('Normalized Flux', fontsize=25)

pmod = os.path.split(modelfn)[1].replace('BT-Settl.spec.7.fits', '_dao').replace('.',',')
f_linelist   = _mod + pmod + '_edited.clineslsd'
#f_linelist = _mod + 'btsettl_CRIRES_daospec.clines'  # using new linelist for now
#f_linelist = _mod + 'btsettl_CRIRES_daospec_refitted.clines'  # using new linelist for now
f_pspec      = _mod + pmod + '.fits'
f_pspec_cont = _mod + pmod + 'C.fits'

tab = Table.read(f_linelist, format='ascii')

for line in tab['col1']:
    ax1.axvline(x=(line*(1.+9e-5))/10000., color='black')



#linelist = Table.read('/Users/bbiller/Data/Doppler_imaging_code/daofiles/lte015.5-5.5-0.0a+0.0.BT-Settl.spec.7_dao_edited.clineslsd', format='ascii')

#for xc in linelist['col1']:
#    ax1.axvline(x=xc/10000.)

    
