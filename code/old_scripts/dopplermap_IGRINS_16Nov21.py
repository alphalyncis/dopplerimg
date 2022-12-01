"""
######################################################

#08-12-20
#Changed the input to the 'mega chip' - all pixels on one row so no NaNs (was previously using a mask to change wavelengths outside of CRIRES range to NaNs)

#17-11-2020
#begin integrations IGRINS capability

#10-11-2020
#pipeline works well with crires data set
#made a copy to keep - EBubb_dopplermap_CRIRES_copy

#29-05-2020
#finished converting to python3
#now tidy and integrate equal area coordinates to finalise


# 18-02-2020
# Emma Bubb - change to run on Python 3

# 29-06-2021
# removed telluric correction to make IGRINS-only version

######################################################

# Original author: Ian Crossfield (Python 2.7)
##################################################
#2013-08-08 16:13 IJMC: 
##################################################

######################################################
"""

# EB: import IC toolkits
import tools3 as tools
import dia3 as dia
import dao3 as dao
import nsdata3 as ns
import analysis3 as an
import maps3 as maps
import EBubb_DI_code_EB_map_class as EBmap
import ELL_map_class as ELL_map
import spec3 as spec
import dime3 as dime # Doppler Imaging & Maximum Entropy 
import phasecurves3 as pc
import modelfitting as mf

#EB: import packages
import numpy as np
import os
from scipy import interpolate
from scipy import signal
from scipy import optimize
from astropy.io import fits
from astropy.table import Table
import glob

import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

import sys

######################################################################################################################################################

"""
#my locations
"""

datdir = '/Users/bbiller/Data/Doppler_imaging_code/doppler_imaging-master-EBubb_DI_code/EBubb_DI_code/'

#set file name of equal area coords_lat_20_40
coords_filename = datdir+'coords_lat_20_40.npy'
nlat, nlon = 20, 40
#nlat, nlon = 40, 80

nobs = 56
    
"""
object rotation parameters
"""

#W1049B
#set period and number of observations
period = 4.87 # period in hours
per = period*3600 # period in seconds
vsini = 29e3 # m/s
inc = 0.3491 # rad

# Specify user-defined options:

LLD = 1.0
alpha = 4500

nk = 103 # get rid of ridiculous loop over nks
#nk = 63
ftol = 0.01 # tolerance for convergence of maximum-entropy

spotfit_1 = False #True #EB: change spotfit to true
spotfit_2 = False #True #True
nstep = 1500

##################################################
# Begin!
##################################################

# load data and get timestamps

filelist = sorted(glob.glob('/Users/bbiller/Data/Doppler_imaging_code/IGRINS_data/W1049B_Kband/20200211/SDCK*_1f.spec.fits'))

fluxes = []
wls = []
mjd = []

for filename in filelist:
    wlname = filename.split('_1f')[0]+'.wave.fits'

    hdu = fits.open(filename.split('_1f')[0]+'.spec.fits')
    flux = fits.getdata(filename)
    wl = fits.getdata(wlname)

    mjd0 = hdu[0].header['MJD-OBS']
    print(filename)
    print(hdu[0].header['MJD-OBS'])
    print(hdu[0].header['JD-OBS'])
    
    
    # trim first and last 100 columns
    flux = flux[:, 100:1948]
    wl = wl[:, 100:1948]

    print(wl.shape)
    
    fluxes.append(flux)
    wls.append(wl)
    mjd.append(mjd0)

    
    
dims = np.array(fluxes).shape
fluxes = np.array(fluxes)
wls = np.array(wls)
mjd = np.array(mjd)

phase = (mjd-mjd[0])*86400 * 2*np.pi/ per

# fix wavelength and flux array to have same shape and remove the first order with all NaNs
wls = wls[:, 1:24, :]
fluxes = fluxes[:, 1:24, :]
    
# load non-broadened model files
chipmodnobroad = fits.getdata('/Users/bbiller/Data/Doppler_imaging_code/IGRINS_W1049B_median_chipmodnobroad.fits')    
fitres = Table.read('/Users/bbiller/Data/Doppler_imaging_code/IGRINS_W1049B_median_fitting_results.txt', format='ascii')
chiplams = fits.getdata('/Users/bbiller/Data/Doppler_imaging_code/IGRINS_W1049B_median_wlmed.fits')

"""##################################
LEAST SQUARES DECONVOLUTION to build the composite line profile
"""#####################################    

# convert to angstroms
chiplams *= 10000.
    

# Set up for Least Squares Deconvolution
_model = './'
pmod = 'output'
npix = chiplams.shape[1]
_mod = '/Users/bbiller/Data/Doppler_imaging_code/daofiles/'
modelfn = _mod+ 'lte015-5.0-0.0a+0.0.BT-Settl.spec.7.fits'  #EB change file name of model

pmod = os.path.split(modelfn)[1].replace('BT-Settl.spec.7.fits', '_dao').replace('.',',')
f_linelist   = _mod + pmod + '_edited.clineslsd'
#f_linelist = _mod + 'btsettl_CRIRES_daospec.clines'  # using new linelist for now
#f_linelist = _mod + 'btsettl_CRIRES_daospec_refitted.clines'  # using new linelist for now
f_pspec      = _mod + pmod + '.fits'
f_pspec_cont = _mod + pmod + 'C.fits'

# Load the LSD files:
(lineloc, lineew, linespec) = dao.getlines(f_linelist)

pspec_cont = fits.getdata(f_pspec_cont) 
hdr_pspec_cont = fits.getheader(f_pspec_cont) 
wspec = hdr_pspec_cont['crval1'] + np.arange(pspec_cont.size)*hdr_pspec_cont['cdelt1']
spline = interpolate.UnivariateSpline(wspec, pspec_cont, s=0.0, k=1.0) #set up interpolation over the continuum measurement

# Compute LSD:
kerns = np.zeros((nobs, 2, nk), dtype=float)
modkerns = np.zeros((nobs, 2, nk), dtype=float)

firstorder=4
  
### In this scheme, run LSD separately for each frame's wavelength solution:
for jj in np.arange(2): # make Doppler map for orders 4 and 5
    for kk in range(nobs): # loop over all observations
#        deltaspec = ns.linespec(lineloc*(1.+fitres['rv'][jj+firstorder]), lineew, chiplams[jj], verbose=False, cont=spline(chiplams[jj]))
        # EB: create delta-function line spectrum from wavelength
#grid, list of line locations, list of equivalent widths


        deltaspec = ns.linespec(lineloc*(1.+8.e-5), lineew, chiplams[jj], verbose=False, cont=spline(chiplams[jj])) # try just a single RV shift
        m,kerns[kk,jj],b,c = dia.dsa(deltaspec, fluxes[kk,jj+firstorder], nk)  # / chipcors[kk, jj], nk)    # EB: DSA=Difference Spectral Analysis. Match given spectra to reference spectra  # removed telluric correction here
        m,modkerns[kk,jj],b,c = dia.dsa(deltaspec, chipmodnobroad[jj+firstorder], nk)# / chipcors[kk, jj], nk)    
			#EB: dia.dsa returns: Response matrix, kernel used in convolution, background offset, chisquared of fit
			#EB: only the kernels are used from here on        

# Compute LSD velocity grid:
dbeta = np.diff(chiplams).mean()/chiplams.mean() # EB: np.diff returns array containing the difference between elements of the array provided. so dbeta will be mean difference/mean value.
dx = -dbeta * np.arange(np.floor(-nk/2.+.5), np.floor(nk/2.+.5)) # EB: edit np.arange
dv = an.c*dx / 1e3 # km/s # EB: an.c is the speed of light


intrinsic_profile = modkerns.mean(0).mean(0)
intrinsic_profile -= intrinsic_profile[0] + (intrinsic_profile[-1] - intrinsic_profile[0])/(nk - 1.) * np.arange(nk)
systematic_rv_offset = (intrinsic_profile==intrinsic_profile.max()).nonzero()[0] - (dv==0).nonzero()[0]
intrinsic_profile = np.interp(np.arange(nk), np.arange(nk) - systematic_rv_offset, intrinsic_profile)
kerns = np.array([[np.interp(np.arange(nk), np.arange(nk) - systematic_rv_offset, kerns[jj, kk]) for kk in np.arange(2)] for jj in range(nobs)])


"""###################################################################################
	# Prepare the Doppler Imaging "response matrix"
"""######################################################################################
	
modIP = 1. - np.concatenate((np.zeros(300), intrinsic_profile, np.zeros(300)))
modDV = - np.arange(np.floor(-modIP.size/2.+.5), np.floor(modIP.size/2.+.5)) * dbeta * an.c / 1e3
flineSpline2 = interpolate.UnivariateSpline(modDV[::-1], modIP[::-1], k=1., s=0.) 

#sys.exit()

observation = 1. - kerns.mean(1).ravel()
observation_norm = observation.copy()
for ii in range(nobs): # EB: frame 0-13
    i0, i1 = ii*nk, (ii+1)*nk # EB: nk is the number of pixels in LSD computation on this run through
    inds = np.concatenate((np.arange(i0, i0+7), np.arange(i1-7, i1)))
    continuumfit = np.polyfit(inds, observation_norm[inds], 1)
    observation_norm[i0:i1] /= np.polyval(continuumfit, np.arange(i0, i1))


"""###################################################################################
	# Start the Imaging Analysis
"""###################################################################################

nx = nlat*nlon
#nx=490 #number of cells in the 20-40 lat-lon equal area division
dime.setup(observation_norm.size, nk)
flatguess = 100*np.ones(nx)
bounds = [(1e-6, 300)]*nx

    
allfits = []

#initialize Doppler map object

mmap = ELL_map.map(nlat=nlat, nlon=nlon, i=inc) #maps.map returns a class object 
#mmap = EBmap.map(filename= coords_filename, inc=inc, deltatheta=0) #EBmap.map returns a class object 
#np.save('mmap_cells_LLD-%i_vsini-%i_alpha-%i_inc-%i' % (LLD, vsini, alpha, inc), mmap.cells)

#np.save(savelocation+'EB_mmap_cells_LLD-%i_vsini-%i_alpha-%i_inc-%i' % (LLD, vsini, alpha, inc*180/np.pi), mmap.cells)
ncell = mmap.ncell
latlon_corners = np.array([c.corners_latlon for c in mmap.cells]) #corner coordinates
phi_corner = latlon_corners[:,1,:]
theta_corner = latlon_corners[:,0,:]

"""###################################################################################
# Compute the "Response Matrix" of Vogt 1987:
"""###################################################################################

Rmatrix = np.zeros((ncell, nobs*dv.size), dtype=np.float32) 

for kk, rot in enumerate(phase):
	speccube = np.zeros((ncell, dv.size), dtype=np.float32) 
	this_map = ELL_map.map(nlat=nlat, nlon=nlon, i=inc, deltaphi=-rot) #maps.map returns a class object 
	#this_map = EBmap.map(filename=coords_filename, inc=inc, deltatheta=-rot) #EBmap.map returns a class object 
	this_doppler = 1. + vsini*this_map.visible_rvcorners.mean(1)/an.c/np.cos(inc)
	good = (this_map.projected_area>0) * np.isfinite(this_doppler)
	for ii in good.nonzero()[0]:
		speccube[ii,:] = flineSpline2(dv + (this_doppler[ii]-1)*an.c/1000.)
		#speccube[ii,:] = flineSpline(lamdv / this_doppler[ii])
		limbdarkening = (1. - LLD) + LLD * this_map.mu
		Rblock = speccube * ((limbdarkening*this_map.projected_area).reshape(this_map.ncell, 1)*np.pi/this_map.projected_area.sum())
		Rmatrix[:,dv.size*kk:dv.size*(kk+1)] = Rblock

flatmodel = dime.normalize_model(np.dot(flatguess, Rmatrix), nk)        

if len(allfits)==0:  # Properly scale measurement weights:
    minfm = flatmodel.min()
    cutoffval = 1. - (1. - minfm) / 22.
    w_observation = (flatmodel < cutoffval).astype(float) / np.median(kerns.mean(1).std(0))**2
        # Scale the observations to match the model's equivalent width:
    out, eout = an.lsq((observation_norm, np.ones(nobs*nk)), flatmodel, w=w_observation)
    sc_observation_norm = observation_norm * out[0] + out[1]
    fitargs = (sc_observation_norm, w_observation, Rmatrix, 0)
    #perfect_fit = optimize.fmin_tnc(dime.entropy_map_norm_sp, flatguess, fprime=dime.getgrad_norm_sp, args=fitargs, maxfun=10000, bounds=[(1e-6, 20)]*nx, ftol=10, disp=5)
    perfect_fit = an.gfit(dime.entropy_map_norm_sp, flatguess, fprime=dime.getgrad_norm_sp, args=fitargs, ftol=ftol, disp=1, maxiter=1e4, bounds=bounds)
    
    perfect_model = dime.normalize_model(np.dot(perfect_fit[0], Rmatrix), nk)
    w_observation /=  w_observation.max() * (sc_observation_norm - perfect_model)[w_observation>0].std()**2

fitargs = (sc_observation_norm, w_observation, Rmatrix, alpha)
bfit = an.gfit(dime.entropy_map_norm_sp, flatguess, fprime=dime.getgrad_norm_sp, args=fitargs, ftol=ftol, disp=1, maxiter=1e4, bounds=bounds)
          #bfit = optimize.fmin_tnc(dime.entropy_map_norm_sp, flatguess, fprime=dime.getgrad_norm_sp, args=fitargs, maxfun=10000, bounds=[(1e-6, 20)]*nx, ftol=1, disp=5)
allfits.append(bfit)
bestparams = bfit[0]
model_observation = dime.normalize_model(np.dot(bestparams, Rmatrix), nk)
metric, chisq, entropy = dime.entropy_map_norm_sp(bestparams, *fitargs, retvals=True)

# reshape into grid
bestparamgrid = np.reshape(bestparams, (-1, nlon))

fits.writeto('map_IGRINS_order4and5_16Nov2021_singleRVshift.fits', bestparamgrid)

#hiphi = np.linspace(0, 2*np.pi, nlat*5)#arange(nlat*5)*2*np.pi/nlat/5 #
#hitheta = np.linspace(0, np.pi, nlon*5)#np.arange(nlon*5)*np.pi/nlon/5 #
#grid_hiphi, grid_hitheta = np.meshgrid(hiphi, hitheta)
#himap = interpolate.griddata(np.concatenate((phi_corner-2*np.pi, phi_corner, phi_corner+2*np.pi)).ravel(), np.concatenate((theta_corner, theta_corner, theta_corner)).ravel(), np.concatenate((np.tile(bestparams, (4,1)).T, np.tile(bestparams, (4,1)).T, np.tile(bestparams, (4,1)).T)).ravel(), (grid_hiphi, grid_hitheta))


fig = plt.figure()
ax = fig.add_subplot(111, projection='mollweide')

lon = np.linspace(-np.pi, np.pi, nlon)
lat = np.linspace(-np.pi/2., np.pi/2., nlat)
Lon,Lat = np.meshgrid(lon,lat)

im = ax.pcolormesh(Lon,Lat,bestparamgrid, cmap=plt.cm.plasma)

plt.show()

fig.savefig('map_WISE1049B_IGRINS_Mollweide_16Nov21.png')


    
sys.exit()
    
