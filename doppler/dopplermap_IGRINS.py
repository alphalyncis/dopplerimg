import dia3 as dia
import dao3 as dao
import nsdata3 as ns
import analysis3 as an
import ELL_map_class as ELL_map
import dime3 as dime # Doppler Imaging & Maximum Entropy 

# import packages
import numpy as np
from scipy import signal
import os
from scipy import interpolate
from astropy.io import fits
import matplotlib.pyplot as plt
import pickle
import glob
import os
homedir = os.path.expanduser('~')

nlat, nlon = 20, 40
datdir = f'{homedir}/workspace/dopplerimg/code/doppler/'
suffix = 'order1-5'
nobs = 14
nchips = 5
firstchip = 1

with open(f'{homedir}/uoedrive/result/CIFIST/IGRINS_W1049B_K_lte015.0-5.0.pickle', 'rb') as f:
    ret = pickle.load(f, encoding="latin1")
obs1 = ret['fobs0']
chiplams = ret['chiplams']
chiplams *= 10000 #convert microns to angstroms
chipmodnobroad = ret['chipmodnobroad']
chipcors = ret['chipcors']
chipmods = ret['chipmods']

# W1049B rotation parameters
#set period and number of observations
period = 4.87 # period in hours
per = period*3600 # period in seconds
vsini = 29e3 # m/s
inc = 0.3491 # rad

# Specify user-defined options:
LLD = 1.0
alpha = 4500
nk = 103 # parameter for LSD
ftol = 0.01 # tolerance for convergence of maximum-entropy
nstep = 1500

##################################################
# Begin!
##################################################

#Create the MJD dates
mjd_0 = np.linspace(0,period,nobs)*(1/24) 
# observation is 5 hours long, 56 intervals, in julian dates, 
# assume intervals are equally spaced across the observation
phase = (mjd_0)*86400 * 2*np.pi/ per # 0 ~ 2*pi
print(mjd_0)
print('phase: ', phase)

##################################
#LEAST SQUARES DECONVOLUTION to build the composite line profile
#####################################    

# Set up for Least Squares Deconvolution
modelfn = f'{datdir}/lte015-5.0-0.0a+0.0.BT-Settl.spec.7.fits'
pmod = os.path.split(modelfn)[1].replace('BT-Settl.spec.7.fits', '_dao').replace('.',',')
f_linelist   = f'{datdir}{pmod}_edited.clineslsd'
# Load the LSD files:
(lineloc, lineew, linespec) = dao.getlines(f_linelist)

f_pspec      = f'{datdir}{pmod}.fits'
f_pspec_cont = f'{datdir}{pmod}C.fits'
pspec_cont = fits.getdata(f_pspec_cont) 
hdr_pspec_cont = fits.getheader(f_pspec_cont)
wspec = hdr_pspec_cont['crval1'] + np.arange(pspec_cont.size)*hdr_pspec_cont['cdelt1']
spline = interpolate.UnivariateSpline(wspec, pspec_cont, s=0.0, k=1.0) #set up interpolation over the continuum measurement

# Compute LSD:
kerns = np.zeros((nobs, nchips, nk), dtype=float)
modkerns = np.zeros((nobs, nchips, nk), dtype=float)
### In this scheme, run LSD separately for each frame's wavelength solution:
for i, jj in enumerate(range(firstchip, firstchip+nchips)): # EB: chip 4,5
    for kk in range(nobs): # EB: frame 0-13
        # EB: create delta-function line spectrum from wavelength grid, list of line locations, list of equivalent widths
        deltaspec = ns.linespec(lineloc*(1.+8e-5), lineew, chiplams[:,jj].mean(0), verbose=False, cont=spline(chiplams[:,jj].mean(0)))
        # lines shifted to rv by 1+9e-5

        # EB: DSA=Difference Spectral Analysis. Match given spectra to reference spectra
        m,kerns[kk,i],b,c = dia.dsa(deltaspec, obs1[kk,jj]/chipcors[kk,jj], nk)
        m,modkerns[kk,i],b,c = dia.dsa(deltaspec, chipmodnobroad[kk,jj]/chipcors[kk,jj], nk)    
		#EB: dia.dsa returns: Response matrix, kernel used in convolution, background offset, chisquared of fit
		#EB: only the kernels are used from here on   

# Compute LSD velocity grid:
dbeta = np.diff(chiplams).mean()/chiplams.mean() # EB: np.diff returns array containing the difference between elements of the array provided. so dbeta will be mean difference/mean value.
dx = -dbeta * np.arange(np.floor(-nk/2.+.5), np.floor(nk/2.+.5)) # EB: edit np.arange
dv = an.c*dx / 1e3 # km/s # EB: an.c is the speed of light

intrinsic_profile = modkerns[:,:].mean(0).mean(0)
intrinsic_profile -= intrinsic_profile[0] + (intrinsic_profile[-1] - intrinsic_profile[0])/(nk - 1.) * np.arange(nk)
systematic_rv_offset = (intrinsic_profile==intrinsic_profile.max()).nonzero()[0] - (dv==0).nonzero()[0]
intrinsic_profile = np.interp(np.arange(nk), np.arange(nk) - systematic_rv_offset, intrinsic_profile)
kerns = np.array([[np.interp(np.arange(nk), np.arange(nk) - systematic_rv_offset, kerns[jj, kk]) for kk in range(nchips)] for jj in range(nobs)])

###################################################################################
# Prepare the Doppler Imaging "response matrix"
######################################################################################

modIP = 1. - np.concatenate((np.zeros(300), intrinsic_profile, np.zeros(300)))
modDV = - np.arange(np.floor(-modIP.size/2.+.5), np.floor(modIP.size/2.+.5)) * dbeta * an.c / 1e3
flineSpline2 = interpolate.UnivariateSpline(modDV[::-1], modIP[::-1], k=1., s=0.)

# XC: concatenate the kerns of 14 obs and average over the 2 orders
observation = 1. - kerns.mean(axis=1).ravel()
observation_norm = observation.copy()
for ii in range(nobs): # EB: frame 0-13
    i0, i1 = ii*nk, (ii+1)*nk # EB: nk is the number of pixels in LSD computation on this run through
    inds = np.concatenate((np.arange(i0, i0+7), np.arange(i1-7, i1)))
    continuumfit = np.polyfit(inds, observation_norm[inds], 1)
    observation_norm[i0:i1] /= np.polyval(continuumfit, np.arange(i0, i1))
plt.figure(figsize=(8,1.5))
plt.plot(observation)
plt.plot(observation_norm)

###################################################################################
# Start the Imaging Analysis
###################################################################################

nx = nlat*nlon
dime.setup(observation_norm.size, nk)
flatguess = 100*np.ones(nx)
bounds = [(1e-6, 300)]*nx
allfits = []

# initialize Doppler map object
mmap = ELL_map.map(nlat=nlat, nlon=nlon, i=inc) #ELL_map.map returns a class object 
ncell = mmap.ncell
latlon_corners = np.array([c.corners_latlon for c in mmap.cells]) #corner coordinates
phi_corner = latlon_corners[:,1,:]
theta_corner = latlon_corners[:,0,:]



Rmatrix = np.zeros((ncell, nobs*dv.size), dtype=np.float32) 

for kk, rot in enumerate(phase):
	speccube = np.zeros((ncell, dv.size), dtype=np.float32) 
	this_map = ELL_map.map(nlat=nlat, nlon=nlon, i=inc, deltaphi=-rot) #ELL_map.map returns a class object 
	
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
    perfect_fit = an.gfit(dime.entropy_map_norm_sp, flatguess, fprime=dime.getgrad_norm_sp, args=fitargs, ftol=ftol, disp=1, maxiter=1e4, bounds=bounds)
    
    perfect_model = dime.normalize_model(np.dot(perfect_fit[0], Rmatrix), nk)
    w_observation /=  w_observation.max() * (sc_observation_norm - perfect_model)[w_observation>0].std()**2


fitargs = (sc_observation_norm, w_observation, Rmatrix, alpha)
bfit = an.gfit(dime.entropy_map_norm_sp, flatguess, fprime=dime.getgrad_norm_sp, args=fitargs, ftol=ftol, disp=1, maxiter=1e4, bounds=bounds)
allfits.append(bfit)
bestparams = bfit[0]
model_observation = dime.normalize_model(np.dot(bestparams, Rmatrix), nk)
metric, chisq, entropy = dime.entropy_map_norm_sp(bestparams, *fitargs, retvals=True)

# reshape into grid
bestparamgrid = np.reshape(bestparams, (-1, nlon))

fits.writeto(f'map_IGRINS_{suffix}.fits', bestparamgrid, overwrite=True)

fig = plt.figure()
ax = fig.add_subplot(111, projection='mollweide')
lon = np.linspace(-np.pi, np.pi, nlon)
lat = np.linspace(-np.pi/2., np.pi/2., nlat)
Lon,Lat = np.meshgrid(lon,lat)
im = ax.pcolormesh(Lon,Lat,bestparamgrid, cmap=plt.cm.plasma)
fig.colorbar(im)
plt.show()
fig.savefig(f'map_IGRINS_{suffix}.png')