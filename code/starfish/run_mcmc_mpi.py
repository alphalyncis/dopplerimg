import sys
import time
import emcee
import numpy as np
from astropy.io import fits
from scipy import signal
import glob
import os
homedir = os.path.expanduser('~')
from schwimmbad import MPIPool
from Starfish.models import SpectrumModel
from Starfish.spectrum import Spectrum
import scipy.stats as st

# get data
filelist = sorted(glob.glob(f'{homedir}/uoedrive/data/IGRINS/SDCK*_1f.spec.fits'))

fluxes = []
wls = []

for filename in filelist:
    wlname = filename.split('_1f')[0]+'.wave.fits'

    flux = fits.getdata(filename)
    wl = fits.getdata(wlname)
    hdr = fits.getheader(filename)
    # trim first and last 100 columns
    flux = flux[:, 100:1948]
    wl = wl[:, 100:1948]
    
    fluxes.append(flux)
    wls.append(wl)
dims = np.array(fluxes).shape
fluxes = np.array(fluxes)
wls = np.array(wls)

obs0 = np.median(fluxes[:, 1:, :], axis=0)*dims[0]  # make mean spectrum
eobs0 = np.median(fluxes[:, 1:, :], axis=0)*np.sqrt(dims[0])  # make mean noise spectrum (assuming just photon noise)
fobs0 = np.vstack([signal.medfilt(obs0[jj], 3) for jj in range(dims[1]-1)])  # smooth spectrum
# fix wavelength array to have same shape
wls = wls[:, 1:24, :]

obs=22
order=5
wl_order = wls[obs, order, :] *10000
flux_order = fobs0[order]
data = Spectrum(wl_order, flux_order, sigmas=None, masks=None, name="Spectrum")


priors = {
    "T": st.uniform(500,2500),
    "vsini": st.uniform(0, 500),
    "vz": st.uniform(10, 100),
    #"cheb:1": st.norm(1, 0.1)
}

model = SpectrumModel(
    emulator="BTSettl_K_emu.hdf5",
    data=data,
    grid_params=[1500, 5.0],
    vsini=25
)
model.load("example_MAP_1.toml")

nwalkers = 100
ndim = len(model.labels)
nsteps = 100

# Initialize gaussian ball for starting point of walkers
scales = {"T": 1, "vsini": 1, "vz": 1}
ball = np.random.randn(nwalkers, ndim)
for i, key in enumerate(model.labels):
    ball[:, i] *= scales[key]
    ball[:, i] += model[key]

def log_prob(P, priors):
    model.set_param_vector(P)
    return model.log_likelihood(priors)
    
def log_prob_example(theta):
    t = time.time() + np.random.uniform(0.005, 0.008)
    while True:
        if time.time() >= t:
            break
    return -0.5*np.sum(theta**2)

with MPIPool() as pool:
    if not pool.is_master():
        pool.wait()
        sys.exit(0)

    backend = emcee.backends.HDFBackend("example_chain.hdf5")
    backend.reset(nwalkers, ndim)

    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, args=(priors,), backend=backend, pool=pool)
        
    start = time.time()
    sampler.run_mcmc(backend.get_last_sample(), nsteps)
    end = time.time()
    print(end - start)