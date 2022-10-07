from operator import mod
import re
import numpy as np
from astropy.table import Table
import modelfitting as mf
from scipy import signal
import pickle
import matplotlib.pyplot as plt

datafile = "../data/fainterspectral-fits_6.pickle"
modelfile = "../data/BT-Settl_lte015-5.0-0.0+0.0_orig.fits"
telluricfile = ""

include_telluric = False
fit_time_averaged = True
NUM_WAVELEN_COEFFS = 4
NUM_CONTINU_COEFFS = 2

class FittingTask():
    def __init__(self) -> None:
        self.data = None
        self.include_telluric = False
        self.fit_time_averaged = True
        self.NUM_WAVELEN_COEFFS = 4
        self.NUM_CONTINU_COEFFS = 2
        self.INIT_PARAMS = [21, 0.3, 9e-5]
        self.lam_obs = None
        self.f_obs = None
        self.error_obs = None
        self.weight_obs = None
        self.lam_template = None
        self.f_template = None
        self.lam_telluric = None
        self.f_telluric = None
        self.npix = None
        self.pix = None
        self.wavelen_coeffs = None
        self.func_wavelen = None
        self.continu_coeffs = None
        self.func_continu = None

    def from_archive(self, datafile, plot=False):
        """Creates a fitting task from archive pickle data"""

        with open(datafile, 'rb') as f:
            self.data = pickle.load(f, encoding='latin1')

        self.lam_obs = self.data['wobs']

        # get time-integrated filtered data
        avg_obs = np.median(self.data['obs0'], axis=0)*14.
        self.f_obs = np.vstack([signal.medfilt(avg_obs[band], 3) for band in range(4)])
        
        # get data for individual timeframes, TODO: try if fits better
        self.f_obs_times = self.data['obs0']

        # get errors of data, TODO: try other error definitions
        self.error_obs = np.median(self.data['obs0'], axis=0)*np.sqrt(14.) 
        self.error_obs /= np.median(self.f_obs, axis=1).reshape(4,1)

        # normalize
        self.f_obs /= np.median(self.f_obs, axis=1).reshape(4,1)

        # obtain weights for chi^2 minimization
        self.weight_obs = 1./self.error_obs**2
        wind = np.isinf(self.weight_obs)
        self.weight_obs[wind] = 0.
        
        self.npix = self.lam_obs.shape(1)
        self.pix = np.arange(self.npix, dtype=float) / self.npix
    
    def plot_obs(self, time_averaged=True):
        if time_averaged:
            fig = plt.figure(figsize=(12,12))
            for band in np.arange(4):
                ax = fig.add_subplot(4,1,band+1)
                ax.plot(self.lam_obs[band, :], self.f_obs[band, :], color='black', 
                    label='observed filtered spectrum (fobs0)')
                ax.plot(self.lam_obs[band, :], self.weight_obs[band, :], color='red', 
                    label='chi^2 weights (wfobs0)')
                if band==0:
                    plt.legend(loc=1)
        else:
            raise NotImplementedError("not yet") # TODO: add individual times plotting
        return

    def load_model(self, modelfile, band=0, plot=False):
        model = Table.read(modelfile, format='fits')
        model['wl'] = model['Wavelength']
        model['flux'] = model['Flux']

        lowerlim = self.lam_obs[band].min() - 0.003
        upperlim = self.lam_obs[band].max() + 0.003
        tind = (model['wl'] > lowerlim) * (model['wl'] < upperlim)
        self.lam_template = model['wl'][tind]
        self.f_template = model['flux'][tind]
        self.f_template /= np.median(self.f_template) # normalize

        if plot:
            fig = plt.figure(figsize=(12,3))
            plt.plot(self.lam_template, self.f_template, label='template spectra')
            plt.legend(loc=1)
        return

    def fit_wavelength_coeffs(self, band=0, plot=False):
        """
        Return coefficients for polynomial of order NPW that convert pixel
        number to wavelength.
        """
        wcoef = np.polyfit(self.pix, self.lam_obs[band], self.NUM_WAVELENTH_COEFFS-1)
        self.func_wavelen = np.poly1d(wcoef)

        if plot:
            fig = plt.figure(figsize=(5,3))
            plt.plot(self.pix[:50], self.lam_obs[band][:50], ".", 
                label="observed wavelength points (first 50)")
            plt.plot(self.pix[:50], self.func_wavelen(self.pix[:50]), "-", 
                label="fitted wavelength solution")
            plt.title("What is wavelength solution")
            plt.xlabel("relative pixel number")
            plt.ylabel("observed wavelengths")
            plt.text(0.003, 2.28755, (
                f"f(x)={wcoef[0]:.2e}x$^3${wcoef[1]:.2e}x$^2$"
                f"+{wcoef[2]:.2e}x+{wcoef[3]:.2e}"))
            plt.legend()
        self.wavelen_coeffs = wcoef
        return

    def fit_continuum_coeffs(self, band=0, plot=False):
        self.continu_coeffs = [-0.1, 1.2/np.median(self.lam_template)]   
        # continuum coefficients, TODO: define these in the beginning somehow

        npix = self.lam_obs.shape[1]
        pix = np.arange(npix, dtype=float)/npix
        
        ind90 = np.sort(self.f_obs[band])[int(0.9*npix)]  # fit top 10% (i.e. unobsorbed)
        ccoef = np.polyfit(pix[self.f_obs[band]>ind90], self.f_obs[band][self.f_obs[band]>ind90], NPC-1)
        self.func_continu = np.poly1d(ccoef)

        if plot:
            fig = plt.figure(figsize=(12,3))
            plt.plot(pix, self.f_obs[band], label="observed spectrum")
            plt.plot(pix[self.f_obs[band]>ind90], self.f_obs[band][self.f_obs[band]>ind90], ".", 
                label="largest 10% data points")
            plt.plot(pix, self.func_continu(pix), label="fitted continuum polynomial")
            plt.text(0.6, 1.3, f"f(x)={ccoef[0]:.4f}x+{ccoef[1]:.4f}", color="tab:green")
            plt.legend()
            plt.show()
        return ccoef

    def chi_square_fit(self, band=0):
        parametric_spectrum_func = mf.modelspec_template
        errfunc_args = (
            parametric_spectrum_func, 
            self.lam_template, self.f_template, 
            self.NUM_WAVELEN_COEFFS, self.NUM_CONTINU_COEFFS, self.npix, 
            self.f_obs[band], self.weight_obs[band]
        )
        self.fit = mf.fmin(
            func=mf.errfunc, 
            x0=np.concatenate((self.INIT_PARAMS, self.wavelen_coeffs, self.continu_coeffs)), 
            args=errfunc_args, 
            full_output=True, 
            disp=True, 
            maxiter=10000, 
            maxfun=10000
        )
        self.bestfit_params = self.fit[0]
        print("Best-fit parameters:", self.bestfit_params)

def get_nobroad_spectrum():
    return


if __name__ == "__main__":
    fitTask = FittingTask.from_archive(datafile="../data/fainterspectral-fits_6.pickle")
    fitTask.plot_obs()
    
    # for band in range(4):
    #     lam_model, f_model = load_model(modelfile, plot=True)

    #     wavelenth_coeffs = fit_wavelength_coeffs(NUM_WAVELENTH_COEFFS)
    #     continuum_coeffs = fit_continuum_coeffs(NUM_CONTINUUM_COEFFS, f_model, f_obs, plot=True)



