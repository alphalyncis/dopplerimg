import numpy as np
from astropy.table import Table
import modelfitting as mf
from scipy import signal
import pickle
import matplotlib.pyplot as plt
from astropy.io import fits
import glob

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
        self.lam_atmo = None
        self.f_atmo = None
        self.npix = None
        self.pix = None
        self.wavelen_coeffs = None
        self.func_wavelen = None
        self.continu_coeffs = None
        self.func_continu = None
        self.comment = None

    @classmethod
    def from_archive(cls, datafile):
        """Creates a fitting task from archive pickle data"""

        self = cls()
        self.comment = f"Fitting data set {datafile}."

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
        
        self.npix = self.lam_obs.shape[1]
        self.pix = np.arange(self.npix, dtype=float) / self.npix

        return self


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

    def load_model(self, band, modelfile, plot=False):
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

    def load_telluric(self,  band, telluricfile, plot=False):
        """Read the telluric spectrum and cutout a section covering the required band."""

        atm0 = fits.getdata(telluricfile)

        lowerlim = self.lam_obs[band].min() - 0.003
        upperlim = self.lam_obs[band].max() + 0.003

        aind = (atm0[:, 0]>lowerlim) * (atm0[:, 0] < upperlim)
        self.lam_atmo = atm0[aind, 0]
        self.f_atmo = atm0[aind, 1]

    def fit_wavelength_coeffs(self, band, plot=False):
        """
        Return coefficients for polynomial of order NUM_WAVELEN_COEFFS that convert 
        pixel number to wavelength.
        """
        wcoef = np.polyfit(self.pix, self.lam_obs[band], self.NUM_WAVELEN_COEFFS-1)
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
        return wcoef

    def fit_continuum_coeffs(self, band, plot=False):
        self.continu_coeffs = [-0.1, 1.2/np.median(self.lam_template)]   
        # continuum coefficients, TODO: define these in the beginning somehow

        npix = self.lam_obs.shape[1]
        pix = np.arange(npix, dtype=float)/npix
        
        ind90 = np.sort(self.f_obs[band])[int(0.9*npix)]  # fit top 10% (i.e. unobsorbed)
        ccoef = np.polyfit(pix[self.f_obs[band]>ind90], self.f_obs[band][self.f_obs[band]>ind90], 
            self.NUM_CONTINU_COEFFS-1)
        self.func_continu = np.poly1d(ccoef)
        self.continu_coeffs = ccoef

        if plot:
            fig = plt.figure(figsize=(12,3))
            plt.plot(pix, self.f_obs[band], label="observed spectrum")
            plt.plot(pix[self.f_obs[band]>ind90], self.f_obs[band][self.f_obs[band]>ind90], ".", 
                label="largest 10% data points")
            plt.plot(pix, self.func_continu(pix), label="fitted continuum polynomial")
            plt.text(0.6, 1.3, f"f(x)={ccoef[0]:.4f}x+{ccoef[1]:.4f}", color="tab:green")
            plt.legend()
            plt.show()

    def chi_square_fit(self, band, include_telluric=False, telluricfile=None):
        """
        Perform the actual chi_sq fitting for one piece of spectrum (specified by band). 
        If include telluric, a telluric file name must be passed.
        """

        if include_telluric:
            try:
                self.load_telluric(telluricfile=telluricfile, band=band)
            except:
                raise FileNotFoundError("Please pass a valid file path for the telluric file.")
            parametric_spec_func = mf.modelspec_tel_template
            errfunc_args = (
                parametric_spec_func, 
                self.lam_template, self.f_template,
                self.lam_atmo, self.f_atmo,
                self.NUM_WAVELEN_COEFFS, self.NUM_CONTINU_COEFFS, self.npix, 
                self.f_obs[band], self.weight_obs[band]
            )
        else:
            parametric_spectrum_func = mf.modelspec_template
            errfunc_args = (
                parametric_spectrum_func, 
                self.lam_template, self.f_template, 
                self.NUM_WAVELEN_COEFFS, self.NUM_CONTINU_COEFFS, self.npix, 
                self.f_obs[band], self.weight_obs[band]
            )
        result = mf.fmin(
            func=mf.errfunc, 
            x0=np.concatenate((self.INIT_PARAMS, self.wavelen_coeffs, self.continu_coeffs)), 
            args=errfunc_args, 
            full_output=True, 
            disp=True, 
            maxiter=10000, 
            maxfun=10000
        )
        return result

    def get_nobroad_spectrum(self):
        return

    def fit_one_band(self, band, modelfile, include_telluric=False, telluricfile=None):
        self.load_model(band=band, modelfile=modelfile)
        self.fit_wavelength_coeffs(band=band)
        self.fit_continuum_coeffs(band=band)
        fit = self.chi_square_fit(band=band, include_telluric=include_telluric, telluricfile=telluricfile)
        return fit
    
    def fit_four_bands(self, modelfile, include_telluric=False, telluricfile=None, plot=True):
        res = {"spec": [], "lam": [], "fit":[]}
        for j in range(4):
            self.load_model(band=j, modelfile=modelfile)
            self.fit_wavelength_coeffs(band=j)
            self.fit_continuum_coeffs(band=j)
            fit = self.chi_square_fit(band=j, include_telluric=include_telluric, telluricfile=telluricfile)
            if include_telluric:
                spec, lam = mf.modelspec_tel_template(
                    fit[0], 
                    self.lam_template, self.f_template, 
                    self.lam_atmo, self.f_atmo,
                    self.NUM_WAVELEN_COEFFS, self.NUM_CONTINU_COEFFS, 
                    self.npix, retlam=True
                )
            else:
                spec, lam = mf.modelspec_template(
                    fit[0], 
                    self.lam_template, self.f_template, 
                    self.NUM_WAVELEN_COEFFS, self.NUM_CONTINU_COEFFS, 
                    self.npix, retlam=True
                )
            res['spec'].append(spec)
            res['lam'].append(lam)
            res['fit'].append(fit)
            # res.append(spec, lam, fit)
        if plot:
            fig = plt.figure(figsize=(15,11))         
            for j in np.arange(4):
                ax = fig.add_subplot(4,1,j+1)
                ax.plot(self.lam_obs[j, :], self.f_obs[j, :], color='black', label="observation")
                if include_telluric:
                    ax.plot(self.lam_obs[j, :], res['spec'][j], color='red', label="best-fit (with telluric")
                else:
                    ax.plot(self.lam_obs[j, :], res['spec'][j], color='red', label="best-fit")
                # ax.plot(self.lam_obs[j, :], chipmodnobroad[j, :], color='red', alpha=0.1, label="best-fit but no broaden")
                ax.plot(self.lam_obs[j, :], self.data['chipmods'][0, j, :], color='green', 
                    alpha=0.4, linewidth=1.5, label="original best fit")
            plt.legend()
        return res

    def fit_many_models(self, models_dir):

        return


if __name__ == "__main__":
    datafile = "../data/fainterspectral-fits_6.pickle"
    modelfile = "../data/BT-Settl_lte015-5.0-0.0+0.0_orig.fits"
    models_dir = "../data/BT-SettlModels"
    telluricfile = "../data/transdata_0,5-14_mic_hires.fits"
    fitTask = FittingTask.from_archive(datafile)
    # fitTask.plot_obs()
    res1 = fitTask.fit_one_band(band=0, modelfile=modelfile, include_telluric=True, telluricfile=telluricfile)
    res4 = fitTask.fit_four_bands(modelfile, include_telluric=True, telluricfile=telluricfile)


    # for band in range(4):
    #     lam_model, f_model = load_model(modelfile, plot=True)

    #     wavelenth_coeffs = fit_wavelength_coeffs(NUM_WAVELENTH_COEFFS)
    #     continuum_coeffs = fit_continuum_coeffs(NUM_CONTINUUM_COEFFS, f_model, f_obs, plot=True)



