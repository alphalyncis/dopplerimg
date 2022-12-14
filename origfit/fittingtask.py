import numpy as np
from astropy.table import Table
import modelfitting as mf
from scipy import signal
import pickle
import matplotlib.pyplot as plt
from astropy.io import fits
import glob
import os

class FittingTask():
    def __init__(self) -> None:
        self.data = None
        self.include_telluric = False
        self.fit_time_averaged = True
        self.NUM_WAVELEN_COEFFS = 4
        self.NUM_CONTINU_COEFFS = 2
        self.INIT_PARAMS = [21, 0.3, 9e-5]
        self.INIT_PARAMS_TEL = [21, 0.3, 9e-5, 42, 1.3]
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
        self.f_obs = np.vstack([signal.medfilt(avg_obs[jj], 3) for jj in range(4)])
        
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
    
    @classmethod
    def from_IGRINS(cls, datadir, band="K"):
        """Creates a fitting task from IGRINS data"""
        self = cls()
        self.comment = f"Fitting data set {datafile}."

        #  Open data
        if band == "K":
            filelist = sorted(glob.glob(f'{datadir}/SDCK*_1f.spec.fits'))
        elif band == "H":
            filelist = sorted(glob.glob(f'{datadir}/SDCH*_1f.spec.fits'))
        
        fluxes = []
        wls = []
        for filename in filelist:
            wlname = filename.split('_1f')[0]+'.wave.fits'

            flux = fits.getdata(filename)
            wl = fits.getdata(wlname)

            # trim first and last 100 columns
            flux = flux[:, 100:1948]
            wl = wl[:, 100:1948]

            fluxes.append(flux)
            wls.append(wl)
        
        # Collapse data over the whole observation into one mean, smoothed spectrum
        dims = np.array(fluxes).shape
        fluxes = np.array(fluxes)
        wls = np.array(wls)
        
        # remove first order with all NaNs when making various arrays
        avg_obs = np.median(fluxes[:, 1:, :], axis=0)*dims[0]  # make mean spectrum
        self.f_obs = np.vstack([signal.medfilt(avg_obs[jj], 3) for jj in range(dims[1]-1)])  # smooth spectrum
        
        # get errors of data
        self.error_obs = np.median(fluxes[:, 1:, :], axis=0)*np.sqrt(dims[0])  # make mean noise spectrum (assuming just photon noise)
        self.error_obs /= np.nanmedian(self.f_obs, 1).reshape(dims[1]-1,1)
        
        # normalize
        self.f_obs /= np.nanmedian(self.f_obs, 1).reshape(dims[1]-1,1)
        self.weight_obs = 1./self.error_obs**2   # make noise spectrum
        #ind = np.where(np.isinf(wfobs0)) # set infinite values in noise spectrum to NaNs
        #wfobs0[ind] = np.nan
        wind = np.isinf(self.weight_obs)  # remove points with infinite values
        self.weight_obs[wind] = 0.

        # fix wavelength array to have same shape
        wls = wls[:, 1:24, :]

        #TODO: to be done



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
        if modelfile.split(".")[-1] == "fits":
            model = fits.getdata(modelfile)
            model = np.array(model.tolist())

        else:
            raise NotImplementedError("can only read fits file.")
            # TODO: read ascii format files

        if model[:,0].min() > 1000: # if unit in angstrom
            model[:,0] /= 10000

        lowerlim = self.lam_obs[band].min() - 0.003
        upperlim = self.lam_obs[band].max() + 0.003
        tind = (model[:,0]>lowerlim) * (model[:,0] < upperlim)
        self.lam_template = model[tind, 0]
        self.f_template = model[tind, 1]
        self.f_template /= np.median(self.f_template) # normalize
        
        # to be compatible with rotational profile convolution kernel
        if len(self.lam_template) < 400:
           new_lam = np.linspace(self.lam_template[0], self.lam_template[-1], 400)
           self.f_template = np.interp(new_lam, self.lam_template, self.f_template)
           self.lam_template = new_lam

        if plot:
            fig = plt.figure(figsize=(12,3))
            plt.plot(self.lam_template, self.f_template, label='template spectra')
            plt.legend(loc=1)
        return self.lam_template, self.f_template

    def load_telluric(self,  band, telluricfile, plot=False):
        """
        Read the (already normalized) telluric spectrum and cutout a segment 
        covering the required band.
        """

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
            guess = np.concatenate((self.INIT_PARAMS_TEL, self.wavelen_coeffs, self.continu_coeffs))
        else:
            parametric_spectrum_func = mf.modelspec_template
            errfunc_args = (
                parametric_spectrum_func, 
                self.lam_template, self.f_template, 
                self.NUM_WAVELEN_COEFFS, self.NUM_CONTINU_COEFFS, self.npix, 
                self.f_obs[band], self.weight_obs[band]
            )
            guess = np.concatenate((self.INIT_PARAMS, self.wavelen_coeffs, self.continu_coeffs))
        result = mf.fmin(
            func=mf.errfunc, 
            x0=guess,
            args=errfunc_args, 
            full_output=True, 
            disp=True, 
            maxiter=10000, 
            maxfun=10000
        )
        return result

    def get_nobroad_spectrum(self):
        return

    def fit_one_band(self, band, modelfile, include_telluric=False, telluricfile=None, plot=True):
        self.load_model(band=band, modelfile=modelfile)
        self.fit_wavelength_coeffs(band=band)
        self.fit_continuum_coeffs(band=band)
        fit = self.chi_square_fit(band=band, include_telluric=include_telluric, telluricfile=telluricfile)
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
        if plot:
            plt.figure(figsize=(15,4))
            plt.plot(self.lam_obs[band, :], self.f_obs[band, :], color='black', label="observed spectrum")
            plt.plot(self.lam_obs[band, :], spec, color='red', label="best-fit spectrum")
            plt.legend()
        return {"spec": [spec], "lam": [lam], "fit": [fit]}
    
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
                    ax.plot(self.lam_obs[j, :], res['spec'][j], color='red', label="best-fit (with telluric)")
                else:
                    ax.plot(self.lam_obs[j, :], res['spec'][j], color='red', label="best-fit")
                # ax.plot(self.lam_obs[j, :], chipmodnobroad[j, :], color='red', alpha=0.1, label="best-fit but no broaden")
                ax.plot(self.lam_obs[j, :], self.data['chipmods'][0, j, :], color='green', 
                    alpha=0.4, linewidth=1.5, label="original best fit")
            plt.legend()
        return res

    def fit_many_models(self, models_dir, include_telluric=False, telluricfile=None):
        result = {}
        for modelfile in glob.glob(f"{models_dir}/*[0-9]*.*"):
            modelname = modelfile.split("/")[-1]
            print(f"fitting model {modelname}...")

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

            result[modelname] = res
        return result


def plot_one_result(ft, result):
    fig = plt.figure(figsize=(15,11))         
    for j in np.arange(4):
        ax = fig.add_subplot(4,1,j+1)
        ax.plot(ft.lam_obs[j, :], ft.f_obs[j, :], color='black', label="observation")
        ax.plot(ft.lam_obs[j, :], result['spec'][j], color='red', label="best-fit (with telluric)")
        ax.plot(ft.lam_obs[j, :], ft.data['chipmods'][0, j, :], color='green', 
            alpha=0.4, linewidth=1.5, label="original best fit from pickle")
    plt.legend()


def print_all_results(result, avg_bands=True):
    for key, item in result.items():
        if avg_bands:
            print(f"model:{key} fmin:{np.array([item['fit'][j][1] for j in range(4)]).mean()}")
        else:
            print(f"model:{key} fmin:{[item['fit'][j][1] for j in range(4)]}")



if __name__ == "__main__":
    homedir = os.path.expanduser('~')
    data_dir = f"{homedir}/uoedrive/data/"
    #data_dir = "/~/workspace/dopplerimg/data/"
    datafile = "CRIRES/fainterspectral-fits_6.pickle"
    modelfile_ian_orig = "BT-Settl_lte015-5.0-0.0+0.0_orig.fits"
    #models_dir = "BT-SettlModels/015-5.0s/bestfits"
    models_dir = "CallieModels"
    telluricfile = "telluric/transdata_0,5-14_mic_hires.fits"
    ft = FittingTask.from_archive(data_dir+datafile)
    # fitTask.plot_obs()
    # res1 = ft.fit_one_band(band=0, modelfile=modelfile, include_telluric=True, telluricfile=telluricfile)
    res4 = ft.fit_four_bands(
        modelfile=data_dir+"CallieModels/t1500g1000nc_m0.0_co1.0.fits", 
        include_telluric=True, 
        telluricfile=data_dir+telluricfile)
    
    # resall = ft.fit_many_models(data_dir+models_dir, include_telluric=True, telluricfile=data_dir+telluricfile)

    # for band in range(4):
    #     lam_model, f_model = load_model(modelfile, plot=True)

    #     wavelenth_coeffs = fit_wavelength_coeffs(NUM_WAVELENTH_COEFFS)
    #     continuum_coeffs = fit_continuum_coeffs(NUM_CONTINUUM_COEFFS, f_model, f_obs, plot=True)



