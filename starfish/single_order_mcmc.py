import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from scipy import signal
import glob
from astropy.table import Table
import modelfitting as mf
import os
import sys
homedir = os.path.expanduser('~')
from Starfish.spectrum import Spectrum
from Starfish.models import SpectrumModel
import emcee
import scipy.stats as st
import arviz as az
import corner

def run(target, band, order, suffix=None, test=False):
    starfishdir = f"{homedir}/dopplerimg/starfish"
    resultdir = f"{starfishdir}/run_IGRINS_{target}_{band}_order{order}"
    if suffix is not None:
        resultdir += suffix
    if not os.path.exists(resultdir):
        os.mkdir(resultdir)
    cwd = os.getcwd()
    print("cwd is", cwd)
    chain_backend_file = f"chain_IGRINS_{target}_{band}_order{order}.hdf5"

    ##############################################################################
    ### Open data
    ##############################################################################
    print(f"Data path is {homedir}/uoedrive/data/IGRINS_{target}/SDC{band}*, order{order}")
    filelist = sorted(glob.glob(f'{homedir}/uoedrive/data/IGRINS_{target}/SDC{band}*_1f.spec.fits'))
    if len(filelist) == 0:
        raise FileNotFoundError("Please check your data path.")
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
    dims = np.array(fluxes).shape
    fluxes = np.array(fluxes)
    wls = np.array(wls)

    obs0 = np.median(fluxes[:, 1:, :], axis=0)*dims[0]  # make mean spectrum
    fobs0 = np.vstack([signal.medfilt(obs0[jj], 5) for jj in range(dims[1]-1)])  # smooth
    # fix wavelength array to have same shape
    wls = wls[:, 1:24, :]
    wl_order = wls[22, order, :]   # take the wl points of obs=22
    flux_order = fobs0[order]

    if test:
        # use Callie's models as synthetic spectra
        modelpath = f"{homedir}/uoedrive/data/CallieModels/t1500g1000f8_m0.0_co1.0.fits"
        model = fits.getdata(modelpath)
        print(f"Testing synthetic spectrum {modelpath}")
        lolim = 1.9
        hilim = 2.5
        tind = (model['wl']>lolim) * (model['wl'] < hilim)
        lam_template = model['wl'][tind]
        template = model['flux'][tind]
        template /= np.median(template)

        npix = 1848
        NPW = 4
        pix = np.arange(npix, dtype=float)/npix
        wcoef = np.polyfit(pix, wls[0, order, :], NPW-1)
        ccoef = [-0.1, 1.2/np.median(template)]
        NPC = len(ccoef)
        ind90 = np.sort(fobs0[order])[int(0.9*npix)]  
        ccoef = np.polyfit(pix[fobs0[order]>ind90], fobs0[order][fobs0[order]>ind90], NPC-1)
        params = [30, 0.9, 20] + list(wcoef) + list(ccoef)
        mymod = mf.modelspec_template(params, lam_template, template, NPW, NPC, npix, retlam=False)
        flux_order = mymod
        # flux_order = np.interp(wl_order, lam_template, template)
    ##############################################################################
    ### Create Starfish Model
    ##############################################################################

    if len(flux_order) == 0:
        raise ValueError("some problem with spectrum data")
        
    data = Spectrum(wl_order*10000, flux_order, sigmas=None, masks=None, name="Spectrum")

    # Create model spectrum
    model = SpectrumModel(
        emulator=f"{starfishdir}/BTSettl_K_emu.hdf5",
        data=data,
        grid_params=[1500, 5.0],
        vsini=25,
        vz=30,
        cheb=[0, 0],
        global_cov=dict(log_amp=0, log_ls=1),
    )
    model.save(f"{resultdir}/model_init.toml")

    # specify priors
    priors = {
        "T": st.uniform(1200,750),
        "logg": st.uniform(4.5,1.0),
        "vsini": st.uniform(0, 500),
        "vz": st.uniform(10, 100),
        "cheb:1": st.uniform(-3, 6),
        "cheb:2": st.uniform(-3, 6),
        "global_cov:log_amp": st.norm(0, 5),
        "global_cov:log_ls": st.uniform(0, 10),
    }
    string = """priors = {
        "T": st.uniform(1200,750),
        "logg": st.uniform(2.5,3.0),
        "vsini": st.uniform(0, 500),
        "vz": st.uniform(10, 100),
        "cheb:1": st.uniform(-3, 6),
        "cheb:2": st.uniform(-3, 6),
        "global_cov:log_amp": st.norm(0, 5),
        "global_cov:log_ls": st.uniform(0, 10),
    }"""
    with open(f"{resultdir}/priors.txt", "w") as f:
        f.write(string)
    print(string)
    print("log_likelihood:", model.log_likelihood(priors))

    # must freeze logg for now otherwise get
    # "ValueError: Querying emulator outside of original parameter range."
    #model.freeze("logg")
    #print("fittable params:", model.labels)  # These are the fittable parameters


    # Numerical optimization before mcmc
    print("Running numerical optimization...")
    model.train(priors, options={"maxiter": 1000})

    model.freeze("global_cov")  # freeze this for speed now
    print("MCMC will fit params:", model.labels)

    model.save(f"{resultdir}/model_optimized.toml")
    model.plot()
    plt.savefig(f"{resultdir}/plot_model.png")

    ##############################################################################
    ### Set up mcmc
    ##############################################################################

    # Set our walkers and dimensionality
    nwalkers = 20
    ndim = len(model.labels)

    # Initialize gaussian ball for starting point of walkers
    ball = np.random.randn(nwalkers, ndim)

    for i, key in enumerate(model.labels):
        ball[:, i] *= model[key]/20 # scale the sigma by 1/10 of param value
        ball[:, i] += model[key]

    # our objective to maximize
    def log_prob(P, priors):
        model.set_param_vector(P)
        return model.log_likelihood(priors)

    # Set up our backend and sampler
    backend = emcee.backends.HDFBackend(f"{resultdir}/{chain_backend_file}")
    backend.reset(nwalkers, ndim)
    sampler = emcee.EnsembleSampler(
        nwalkers, ndim, log_prob, args=(priors,), backend=backend
    )

    ##############################################################################
    ### Set mcmc run
    ##############################################################################

    # Start sampler with a max burn-in of 500 samples
    max_n = 1000

    # We'll track how the average autocorrelation time estimate changes
    index = 0
    autocorr = np.empty(max_n)
    old_tau = np.inf  # will be useful to testing convergence

    # Now we'll sample for up to max_n steps
    print("Running mcmc burn-in...")
    for sample in sampler.sample(ball, iterations=max_n, progress=True, store=True):
        # Only check convergence every 10 steps
        if sampler.iteration % 10:
            continue

        # Compute the autocorrelation time so far
        # Using tol=0 means that we'll always get an estimate even
        # if it isn't trustworthy
        tau = sampler.get_autocorr_time(tol=0)
        autocorr[index] = np.mean(tau)
        index += 1
        # skip math if it's just going to yell at us
        if np.isnan(tau).any() or (tau == 0).any():
            continue
        # Check convergence
        converged = np.all(tau * 10 < sampler.iteration)
        converged &= np.all(np.abs(old_tau - tau) / tau < 0.02)
        if converged:
            print(f"Converged at sample {sampler.iteration}")
            break
        old_tau = tau

    # After our model has converged, let's take a few extra samples
    print("Running mcmc samples...")
    sampler.run_mcmc(backend.get_last_sample(), 500, progress=True, store=True)


    ##############################################################################
    ### MCMC Chain analysis
    ##############################################################################

    reader = emcee.backends.HDFBackend(f"{resultdir}/{chain_backend_file}")
    full_data = az.from_emcee(reader, var_names=model.labels)

    # Plot traces
    az.plot_trace(full_data, figsize=(10,9))
    plt.xlabel("step number")
    plt.tight_layout()
    plt.savefig(f"{resultdir}/plot_traces.png")

    # Discard burn-in samples
    try:
        tau = reader.get_autocorr_time(tol=0)
        print("tau:", tau)
        burnin = int(tau.max())
        thin = int(0.3 * np.min(tau))
    except:
        print("At least 1 param did not converge")
        burnin = 100
        thin = 30
    burn_samples = reader.get_chain(discard=burnin, thin=thin)

    dd = dict(zip(model.labels, burn_samples.T))
    burn_data = az.from_dict(dd)

    # Plot traces after burn-in
    az.plot_trace(burn_data, figsize=(10,9));
    plt.xlabel("step number")
    plt.tight_layout()
    plt.savefig(f"{resultdir}/plot_traces_burned.png")

    # Save summary table
    statistics = az.summary(burn_data)
    statistics.to_csv(f"{resultdir}/statistics.csv")

    # Save posterior plot
    az.plot_posterior(burn_data, ["vsini", "vz", "T", "logg", "cheb:1", "cheb:2"])
    plt.savefig(f"{resultdir}/plot_posterior.png")

    # Save corner plot
    sigmas = ((1 - np.exp(-0.5)), (1 - np.exp(-2)))
    corner.corner(
        burn_samples.reshape((-1, 6)),
        labels=model.labels,
        quantiles=(0.05, 0.16, 0.84, 0.95),
        levels=sigmas,
        show_titles=True
    )
    plt.savefig(f"{resultdir}/plot_corner.png")

    # Save best-fit model
    bestfit_params = dict(az.summary(burn_data)["mean"])
    model.set_param_dict(bestfit_params)

    model.plot()
    plt.savefig(f"{resultdir}/plot_bestfit_model.png")
    model.save(f"{resultdir}/model_bestfit.toml")

if __name__ == "__main__":
    band = "K"
    target = "W1049B"
    try:
        order = int(sys.argv[1]) 
    except:
        raise NameError("Must provide an order number to fit (start from 0)")
    try:
        testflag = sys.argv[2]
        if testflag == "test":
            test = True
        else:
            test = False
    except:
        testflag = ''
        test = False
    run(target=target, band=band, order=order, suffix=f"_{testflag}_vz30", test=test)
