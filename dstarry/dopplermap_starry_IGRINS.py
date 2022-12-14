import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import pickle
import starry
import pymc3 as pm
import pymc3_ext as pmx
import theano.tensor as tt
from tqdm import tqdm
from scipy.signal import savgol_filter
import os
homedir = os.path.expanduser('~')

test = False
if test:
    testflag="testfixed_"
else:
    testflag=""
firstchip = 4
nobs, nchip, npix = 14, 2, 1848
chips = range(firstchip, firstchip+nchip)
#nchip, chips = 10, range(0,20,2) # even chips
#nchip, chips = 10, range(1, 20, 2) # odd chips
resultdir = f'starry_IGRINS/{testflag}order{firstchip}+{nchip}_nocors_'

# Load the dataset
with open(f'{homedir}/uoedrive/result/CIFIST/IGRINS_W1049B_K_lte015.0-5.0.pickle', "rb") as f:
    data = pickle.load(f, encoding="latin1")

# Timestamps (from Ian Crossfield)
t = np.array(
    [
        0.0,
        0.384384,
        0.76825872,
        1.15339656,
        1.5380364,
        1.92291888,
        2.30729232,
        2.69268336,
        3.07654752,
        3.46219464,
        3.8473428,
        4.23171624,
        4.61583456,
        4.99971048,
    ]
)

# Use the first epoch's wavelength array
lams = data["chiplams"][0]

# Interpolate onto that array
observed = np.empty((nobs, nchip, npix))
template = np.empty((nobs, nchip, npix))
broadened = np.empty((nobs, nchip, npix))
for k in range(nobs):
    for i, c in enumerate(chips):
        observed[k][i] = np.interp(
            lams[c],
            data["chiplams"][k][c],
            data["fobs0"][k][c] #/ data["chipcors"][k][c],
        )
        template[k][i] = np.interp(
            lams[c],
            data["chiplams"][k][c],
            data["chipmodnobroad"][k][c] #/ data["chipcors"][k][c],
        )
        if test:
            observed[k][i] = np.interp(
            lams[c],
            data["chiplams"][k][c],
            data["chipcors"][6][c],
            )
            template[k][i] = np.interp(
            lams[c],
            data["chiplams"][k][c],
            data["chipcors"][k][c],
            )

# Smooth the data and compute the median error from the MAD
smoothed = np.array(
    [savgol_filter(np.mean(observed[:, c], axis=0), 19, 3) for c in range(nchip)]
)
resid = observed - smoothed
error = 1.4826 * np.median(np.abs(resid - np.median(resid)))

# Clip outliers aggressively
level = 4
mask = np.abs(resid) < level * error
mask = np.min(mask, axis=0)

# Manually clip problematic regions
#mask[0][np.abs(lams[0] - 2.290) < 0.00015] = False
#mask[1][np.abs(lams[1] - 2.310) < 0.0001] = False
#mask[3][np.abs(lams[3] - 2.33845) < 0.0001] = False
#mask[3][np.abs(lams[3] - 2.340) < 0.0004] = False

# Get the final data arrays
pad = 100
wav = [None for n in range(nchip)]
wav0 = [None for n in range(nchip)]
flux = [None for n in range(nchip)]
mean_spectrum = [None for n in range(nchip)]
for c in range(nchip):
    wav[c] = lams[c][:][pad:-pad]
    flux[c] = observed[:, c][:, :][:, pad:-pad]
    wav0[c] = lams[c][:]
    mean_spectrum[c] = np.mean(template[:, c][:, :], axis=0) # obs-averaged nobroad modelspec

from astropy.io import fits
telluricfile = f"{homedir}/uoedrive/data/telluric/transdata_0,5-14_mic_hires.fits"
atm0 = fits.getdata(telluricfile)
lolim=2.29
hilim=2.36
aind = (atm0[:, 0]>lolim) * (atm0[:, 0] < hilim)
lam_atmo = atm0[aind, 0]
atmo = atm0[aind, 1]

plt.figure(figsize=(20,4))
k = 6
for i, c in enumerate(chips):
    plt.plot(lams[c], observed[k,i,:], color="tab:blue", label="observed")
    plt.plot(lams[c], template[k,i,:], color="tab:orange", linewidth=1,label="template")
    plt.plot(lams[c], resid[k,i,:], color="tab:green", linewidth=0.5, label="error")
    plt.plot(lam_atmo, atmo, color="k", linewidth=0.1, label="atmo")
    plt.plot(lams[c], data["chipcors"][k][c], color="r", linewidth=0.2,label="chipcors")
plt.legend(loc=2, bbox_to_anchor=(1,1))

# Set up a pymc3 model so we can optimize
with pm.Model() as model:

    # Dimensions
    ydeg = 10
    udeg = 1
    nc = 1
    nt = 14

    # Regularization
    flux_err = 0.02  # Eyeballed; doesn't matter too much since we're not doing posterior inference
    spectrum_sigma = 5e-2  # Hand-tuned; increase this to allow more departure from the template spectrum

    # Fixed
    vsini_max = 30000.0
    period = pm.math.constant(4.9)  # Crossfield et al. (2014)
    inc = pm.math.constant(70.0)  # Crossfield et al. (2014)

    # Free
    vsini = pm.Normal("vsini", mu=26100.0, sd=200)  # Crossfield et al. (2014)
    u1 = pm.Uniform("u1", 0.0, 1.0)  # Crossfield et al. (2014)

    # Deterministic
    veq = vsini / tt.sin(inc * np.pi / 180)
    theta = 360.0 * t / period

    # The surface map
    A = starry.DopplerMap(ydeg=ydeg).sht_matrix(smoothing=0.075)
    npix = A.shape[1]

    p = pm.Uniform("p", lower=0.0, upper=1.0, shape=(npix,))
    y = tt.dot(A, p)

    # The spectrum in each channel
    spectrum = [None for c in range(nchip)]
    for c in range(nchip):
        spectrum[c] = pm.Normal(
            f"spectrum{c}",
            mu=mean_spectrum[c], # the obs-avged chipmodnobroad is used as prior for spectrum
            sd=spectrum_sigma, 
            shape=mean_spectrum[c].shape,
        )

    # The continuum renormalization factor
    baseline = pm.Uniform("baseline", 0.3, 3.0, shape=(nt,), testval=0.65)

    # A small per-channel baseline offset
    offset = pm.Normal("offset", 0.0, 0.1, shape=(nchip,))

    # Compute the model & likelihood for each channel
    map = [None for c in range(nchip)]
    flux_model = [None for c in range(nchip)]
    for c in range(nchip):

        # Instantiate a map
        map[c] = starry.DopplerMap(
            ydeg=ydeg,
            udeg=udeg,
            nc=nc,
            veq=veq,
            inc=inc,
            nt=nt,
            wav=wav[c],
            wav0=wav0[c],
            lazy=True,
            vsini_max=vsini_max,
        )
        map[c][1] = u1
        map[c][:, :] = y
        map[c].spectrum = spectrum[c]

        # Compute the model
        flux_model[c] = offset[c] + tt.reshape(baseline, (-1, 1)) * map[
            c
        ].flux(theta=theta, normalize=False)

        # Likelihood term
        pm.Normal(
            f"obs{c}",
            mu=tt.reshape(flux_model[c], (-1,)),
            sd=flux_err,
            observed=flux[c].reshape(
                -1,
            ),
        )

# Optimization settings
lr = 1e-3
niter = 1000

# Optimize!
loss = []
best_loss = np.inf
map_soln = model.test_point
iterator = tqdm(
    pmx.optim.optimize_iterator(pmx.optim.Adam(lr=lr), niter, start=map_soln),
    total=niter,
    disable=os.getenv("CI", "false") == "true",
)
with model:
    for obj, point in iterator:
        iterator.set_description(
            "loss: {:.3e} / {:.3e}".format(obj, best_loss)
        )
        loss.append(obj)
        if obj < best_loss:
            best_loss = obj
            map_soln = point

# Plot the loss
fig, ax = plt.subplots(1, figsize=(5, 5))
ax.plot(loss)
ax.set_xlabel("iteration number")
ax.set_ylabel("loss")
#fig.savefig(f"{resultdir}luhman16b_loss.png", bbox_inches="tight")


# Plot the MAP map
with model:
    y_map = pmx.eval_in_model(y, point=map_soln)
    inc_map = pmx.eval_in_model(inc, point=map_soln)
map_map = starry.Map(ydeg, inc=inc_map)
map_map[:, :] = y_map

map_map.show(projection="moll", colorbar=True)
map_map.show(projection="moll", colorbar=True, file=f"{resultdir}luhman16b_map.png")

# Save the MAP map (just in case we need it later)
np.savez(f"{resultdir}luhman16b_map.npz", y_map=y_map, inc_map=inc_map)



'''
# To plot the other way:
oldmap = np.load('starry_IGRINS/order0+20_nocors_luhman16b_map.npz')
y_map = oldmap['y_map']
inc_map = oldmap['inc_map']
map_map = starry.Map(ydeg, inc=inc_map)
map_map[:, :] = y_map
times = np.array([0.0, 0.8, 1.6, 2.4, 3.2, 4.1])
thetas = 360 * times / period
fig = plt.figure(figsize=(8, 8))
f = 1 / 0.64
ax = [
    plt.axes([0.3225 * f, 0.34 * f, 0.3125, 0.3125]),
    plt.axes([0.44 * f, 0.17 * f, 0.3125, 0.3125]),
    plt.axes([0.3225 * f, 0.0, 0.3125, 0.3125]),
    plt.axes([0.1175 * f, 0.0, 0.3125, 0.3125]),
    plt.axes([0.0, 0.17 * f, 0.3125, 0.3125]),
    plt.axes([0.1175 * f, 0.34 * f, 0.3125, 0.3125]),
]
for n, axis in enumerate(ax):
    map_map.show(ax=axis, theta=thetas[n])
    axis.invert_yaxis()
    axis.invert_xaxis()
    angle = np.pi / 3 * (1 - n)
    plt.annotate(
        "{:.1f} hr".format(times[n]),
        xy=(0.515, 0.43),
        xycoords="figure fraction",
        xytext=(80 * np.cos(angle), 65 * np.sin(angle)),
        textcoords="offset points",
        va="center",
        ha="center",
        fontsize=15,
    )
'''