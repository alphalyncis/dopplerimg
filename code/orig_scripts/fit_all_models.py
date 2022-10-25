import numpy as np
from astropy.table import Table
import modelfitting as mf
from scipy import signal
import pickle
import matplotlib.pyplot as plt

# definitely replace with glob here
modelfns = ['lte010-4.0-0.0a+0.0.BT-Settl.spec.7.bz2','lte010-4.5-0.0a+0.0.BT-Settl.spec.7.bz2','lte010-5.0-0.0a+0.0.BT-Settl.spec.7.bz2','lte010-5.5-0.0a+0.0.BT-Settl.spec.7.bz2','lte010.5-4.0-0.0a+0.0.BT-Settl.spec.7.bz2','lte010.5-4.5-0.0a+0.0.BT-Settl.spec.7.bz2','lte010.5-5.0-0.0a+0.0.BT-Settl.spec.7.bz2','lte010.5-5.5-0.0a+0.0.BT-Settl.spec.7.bz2','lte011-4.0-0.0a+0.0.BT-Settl.spec.7.bz2','lte011-4.5-0.0a+0.0.BT-Settl.spec.7.bz2','lte011-5.0-0.0a+0.0.BT-Settl.spec.7.bz2','lte011-5.5-0.0a+0.0.BT-Settl.spec.7.bz2','lte011.5-4.0-0.0a+0.0.BT-Settl.spec.7.bz2','lte011.5-4.5-0.0a+0.0.BT-Settl.spec.7.bz2','lte011.5-5.0-0.0a+0.0.BT-Settl.spec.7.bz2','lte011.5-5.5-0.0a+0.0.BT-Settl.spec.7.bz2','lte012-4.0-0.0a+0.0.BT-Settl.spec.7.bz2','lte012-4.5-0.0a+0.0.BT-Settl.spec.7.bz2','lte012-5.0-0.0a+0.0.BT-Settl.spec.7.bz2','lte012-5.5-0.0a+0.0.BT-Settl.spec.7.bz2','lte012.5-4.0-0.0a+0.0.BT-Settl.spec.7.bz2','lte012.5-4.5-0.0a+0.0.BT-Settl.spec.7.bz2','lte012.5-5.0-0.0a+0.0.BT-Settl.spec.7.bz2','lte012.5-5.5-0.0a+0.0.BT-Settl.spec.7.bz2','lte013-4.0-0.0a+0.0.BT-Settl.spec.7.bz2','lte013-4.5-0.0a+0.0.BT-Settl.spec.7.bz2','lte013-5.0-0.0a+0.0.BT-Settl.spec.7.bz2','lte013-5.5-0.0a+0.0.BT-Settl.spec.7.bz2','lte013.5-4.0-0.0a+0.0.BT-Settl.spec.7.bz2','lte013.5-4.5-0.0a+0.0.BT-Settl.spec.7.bz2','lte013.5-5.0-0.0a+0.0.BT-Settl.spec.7.bz2','lte013.5-5.5-0.0a+0.0.BT-Settl.spec.7.bz2','lte014-4.0-0.0a+0.0.BT-Settl.spec.7.bz2','lte014-4.5-0.0a+0.0.BT-Settl.spec.7.bz2','lte014-5.0-0.0a+0.0.BT-Settl.spec.7.bz2','lte014-5.5-0.0a+0.0.BT-Settl.spec.7.bz2','lte014.5-4.0-0.0a+0.0.BT-Settl.spec.7.bz2','lte014.5-4.5-0.0a+0.0.BT-Settl.spec.7.bz2','lte014.5-5.0-0.0a+0.0.BT-Settl.spec.7.bz2','lte014.5-5.5-0.0a+0.0.BT-Settl.spec.7.bz2','lte015-4.0-0.0a+0.0.BT-Settl.spec.7.bz2','lte015-4.5-0.0a+0.0.BT-Settl.spec.7.bz2','lte015-5.0-0.0a+0.0.BT-Settl.spec.7.bz2','lte015-5.5-0.0a+0.0.BT-Settl.spec.7.bz2','lte015.5-4.0-0.0a+0.0.BT-Settl.spec.7.bz2','lte015.5-4.5-0.0a+0.0.BT-Settl.spec.7.bz2','lte015.5-5.0-0.0a+0.0.BT-Settl.spec.7.bz2','lte015.5-5.5-0.0a+0.0.BT-Settl.spec.7.bz2','lte016-4.0-0.0a+0.0.BT-Settl.spec.7.bz2','lte016-4.5-0.0a+0.0.BT-Settl.spec.7.bz2','lte016-5.0-0.0a+0.0.BT-Settl.spec.7.bz2','lte016-5.5-0.0a+0.0.BT-Settl.spec.7.bz2']


_mod = '/Users/ianc/proj/pcsa/data/model/'
modelfns = [_mod + modelfn.replace('.bz2', '.fits') for modelfn in np.sort(modelfns)]
allfits = []
allmods = []
alllams = []
allcors = []

for modelfn in modelfns:
  model = pyfits.getdata(modelfn)
  NPW = 4
  npix = wobs.shape[1]
  pix = np.arange(npix, dtype=float)/npix
  chipfits = []
  chipmods = np.zeros((4, 1024), dtype=float)
  chipcors = np.zeros((4, 1024), dtype=float)
  chiplams = np.zeros((4, 1024), dtype=float)

  for jj in range(4):
    lolim = wobs[jj].min() - 0.003
    hilim = wobs[jj].max() + 0.003
    tind = (model[0]>lolim) * (model[0] < hilim)
    aind = (atm0[0]>lolim) * (atm0[0] < hilim)
    lam_template = model[0,tind]
    template = model[1,tind]
    template /= np.median(template)
    lam_atmo = atm0[0,aind]
    atmo = atm0[1,aind]
    wcoef = np.polyfit(pix, wobs[jj], NPW-1)
    ccoef = [-0.1, 1.2/np.median(template)]
    NPC = len(ccoef)
    ind90 = np.sort(fobs0[jj])[int(0.9*npix)]
    ccoef = np.polyfit(pix[fobs0[jj]>ind90], fobs0[jj][fobs0[jj]>ind90], NPC-1)

    guess = np.concatenate(([21, 0.3, 9e-5, 42, 1.3], wcoef, ccoef))
    #mygmod, mygw = fit_atmo.modelspec_tel_template(guess, lam_template, template, lam_atmo, atmo, NPW, NPC, npix, retlam=True)
    fitargs = (fit_atmo.modelspec_tel_template, lam_template, template, lam_atmo, atmo, NPW, NPC, npix, fobs0[jj], wfobs0[jj])
    fit = an.fmin(pc.errfunc, guess, args=fitargs, full_output=True, disp=False)
    mymod, myw = fit_atmo.modelspec_tel_template(fit[0], *fitargs[1:-2], retlam=True)
    mycor = fit_atmo.modelspec_tel_template(fit[0], lam_template, np.ones(template.size), *fitargs[3:-2], retlam=False)
    chipfits.append(fit)
    chipmods[jj] = mymod
    chiplams[jj] = myw
    chipcors[jj] = mycor

  allfits.append(chipfits)
  allmods.append(chipmods)
  alllams.append(chiplams)
  allcors.append(chipcors)

a = dict(allfits=allfits, allmods=allmods, alllams=alllams, allcors=allcors, wrough=wobs, dats=datsR, datsi=datsRi, fdats=fobs0, modelfns=modelfns)
tools.savepickle(a, 'crires_combined_spectral_modeling_%s_NPW=%i_latest2.pickle' % (name, NPW) )

allparams = array([[af[0] for af in aff] for aff in allfits])

allchisq = np.array([[ff[1] for ff in af] for af in allfits])
totchisq = allchisq.sum(1)
bestind = (totchisq==totchisq.min()).nonzero()[0][0]