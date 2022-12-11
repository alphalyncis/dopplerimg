2013-05-05 16:39 IJMC: Beginning reduction of CRIRES data for Luhman
16 Brown Dwarf binary.  Trying to follow the steps in Section 8 of the
CRIRES pipeline user's manual.


rexfile = 'rexfile_cal'
darkname = 'crires_spec_dark_cal.fits'
darklist = 'input_dark_cal.txt'
flatname = 'crires_spec_flat_cal_set01.fits'  # must include set01!
bpmname = 'crires_spec_flat_cal_bpm_set01.fits'  # must include set01!
flatlist = 'input_flat_cal.txt'
calname = 'crires_util_combine_cal' # no ".fits" !!!
callist = 'input_cal.txt'

import pyfits
import os
from pylab import *
fits = []
lis = os.listdir('.')
for fn in lis:
    if fn[-5:]=='.fits' and (fn.find('CRIRE.2013-05-04T2')>-1 or fn.find('CRIRE.2013-05-05T0')>-1): fits.append(fn)

fits = np.array(fits)
hdrs = map(pyfits.getheader, fits)
object = array([h['object'] for h in hdrs])
filter = array([h['HIERARCH ESO INS FILT1 NAME'] for h in hdrs])
ndit = array([h['HIERARCH ESO DET NDIT'] for h in hdrs])
dit = array([h['HIERARCH ESO DET DIT'] for h in hdrs])
ny = array([h['HIERARCH ESO DET WINDOW NY'] for h in hdrs])
mjd0 = np.array([h['MJD-OBS'] for h in hdrs])
airmass0 = np.array([h['HIERARCH ESO TEL AIRM START'] for h in hdrs])
airmass1 = np.array([h['HIERARCH ESO TEL AIRM END'] for h in hdrs])
airmass = 0.5*(airmass0 + airmass1)
seeing = np.array([h['HIERARCH ESO TEL IA FWHM'] for h in hdrs])
humid  = np.array([h['HIERARCH ESO TEL AMBI RHUM'] for h in hdrs])
pa_start  = np.array([h['HIERARCH ESO TEL PARANG START'] for h in hdrs])
pa_end  = np.array([h['HIERARCH ESO TEL PARANG END'] for h in hdrs])
utc  = np.array([h['DATE-OBS'] for h in hdrs])
pa = 0.5 * (pa_start + pa_end)

#fmtstr = '%26s   %1.2f   %1.2f   %1.1f   %i'
#for ii in range(3,59):
#  print fmtstr % (utc[ii].replace('T', ' ')[0:-3], seeing[ii], airmass[ii], pa[ii], humid[ii])

exptime = ndit*dit

# manual says that DARK and FLAT must have the same DIT value:
dark4flat_ind = (object=='DARK') * (dit==2) * (ny==512)
flat_ind = (object=='FLAT') * (dit==2) * (exptime==60) #* (filter=='Ks') * (ny==512)
std_ind = (object=='STD') * (ny==512)
dark4flat_fns = fits[dark4flat_ind]
flat_fns = fits[flat_ind]
cal_fns = fits[std_ind]

# Write the DARKS file:
f = open(darklist, 'w')
for dfn in dark4flat_fns:
    f.write('%s  CAL_DARK\n' % dfn)

f.close()

# Write the FLATS file:
f = open(flatlist, 'w')
for dfn in flat_fns:
    f.write('%s  CAL_FLAT\n' % dfn)

f.write('%s CALPRO_DARK\n' % darkname)
f.close()

# Write the OBJECTS file:
f = open(callist, 'w')
for dfn in cal_fns:
    f.write('%s  OBS_NOD_JIT\n' % dfn)

f.write('%s CALPRO_PBM\n' % bpmname)
f.write('%s CALPRO_FLAT\n' % flatname)
f.close()

# Write the REX master control script:
f = open(rexfile, 'w')
f.write('esorex crires_spec_dark %s\n' % darklist)
f.write('mv crires_spec_dark.fits %s\n' % darkname)
f.write('esorex crires_spec_flat %s\n' % (flatlist))
f.write('mv crires_spec_flat_set01.fits %s\n' % flatname)
f.write('mv crires_spec_flat_set01_bpm.fits %s\n' % bpmname)
f.write('esorex crires_util_combine -onlyA -onlyB %s\n' % callist)
f.write('esorex crires_util_extract crires_util_combine_comb.fits crires_util_combine_contrib.fits\n')
f.write('mv crires_util_extract_extracted.fits %s_all.fits \n' % calname )
f.write('esorex crires_util_extract crires_util_combine_comb_noddedA.fits crires_util_combine_contrib_noddedA.fits \n')
f.write('mv crires_util_extract_extracted.fits %s_A.fits\n' % calname )
f.write('esorex crires_util_extract crires_util_combine_comb_noddedB.fits crires_util_combine_contrib_noddedB.fits \n')
f.write('mv crires_util_extract_extracted.fits %s_B.fits\n'  % calname)
f.close()



##################################################
# Calibrate the Science Data
##################################################

import pyfits
import os
from pylab import *

rexfile = 'rexfile_sci'
darkname = 'crires_spec_dark_sci.fits'
darklist = 'input_dark_sci.txt'
flatname = 'crires_spec_flat_cal_set01.fits'  # must include set01!
bpmname = 'crires_spec_flat_cal_bpm_set01.fits'  # must include set01!
sciname = 'crires_util_combine_sci' # no ".fits" !!!
scilist = 'input_sci.txt'


fits = []
lis = os.listdir('.')
for fn in lis:
    if fn[-5:]=='.fits' and fn[0:5]=='CRIRE': fits.append(fn)

fits = np.array(fits)
hdrs = map(pyfits.getheader, fits)
object = array([h['object'] for h in hdrs])
filter = array([h['HIERARCH ESO INS FILT1 NAME'] for h in hdrs])
ndit = array([h['HIERARCH ESO DET NDIT'] for h in hdrs])
dit = array([h['HIERARCH ESO DET DIT'] for h in hdrs])
ny = array([h['HIERARCH ESO DET WINDOW NY'] for h in hdrs])

exptime = ndit*dit

# manual says that DARK and FLAT must have the same DIT value:
dark4sci_ind = (object=='DARK') * (dit==300) * (ny==512)
sci_ind = (object=='WISE104915.57-531906.1') * (ny==512)
dark4flat_fns = fits[dark4flat_ind]
sci_fns = fits[sci_ind]

# Write the DARKS file:
f = open(darklist, 'w')
for dfn in dark4flat_fns:
    f.write('%s  CAL_DARK\n' % dfn)

f.close()

for kk in range(14):
  this_scilist = scilist.replace('.txt', '%i.txt'%kk)
  this_rexfile = rexfile + '_%i' % kk
  # Write the OBJECTS file:
  f = open(this_scilist, 'w')
  for dfn in sci_fns[9:][kk*4:(kk+1)*4]:
      f.write('%s  OBS_NOD_JIT\n' % dfn)
  f.write('%s CALPRO_PBM\n' % bpmname)
  f.write('%s CALPRO_FLAT\n' % flatname)
  f.close()
  # Write the REX master control script:
  f = open(this_rexfile, 'w')
  if kk==0:
    f.write('esorex crires_spec_dark %s\n' % darklist)
    f.write('mv crires_spec_dark.fits %s\n' % darkname)
  f.write('esorex crires_util_combine -onlyA -onlyB %s\n' % this_scilist)
  f.write('mv crires_util_combine_comb.fits crires_util_combine_comb_%i.fits\n' % kk)
  f.write('mv crires_util_combine_comb.fits crires_util_combine_comb_%i.fits\n' % kk)

  f.write('mv crires_util_combine_comb_noddedA.fits crires_util_combine_comb_noddedA_%i.fits\n' % kk)
  f.write('mv crires_util_combine_comb_noddedB.fits crires_util_combine_comb_noddedB_%i.fits\n' % kk)
  f.write('esorex crires_util_extract crires_util_combine_comb%i.fits crires_util_combine_contrib.fits\n' %kk )
  f.write('mv crires_util_extract_extracted.fits %s_%i_all.fits \n' % (sciname, kk) )
  f.write('esorex crires_util_extract crires_util_combine_comb_noddedA.fits crires_util_combine_contrib_noddedA.fits \n')
  f.write('esorex crires_util_extract crires_util_combine_comb_noddedA_%i.fits crires_util_combine_contrib_noddedA.fits \n' % kk )
  f.write('mv crires_util_extract_extracted.fits %s_%i_A.fits\n' % (sciname, kk) )
  f.write('esorex crires_util_extract crires_util_combine_comb_noddedB_%i.fits crires_util_combine_contrib_noddedB.fits \n' %kk )
  f.write('mv crires_util_extract_extracted.fits %s_%i_B.fits\n'  % (sciname, kk))
  f.close()

f = open('sciby4', 'w')
f.writelines(['./rexfile_sci_%i\n' % kk for kk in range(14)])
f.close()
os.system('chmod u+x rexfile_sci_*')
os.system('chmod u+x sciby4')
os.system('./sciby4')

##################################################
2013-05-07 09:16 IJMC: Another way of extracting spectra:

from pyraf import iraf
from pylab import *
from scipy import signal
import pyfits, os
import numpy as np
import nsdata as ns
import analysis as an

iraf.echelle.dispaxis=1
iraf.onedspec.dispaxis=1
iraf.apall.upper = 5
iraf.apall.lower = -5
bsampL = '-43:-23'
bsampR = '23:43'
itime = 300
file = pyfits.open('crires_util_combine_comb_all.fits')
nframes = 56

for kk in range(1,14):
  file = pyfits.open('crires_util_combine_comb_%i.fits' % kk )
  nframes = 4
  for ii in range(4):
      thisfn = 'chip%i.fits' % (ii+1)
      sdata = signal.medfilt2d(file[ii+1].data, 3)
      noiselevel = an.dumbconf(((sdata - file[ii+1].data)).ravel(), .683)[0]
      badpix = np.abs(sdata - file[ii+1].data) > 10*noiselevel
      cleandata = ns.bfixpix(file[ii+1].data, badpix, n=4, retdat=True)
      pyfits.writeto(thisfn, cleandata * itime * nframes, clobber=True) # convert to ADU

      iraf.apall(thisfn, output='specL_chip%i_frame%i' % (ii+1,kk), nfind=1, minsep=9, background='median', weights='variance', pfit='fit1d', readnoise=10, gain=7, t_nlost=999, nsum=-11, t_nsum=9, b_niterate=1, b_order=1, b_sample=bsampL, t_order=3, t_naverage=-5, t_niterate=100, t_grow=100, t_sample='10:1014', find=True, recenter='yes', trace='yes', edit='yes', clean=True, review=False, ylevel='INDEF', peak='yes', upper=5, lower=-5, resize=False)
      iraf.apall(thisfn, output='specR_chip%i_frame%i' % (ii+1,kk), nfind=1, minsep=9, background='median', weights='variance', pfit='fit1d', readnoise=10, gain=7, t_nlost=999, nsum=-11, t_nsum=9, b_niterate=1, b_order=1, b_sample=bsampR, t_order=3, t_naverage=-5, t_niterate=100, t_grow=100, t_sample='10:1014', clean=True, review=False, ylevel='INDEF', peak='yes', upper=5, lower=-5, find=True, recenter='yes', trace='yes', edit='yes', resize=False)
      os.remove(thisfn)
      print
  file.close()

##################################################



calfile = pyfits.open('crires_util_combine_cal_all.fits')
for jj in range(4):
    cal = calfile[jj+1].data['Extracted_OPT']
    pyfits.writeto('cal_chip%i.fits' % (jj+1), cal)

calfile.close()

telluric_list = '/Users/ianc/proj/pcsa/data/atmo/telluric_hr_K.dat'
for jj in range(1,4):
    fn = 'cal_chip%i.fits' % (jj+1)
    iraf.identify(fn, coordlist=telluric_list, ftype='absorption', fwidth='10', niterate=3, low=5, high=5, order=3, func='chebyshev')
    iraf.dispcor(fn, output='', linearize=False, log=False,flux=True)
    iraf.wspectext(fn, fn.replace('.fits', '.txt'), header=False)

w0 = np.loadtxt('temp.txt')
wobs = w0[:,0]/1e4


spectra = []
sspectra = []
for kk in range(1,14):
    file0 = pyfits.open('crires_util_combine_sci_%i_all.fits' % kk )
    spectra.append(np.concatenate([f00.data['Extracted_OPT'] for f00 in file0[1:]]))
    sspectra.append(signal.medfilt(spectra[-1], 3))
    file0.close()

spectra = np.array(spectra)
sspectra = np.array(sspectra)




mm = spectra/median(spectra,1).reshape(len(spectra),1)
mm0 = spectra[:,0:1024]/cal/median(spectra[:,0:1024]/cal,1).reshape(len(spectra),1)



modelfn = '/Users/ianc/proj/pcsa/data/model/lte014-5.0-0.0a+0.0.BT-Settl.spec.7.fits'
modelfn = '/Users/ianc/proj/pcsa/data/model/lte014-5.5-0.0a+0.0.BT-Settl.spec.7.fits'
modelfn = '/Users/ianc/proj/pcsa/data/model/lte014.5-5.0-0.0a+0.0.BT-Settl.spec.7.fits'
model = pyfits.getdata(modelfn)

mod = interp(wobs*(1.-30/3e5), model[0], model[1])
obs0 = median(mm,0)[0:1024]
obs = median(mm0,0)[0:1024]/cal[0:1024]
mod /= median(mod)
obs /= median(obs)

weights = np.ones(obs.shape)
weights[0:500] = 0.
weights[880:] = 0.
xkern = arange(-50, 50.)
kern = an.gaussian([1., 7, 0., 0], xkern)
kern /= kern.sum()
mod1 = np.convolve(mod, kern, 'same')

fobs = signal.medfilt(obs, 3)
wind = (wobs>2.293) * (wobs<2.3)
wfobs[True-wind] = 0
lin = np.linspace(0,1,fobs.size)
dc = np.ones(fobs.size)
coefs = an.lsq((mod1[wind], lin[wind], dc[wind]), fobs[wind], wfobs[wind])[0]
gaussmod = coefs[0]*mod1 + coefs[1]*lin + coefs[2]
himod    = coefs[0]*mod  + coefs[1]*lin + coefs[2]



cmod, dkern, back, chisq  = dia.dsa(himod[wind], obs[wind], 21, w=wfobs[wind])



##################################################
# FIT TELLURIC-AFFECTED SPECTRUM:
##################################################
import pyfits, fit_atmo
import tools
import analysis as an
from scipy import signal
import phasecurves as pc
import nsdata as ns
import os
import modelstar as ms

fnsR = ['specR_chip%i.fits' % ii for ii in range(1,5)]
fnsL = ['specL_chip%i.fits' % ii for ii in range(1,5)]
datsR = array(map(pyfits.getdata,fnsR)).squeeze()
datsL = array(map(pyfits.getdata,fnsL)).squeeze()

fnsRi = ([['specR_chip%i_frame%i.fits' % (ii,kk) for ii in range(1,5)] for kk in range(14)])
datsRi = array([map(pyfits.getdata,fnsRi[kk]) for kk in range(14)]).squeeze()
fnsLi = ([['specL_chip%i_frame%i.fits' % (ii,kk) for ii in range(1,5)] for kk in range(14)])
datsLi = array([map(pyfits.getdata,fnsLi[kk]) for kk in range(14)]).squeeze()



wobs = np.zeros((4, 1024), dtype=float)
for jj in range(4):
    wobs[jj] = np.loadtxt('cal_chip%i.txt' % (jj+1))[:,0]/1e4

obs0 = np.median(datsRi[:,:,0,:], axis=0)*len(fnsRi)
eobs0 = np.median(datsRi[:,:,3,:], axis=0)*np.sqrt(len(fnsRi))
name = 'brighter'

obs0 = np.median(datsLi[:,:,0,:], axis=0)*len(fnsLi)
eobs0 = np.median(datsLi[:,:,3,:], axis=0)*np.sqrt(len(fnsLi))
name = 'fainter'


fobs0 = np.vstack([signal.medfilt(obs0[jj], 3) for jj in range(4)])
eobs0 /= np.median(fobs0, 1).reshape(4,1)
fobs0 /= np.median(fobs0, 1).reshape(4,1)
wfobs0 = 1./eobs0**2
#wind = (wobs>2.288) * (wobs<2.35)
wind = ((wobs > wobs.min(1).reshape(4,1)+0.0005) * (wobs < wobs.max(1).reshape(4,1)-0.0005))
wfobs0[True-wind] = 0


telfn = '/Users/ianc/proj/pcsa/data/atmo/transdata_2,28-2,37_mic_hires.fits'
atm0 = pyfits.getdata(telfn)

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

#### FIGURES ####
figure()
ax1=subplot(121, position=[.3, .13, .35, .8])
semilogx(allchisq, range(len(allfits)), '-', linewidth=1)
semilogx(totchisq, range(len(allfits)), 'k-', linewidth=3.5)
yticks(range(len(allfits)), [fn[1] for fn in map(os.path.split, modelfns)], fontsize=8)
xlabel('Chi-squared', fontsize=15)
minorticks_on()
xlim([allchisq.min()/1.5, totchisq.max()*1.5])
grid()
title(name)
ylim([-0.5, len(modelfns)-1])

ax2=subplot(122, position=[.65, .13, .25, .8])
ns.imshow(allchisq)
colorbar()
ax2.set_position([.66, .13, .17, .8])
xlabel('\nChip Number', fontsize=15)
xticks([0,1,2,3], [1,2,3,4])
yticks(ax1.get_yticks(), [])
ylim(ax1.get_ylim())

teff, logg, mh = np.array([ms.btsettl_param(os.path.split(fn)[1]) for fn in modelfns]).T
teff0 = np.sort(np.unique(teff))
logg0 = np.sort(np.unique(logg))
scaled_chisq = allchisq.sum(1) 
scaled_chisq *= wind.sum()/scaled_chisq.min()
cmap = griddata(teff, logg, scaled_chisq, teff0, logg0)
vmap = griddata(teff, logg, allparams[:,:,0].mean(1), teff0, logg0)

figure()
contourf(teff0, logg0, allchisq.sum(1).reshape(13,4).T, 21)
title('log10(scaled chi-squared; %s component)' %name )
xlabel('T_eff [K]')
ylabel('log10 g [cgs]')
xticks(teff0)
yticks(logg0)
minorticks_on()
colorbar()


figure(tools.nextfig(), [16, 7.5])
ax1=subplot(121, position=[.1, .63, .85, .25])
plot(alllams[bestind].T, fobs0.T/allcors[bestind].T, 'k-', linewidth=2)
ylabel('Flux', fontsize=16)
ax1.text(.9, .85, 'Telluric-Corrected Calibrated Spectrum', fontsize=16, horizontalalignment='right', transform=ax1.transAxes)
ax2=subplot(121, position=[.1, .38, .85, .25])
plot(alllams[bestind].T, fobs0.T, 'k-', linewidth=2)
plot(alllams[bestind].T, allmods[bestind].T, 'r-', linewidth=2)
ylabel('Flux', fontsize=16)
ax2.text(.9, .85, 'Uncalibrated and Best-fit Model Spectra', fontsize=16, horizontalalignment='right', transform=ax2.transAxes)

ax2.text(.1, .15, 'Best-fit Model: '+os.path.split(modelfns[bestind])[1], fontsize=16, horizontalalignment='left', transform=ax2.transAxes)
ax3=subplot(121, position=[.1, .13, .85, .25])
plot(alllams[bestind].T, fobs0.T - allmods[bestind].T, 'k-', linewidth=2)
ylabel('Residuals', fontsize=16)
ax3.text(.9, .85, 'Best-fit Model Residuals', fontsize=16, horizontalalignment='right', transform=ax3.transAxes)
xlabel('Wavelength [$\mu m$]', fontsize=16)
yl = an.dumbconf((fobs0.T - allmods[bestind].T).ravel(), .683)[0]*6
ylim([-yl, yl])

[ax.set_xlim([alllams[bestind].min()-0.0005, alllams[bestind].max()+0.0005]) for ax in [ax1, ax2, ax3]]
[ax.set_ylim([0, max(ax.get_ylim())*1.1]) for ax in [ax1, ax2]]
[ax.set_xticklabels([]) for ax in [ax1, ax2]]
[ax.minorticks_on() for ax in [ax1, ax2, ax3]]


#
##################################################
##################################################
2013-05-10 15:55 IJMC: Now look at individual files (or at least,
co-added into 14 sets of 4 frames each):

import pyfits, fit_atmo
import tools
import analysis as an
from scipy import signal
import phasecurves as pc
import nsdata as ns
import os
from pylab import *
_mod = os.path.expanduser('~') + '/proj/pcsa/data/model/'
_dat = os.path.expanduser('~') + '/proj/pcsa/data/raw/20130505/'

fnsRi = ([[_dat+'specR_chip%i_frame%i.fits' % (ii,kk) for ii in range(1,5)] for kk in range(14)])
datsRi = array([map(pyfits.getdata,fnsRi[kk]) for kk in range(14)]).squeeze()
fnsLi = ([[_dat+'specL_chip%i_frame%i.fits' % (ii,kk) for ii in range(1,5)] for kk in range(14)])
datsLi = array([map(pyfits.getdata,fnsLi[kk]) for kk in range(14)]).squeeze()



dats = datsLi.copy()
modelfn = _mod + 'lte015-5.0-0.0a+0.0.BT-Settl.spec.7.fits' # 'lte014.5-5.0-0.0a+0.0.BT-Settl.spec.7.fits' #
name = 'fainter'

dats = datsRi.copy()
modelfn = _mod + 'lte015-5.0-0.0a+0.0.BT-Settl.spec.7.fits'
name = 'brighter'

wobs = np.zeros((4, 1024), dtype=float)
for jj in range(4):
    wobs[jj] = np.loadtxt('cal_chip%i.txt' % (jj+1))[:,0]/1e4

obs0 = dats[:,:,0,:].copy()
eobs0 = dats[:,:,3,:].copy()
efobs0 = dats[:,:,3,:].copy()
nobs = obs0.shape[0]
fobs0 = np.array([[signal.medfilt(obs0[kk,jj], 3) for jj in range(4)] for kk in range(nobs)])
efobs0 /= np.median(fobs0, 2).reshape(nobs,4,1)
fobs0 /= np.median(fobs0, 2).reshape(nobs,4,1)
wfobs0 = 1./efobs0**2
eobs0 /= np.median(obs0, 2).reshape(nobs,4,1)
obs0 /= np.median(obs0, 2).reshape(nobs,4,1)

# Clean up the observations:
dif = obs0 - obs0.mean(0)
obs1 = np.zeros(obs0.shape, dtype=float)
badobs1 = np.zeros(obs0.shape, dtype=bool)
for jj in range(4):
    badobs1[:,jj] = np.abs(dif[:,jj]) > 1.666*an.dumbconf(dif[:,jj].ravel(), .9973)[0]
    obs1[:,jj] = ns.bfixpix(obs0[:,jj], badobs1[:,jj], n=4, retdat=True)


wobs0 = 1./eobs0**2
#wind = (wobs>2.288) * (wobs<2.35)
wind = ((wobs > wobs.min(1).reshape(4,1)+0.0005) * (wobs < wobs.max(1).reshape(4,1)-0.0005))
wfobs0[:,True-wind] = 0
wobs0[:,True-wind] = 0


telfn = os.path.expanduser('~') + '/proj/pcsa/data/atmo/transdata_2,28-2,37_mic_hires.fits'
atm0 = pyfits.getdata(telfn)

model = pyfits.getdata(modelfn)
NPW = 4
npix = wobs.shape[1]
pix = np.arange(npix, dtype=float)/npix
individual_fits = []
chipmods = np.zeros((nobs, 4, 1024), dtype=float)
chipmodnobroad = np.zeros((nobs, 4, 1024), dtype=float)
chipcors = np.zeros((nobs, 4, 1024), dtype=float)
chiplams = np.zeros((nobs, 4, 1024), dtype=float)


for kk in range(nobs):
  chipfits = []
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
    ind90 = np.sort(obs1[kk,jj])[int(0.9*npix)]  
    ccoef = np.polyfit(pix[obs1[kk,jj]>ind90], obs1[kk,jj][obs1[kk,jj]>ind90], NPC-1)

    guess = np.concatenate(([25, 0.3, 9e-5, 42, 1.3], wcoef, ccoef))
    #mygmod, mygw = fit_atmo.modelspec_tel_template(guess, lam_template, template, lam_atmo, atmo, NPW, NPC, npix, retlam=True)
    fitargs = (fit_atmo.modelspec_tel_template, lam_template, template, lam_atmo, atmo, NPW, NPC, npix, obs1[kk,jj], wobs0[kk,jj])
    fit = an.fmin(pc.errfunc, guess, args=fitargs, full_output=True, disp=False)
    #testpos = tools.get_emcee_start(fit[0], np.abs(fit[0]/1000), 10, np.inf, fitargs, retchisq=False)
    #testfits = [an.fmin(pc.errfunc, testpos0, args=fitargs, full_output=True, disp=False) for testpos0 in testpos]
    #testchi = np.array([testfit0[1] for testfit0 in testfits])
    #if testchi.min() < fit[1]:
    #  fit = testfits[(testchi==testchi.min()).nonzero()[0][0]]
    
    mymod, myw = fit_atmo.modelspec_tel_template(fit[0], *fitargs[1:-2], retlam=True)
    mycor = fit_atmo.modelspec_tel_template(fit[0], lam_template, np.ones(template.size), *fitargs[3:-2], retlam=False)
    fitnobroad = fit[0].copy()
    fitnobroad[0] = 0.
    mymodnobroad = fit_atmo.modelspec_tel_template(fitnobroad, *fitargs[1:-2])
    chipfits.append(fit)
    chipmodnobroad[kk,jj] = mymodnobroad
    chipmods[kk,jj] = mymod
    chiplams[kk,jj] = myw
    chipcors[kk,jj] = mycor
  individual_fits.append(chipfits)


saveout = dict(chipmods=chipmods, chiplams=chiplams, chipcors=chipcors, obs1=obs1, obs0=obs0, wobs0=wobs0, wobs=wobs, modelfn=modelfn, individual_fits=individual_fits, chipmodnobroad=chipmodnobroad)
tools.savepickle(saveout, name+'spectral-fits_6.pickle')

#allfits = array([[ii[0] for ii in jj] for jj in individual_fits])
#_dat = os.path.expanduser('~') + '/proj/pcsa/data/raw/20130505/'
#dat = tools.loadpickle(_dat + 'fainterspectral-fits_6.pickle')
#nobs = 14
#for key in dat.keys():
#  exec('%s = dat["%s"]' % (key, key))

masterlam = chiplams.mean(0)
obs1_interp = np.array([[np.interp(masterlam[jj], chiplams[ii,jj], obs1[ii,jj], left=0, right=0) for jj in range(4)] for ii in range(nobs)])
chipcors_interp = np.array([[np.interp(masterlam[jj], chiplams[ii,jj], chipcors[ii,jj], left=0, right=0) for jj in range(4)] for ii in range(nobs)])
masterspec_all = obs1_interp/chipcors_interp
masterspec = masterspec_all.mean(0)
emasterspec = masterspec_all.std(0)/np.sqrt(nobs)
f = open('masterspec_%s.dat' % name, 'w')
for ii in xrange(masterspec.size):
  f.write('%1.7f   %1.6f   %1.6f\n' % (chiplams.mean(0).ravel()[ii], masterspec.ravel()[ii], emasterspec.ravel()[ii]))

f.close()


##################################################
##################################################
##################################################



import dia, dao, pyfits
from scipy import interpolate
from scipy import signal
import nsdata as ns
import analysis as an

_model = './'
pmod = 'output'
npix = wobs.shape[1]

#modelfn =  'lte013.5-4.5-0.0a+0.0.BT-Settl.spec.7.fits'
#modelfn =  '/Users/ianc/proj/pcsa/data/model/lte015-5.0-0.0a+0.0.BT-Settl.spec.7.fits'
#modelfn =  '/Users/ianc/proj/pcsa/data/model/lte015-5.5-0.0a+0.0.BT-Settl.spec.7.fits'
_mod = os.path.expanduser('~') + '/proj/pcsa/data/model/'
modelfn = modelfn.replace('/Users/ianc', os.path.expanduser('~'))
model = pyfits.getdata(modelfn)
tind = (model[0]>2.27) * (model[0] < 2.355)
lam_template = model[0,tind]
template = model[1,tind]
template /= np.median(template)

#pmod = 'lte015-5.0-0.0a+0.0.BT-Settl.spec.7_dao'
#pmod = 'lte013.5-4.5-0.0a+0.0.BT-Settl.spec.7_dao'
#pmod = os.path.split(modelfn)[1].replace('.fits', '_dao')

pmod = os.path.split(modelfn)[1].replace('BT-Settl.spec.7.fits', '_dao').replace('.',',')
#dao.runspec(lam_template*1e4, template, pmod, [lam_template.min()*1e4, lam_template.max()*1e4], 8, order=3)
print "You should really examine those automatically-measured line strengths by hand...!"



#f_linelist   = _mod + pmod + '.clineslsd'
f_linelist   = _mod + pmod + '_edited.clineslsd'
f_pspec      = _mod + pmod + '.fits'
f_pspec_cont = _mod + pmod + 'C.fits'


(lineloc, lineew, linespec) = dao.getlines(f_linelist)  
pspec_cont = pyfits.getdata(f_pspec_cont) 
hdr_pspec_cont = pyfits.getheader(f_pspec_cont) 
wspec = hdr_pspec_cont['crval1'] + np.arange(pspec_cont.size)*hdr_pspec_cont['cdelt1']
spline = interpolate.UnivariateSpline(wspec, pspec_cont, s=0.0, k=1.0)
deltaspec0 = ns.linespec(lineloc, lineew, wspec, verbose=False, cont=pspec_cont)



m0,k0,b0,c0 = dia.dsa(deltaspec0, template, 123)
dbeta0 = np.diff(wspec).mean()/wspec.mean()
dx0 = -np.arange(-k0.size/2, k0.size/2) * dbeta0
dv0 = an.c*dx0 / 1e3 # km/s


nk = 123
kerns = np.zeros((nobs, 4, nk), dtype=float)
modkerns = np.zeros((nobs, 4, nk), dtype=float)
### In this scheme, run LSD separately for each frame's wavelength solution:
for jj in range(4):
  for kk in range(14):
    deltaspec = ns.linespec(lineloc*(1.+9e-5), lineew, chiplams[kk,jj]*1e4, verbose=False, cont=spline(chiplams[kk,jj]*1e4))
    m,kerns[kk,jj],b,c = dia.dsa(deltaspec, obs1[kk,jj]/chipcors[kk,jj], nk)    
    m,modkerns[kk,jj],b,c = dia.dsa(deltaspec, chipmodnobroad[kk,jj]/chipcors[kk,jj], nk)    

dbeta = np.diff(chiplams).mean()/chiplams.mean()
dx = -dbeta * arange(np.floor(-nk/2.+.5), np.floor(nk/2.+.5))
dv = an.c*dx / 1e3 # km/s



### In this scheme, run LSD for all chips after interpolating to a
###    single wavelength grid:

hisamp = 3
dbeta = 2e-6
nk = 93 * hisamp
kerns = np.zeros((nobs, 4, nk), dtype=float)
modkerns = np.zeros((nobs, 4, nk), dtype=float)
wobs_hi = np.zeros((4, hisamp*npix), dtype=float)
deltaspec_hi = np.zeros((4, hisamp*npix), dtype=float)
obscor_hi = np.zeros((nobs, 4, hisamp*npix), dtype=float)

for jj in range(4):
    wobs_hi[jj] = np.arange(-hisamp*npix/2, hisamp*npix/2, dtype=float) * dbeta * chiplams[:,jj,:].mean() + chiplams[:,jj,:].mean()
    #wobs_hi[jj] = np.linspace(chiplams[:,jj,:].min(), chiplams[:,jj,:].max(), npix*hisamp)
    deltaspec_hi[jj] = ns.linespec(lineloc*(1.+9e-5), lineew, wobs_hi[jj]*1e4, verbose=False, cont=spline(wobs_hi[jj]*1e4))
    for kk in range(nobs):
        obscor_hi[kk,jj] = ns.interp(wobs_hi[jj], chiplams[kk,jj], obs1[kk,jj]/chipcors[kk,jj], left=0., right=0.)
        temp = ns.interp(wobs_hi[jj], chiplams[kk,jj], chipmodnobroad[kk,jj]/chipcors[kk,jj], left=0., right=0.)
        thisind = obscor_hi[kk,jj]>0 # avoid out-of-range values
        m,kerns[kk,jj],b,c = dia.dsa(deltaspec_hi[jj][thisind], obscor_hi[kk,jj][thisind], nk)
        mmm,modkerns[kk,jj],bbb,ccc = dia.dsa(deltaspec_hi[jj][thisind], temp[thisind], nk)
	

dx = -dbeta * arange(np.floor(-nk/2.+.5), np.floor(nk/2.+.5)) #np.arange(-nk/2, nk/2) * dbeta
dv = an.c*dx / 1e3 # km/s




kernstack = (kerns.mean(1) - kerns.mean(1).mean(0))
fkernstack = array([signal.medfilt(kk, 3) for kk in kernstack])
intrinsic_profile = modkerns.mean(0).mean(0)
intrinsic_profile -= intrinsic_profile[0] + (intrinsic_profile[-1] - intrinsic_profile[0])/(nk - 1.) * np.arange(nk)

systematic_rv_offset = (intrinsic_profile==intrinsic_profile.max()).nonzero()[0] - (dv==0).nonzero()[0]
intrinsic_profile = np.interp(np.arange(nk), np.arange(nk) - systematic_rv_offset, intrinsic_profile)
kerns = np.array([[np.interp(np.arange(nk), np.arange(nk) - systematic_rv_offset, kerns[jj, kk]) for kk in range(4)] for jj in range(nobs)])




mean_profiles = kerns.mean(0)
emean_profiles = kerns.std(0)
mean_profile = mean_profiles.mean(0)
#fft_profiles = np.array([[fft(kerns[jj,kk]) for kk in range(4)] for jj in range(nobs)])

# Fit for vsini using the LSD profile:
#rguess = [24, 0.6, 0.0]
#rot_model = spec.rotationalProfile(rguess, dv)
#gfwhm_pix = allfits[:,:,3].mean() * np.diff(lam_atmo).mean() / np.diff(chiplams,axis=2).mean()
#xkern = np.arange(int(-5*gfwhm_pix), int(5*gfwhm_pix))
#gkern = an.gaussian([1., gfwhm_pix/2.3548, 0, 0], xkern)
#conv_rot_model = np.convolve(rot_model, gkern, 'same')



from nsdata import imshow

xl = [-90, 90]
fs = 16

figure()
plot(dv, mean_profiles.T)
leg = legend(['Chip %i' % (ii+1) for ii in range(4)])
tools.legc(leg)
minorticks_on()
xlabel('Velocity [km/s]', fontsize=fs)
ylabel('Intensity', fontsize=fs)
title('Mean Line Profiles: %s component'  % name, fontsize=fs*1.2)
xlim(xl)


figure()
plot(dv, intrinsic_profile, 'k', linewidth=2)
tools.legc(leg)
minorticks_on()
xlabel('Velocity [km/s]', fontsize=fs)
ylabel('Intensity', fontsize=fs)
title('Intrinsic Model Profile: %s component'  % name, fontsize=fs*1.2)
xlim(xl)
text(0.5, 0.06, os.path.split(modelfn)[1], horizontalalignment='center', transform=gca().transAxes)


figure()
for ii in range(nobs):
    plot(dv, kerns[ii].mean(0) - ii/300., 'k')

minorticks_on()
xlim(xl)
xlabel('Velocity [km/s]', fontsize=fs)
ylabel('Intensity - Constant', fontsize=fs)
title('Time Series of Mean Line Profiles: %s component'  % name, fontsize=fs*1.2)

figure()
for ii in range(nobs):
    plot(dv, kernstack[ii] - ii/300., 'k')
    plot(xl, [0 - ii/300.]*2, '--k')

minorticks_on()
xlim(xl)
xlabel('Velocity [km/s]', fontsize=fs)
ylabel('Intensity - Constant', fontsize=fs)
title('Mean-Subtracted Time Series: %s component'  % name, fontsize=fs*1.2)




figure()
imshow(kernstack, x=dv, y=hr)
ylabel('Elapsed Time [hr]', fontsize=fs*1.3)
xlabel('Velocity [km/s]', fontsize=fs*1.3)
minorticks_on()
xlim(xl)
yl = ylim()
#ylim(ylim()[::-1])
#title('Mean-Subtracted Time Series: %s component'  % name, fontsize=fs*1.2)
plot([26.1]*2, ylim(), '--k', linewidth=3)
plot([-26.1]*2, ylim(), '--k', linewidth=3)
ylim(yl)
clim([-0.0025, .0025])
colorbar()


figure(tools.nextfig(), [10,10])
for ii in range(4):
    subplot(2,2,ii+1)
    imshow(kerns[:,ii]-kerns[:,ii].mean(0), x=dv)
    ylabel('Time', fontsize=fs)
    xlabel('Velocity [km/s]', fontsize=fs)
    minorticks_on()
    xlim(xl)
    title('Chip %i: %s component'  % (ii+1,name), fontsize=fs*1.2)
    clim([-0.004, .004])


tools.printfigs('/Users/ianc/python/figures/crires_%s_initial_plots4.pdf' % name)

2013-06-25 08:09 IJMC: Back to work!

The code in the several sections above suffices to adequately and
repeatably generate the mean LSD line profiles necessary for an
initial attempt at a Doppler Imaging analysis.  The next step is to
create the matrix framework necessary for this task, using the
approach of Vogt 1987.


2013-06-27 09:47 IJMC: 
# Add the timestamps into the data (I forgot earlier!):

for kk in range(14):
  file = pyfits.open('crires_util_combine_comb_%i.fits' % kk )
  for ii in range(4):
      datfn = 'specL_chip%i_frame%i.fits' % (ii+1, kk)
      data = pyfits.getdata(datfn)
      hdr = pyfits.getheader(datfn)
      hdr['MJD-OBS'] = file[0].header['MJD-OBS']
      hdr['DATE'] = file[0].header['DATE']
      pyfits.writeto(datfn, data, hdr, clobber=True)
  file.close()
  file = pyfits.open('crires_util_combine_comb_%i.fits' % kk )
  for ii in range(4):
      datfn = 'specR_chip%i_frame%i.fits' % (ii+1, kk)
      data = pyfits.getdata(datfn)
      hdr = pyfits.getheader(datfn)
      hdr['MJD-OBS'] = file[0].header['MJD-OBS']
      hdr['DATE'] = file[0].header['DATE']
      pyfits.writeto(datfn, data, hdr, clobber=True)
  file.close()

##################################################
# 2013-06-28 10:59 IJMC: 
##################################################
## Follow all the above steps to get the mean line profiles.  Then: to business!

import maps
import pyfits
import spec

home = os.path.expanduser('~')
#modelfn = home + '/proj/pcsa/data/model/lte015-5.5-0.0a+0.0.BT-Settl.spec.7.fits'

# Load the timestamps:
mjd = np.array([pyfits.getval(fn[0], 'MJD-OBS') for fn in fnsLi])
per = 4.87*3600 # seconds
phase = (mjd-mjd[0])*86400 * 2*np.pi/ per
vsini = 30e3 # m/s


# Initialize some important variables & constants:
inc = 0.4 
radius = an.rjup
nlat, nlon = 20, 40

# Load a model spectrum:
wlim = [2.15, 2.45]
lam0, fstar0 = pyfits.getdata(modelfn)
msinds = [find((lam0 - wlim[0]/1.002)>0)[0], find((lam0 - wlim[1]*1.002)<0)[-1]]
lam0 = lam0[msinds[0]:msinds[1]]
fstar0 = fstar0[msinds[0]:msinds[1]]
fstarSpline = interpolate.UnivariateSpline(lam0, fstar0, k=3., s=0.)

# Build a model of the intrinsic (unbroadened) line profile:
linelim0 = 2.3021274032
linelim = (2.30174868, 2.30250613)
linelim = (linelim0 * (1. + 1e3*dv.min()/an.c), linelim0 * (1. + 1e3*dv.max()/an.c))
lineind = (lam0 > linelim[0]) * (lam0 <= linelim[1])
lamdv = (1. + np.sort(dv)*1e3/an.c) * linelim0

deltafcn = 1. - np.convolve((np.abs(lam0-linelim0)==np.abs(lam0-linelim0).min()), [0.2, 0.2, 0.2, 0.2, 0.2], 'same')
gaussKern = an.gaussian([-1, 5.8, 0, 1], (lam0-linelim0)/np.diff(lam0[lineind]).mean())
flineSpline = interpolate.UnivariateSpline(lam0, deltafcn, k=1., s=0.)
flineSpline = fstarSpline
flineSpline = interpolate.UnivariateSpline(lam0, gaussKern, k=1., s=0.)
flineSpline2 = interpolate.UnivariateSpline(dv[::-1], intrinsic_profile[::-1], k=1., s=0.) 

modIP = 1. - np.concatenate((np.zeros(300), intrinsic_profile, np.zeros(300)))
modDV = - np.arange(np.floor(-modIP.size/2.+.5), np.floor(modIP.size/2.+.5)) * dbeta * an.c / 1e3
flineSpline2 = interpolate.UnivariateSpline(modDV[::-1], modIP[::-1], k=1., s=0.) 


# Define the 'map object':
map = maps.map(nlat=nlat, nlon=nlon, i=inc)
ncell = map.ncell
latlon_corners = np.array([c.corners_latlon for c in map.cells])
phi_corner = latlon_corners[:,1,:]
theta_corner = latlon_corners[:,0,:]

LLD = 1.0
# Compute the "Response Matrix" of Vogt 1987:
Rmatrix = np.zeros((ncell, nobs*dv.size), dtype=float32)
for kk, rot in enumerate(phase):
    speccube = np.zeros((ncell, dv.size), dtype=float32)
    this_map = maps.map(nlat=nlat, nlon=nlon, i=inc, deltaphi=-rot)
    this_doppler = 1. + vsini*this_map.visible_rvcorners.mean(1)/an.c/np.cos(inc)
    good = (this_map.projected_area>0) * np.isfinite(this_doppler)
    for ii in good.nonzero()[0]:
        speccube[ii,:] = flineSpline2(dv + (this_doppler[ii]-1)*an.c/1000.)
        #speccube[ii,:] = flineSpline(lamdv / this_doppler[ii])
    limbdarkening = (1. - LLD) + LLD * this_map.mu
    Rblock = speccube * ((limbdarkening*this_map.projected_area).reshape(this_map.ncell, 1)*np.pi/this_map.projected_area.sum())
    Rmatrix[:,dv.size*kk:dv.size*(kk+1)] = Rblock

# Add some weighted regularization: 
#  In this case, we add the constraint that each cell is roughly equal
#  to theneighboring cells directly (a) West and (b) South of it.
R0 = Rmatrix.copy()
for ii in range(map.ncell):
    ind1 = (phi_corner%(2*np.pi))==phi_corner[ii].min()
    ind2 = theta_corner==theta_corner[ii]
    ii_latadjacent = list(set((ind1*ind2).sum(1).nonzero()[0]) - set([ii]))[0]
    newvec = np.zeros(map.ncell) 
    newvec[ii] = 1
    newvec[ii_latadjacent] = -1
    ind3 = (phi_corner%(2*np.pi))==phi_corner[ii]
    ind4 = theta_corner.max(1)==theta_corner[ii].min()
    ii_lonadjacent = (ind3.all(1)*ind4).nonzero()[0]
    R0 = np.hstack((R0, newvec.reshape(ncell,1) ))
    if ii_lonadjacent.size>0:
        newvec2 = np.zeros(map.ncell) # cells at same longitude
        newvec2[ii] = 1
        newvec2[ii_lonadjacent] = -1
        R0 = np.hstack((R0, newvec2.reshape(ncell,1) ))
            
observation = 1. - kerns.mean(1).ravel()
observation_constraint = observation.copy()
observation_constraint = np.concatenate((observation_constraint, np.zeros(R0.shape[1] - Rmatrix.shape[1])))

# Continuum-normalize each observed line profile:
observation_cnorm = observation_constraint.copy()
for ii in range(nobs):
    i0, i1 = ii*nk, (ii+1)*nk
    inds = np.concatenate((np.arange(i0, i0+7), np.arange(i1-7, i1)))
    continuumfit = np.polyfit(inds, observation_cnorm[inds], 1)
    observation_cnorm[i0:i1] /= np.polyval(continuumfit, np.arange(i0, i1))

observation_norm = observation.copy()
for ii in range(nobs):
    i0, i1 = ii*nk, (ii+1)*nk
    inds = np.concatenate((np.arange(i0, i0+7), np.arange(i1-7, i1)))
    continuumfit = np.polyfit(inds, observation_norm[inds], 1)
    observation_norm[i0:i1] /= np.polyval(continuumfit, np.arange(i0, i1))




# Prepare the weighting matrix & other arrays for weighted LSQ:
e_observation = kerns.mean(1).std(0)
e_observation = np.tile(np.median(kerns.mean(1).std(0)), nk)
w_observation = 1./np.tile(e_observation, nobs)**2

w_constraints = 0.001*median(w_observation) * np.ones(R0.shape[1] - Rmatrix.shape[1]) 
W_matrix = np.diag(np.concatenate((w_observation, w_constraints)))
R0W = np.dot(R0, W_matrix)
R0WRTinv = np.linalg.inv(np.dot(R0W, R0.T))

RW = np.dot(Rmatrix, np.diag(w_observation))
RWRTinv = np.linalg.pinv(np.dot(RW, Rmatrix.T))
rm = np.dot(RWRTinv, np.dot(RW, observation_norm)).reshape(nlat, nlon)
mo = np.dot(rm.ravel(), Rmatrix)


# Invert the line profiles into maps:
recovered_map0 = np.dot(R0WRTinv, np.dot(R0W, observation_constraint)).reshape(nlat, nlon)
recovered_map = np.dot(R0WRTinv, np.dot(R0W, observation_cnorm)).reshape(nlat, nlon)


model_observation0 = np.dot(recovered_map0.ravel(), R0)
model_observation = np.dot(recovered_map.ravel(), R0)

chisq_data = ((observation_cnorm - model_observation)[0:Rmatrix.shape[1]]**2 * w_observation).sum()
chisq_constraints = ((observation_cnorm - model_observation)[Rmatrix.shape[1]:]**2 * w_constraints).sum()
entropic_map = recovered_map.copy()
entropic_map[entropic_map<=0] = 1e-3
map_entropy = -(entropic_map/entropic_map.sum() * np.log(entropic_map/entropic_map.sum())).sum()

#tools.plotc(map.visible_corners[:,1,:].mean(1), map.visible_corners[:,2,:].mean(1), 1-recovered_map.ravel().astype(float), zmin=-25, zmax=25)

### GOTO (2013-08-06 17:28) below.



##################################################
# 2013-07-01 09:22 IJMC: Start to think about maximum-entropy considerations:
##################################################
import phasecurves as pc
from mpl_toolkits.basemap import Basemap

e_observation = np.ones(nk)*np.median(kerns.mean(1).std(0))
e_observation[0:80] = 1e9
e_observation[-80:] = 1e9
w_observation = 1./np.tile(e_observation, nobs)**2

RW = np.dot(Rmatrix, np.diag(w_observation))
RWRTinv = np.linalg.pinv(np.dot(RW, Rmatrix.T))
recovered_map = np.dot(RWRTinv, np.dot(RW, observation_cnorm[0:nk*nobs])).reshape(nlat, nlon)
model_observation = np.dot(recovered_map.ravel(), Rmatrix)
# Re-scale the weights to give a more 'reasonable' chi-squared value (i.e., = N)
chisq_firsttry = (w_observation*(observation_cnorm[0:nk*nobs] - model_observation)**2).sum()
npts = (w_observation>1e-10).sum()
w_observation *= 1.0 * npts / chisq_firsttry

j = np.arange(Rmatrix.shape[1])
k = (np.floor(1.0*j/nk)*nk).astype(int)
l = (np.ceil((j+1.0)/nk)*nk - 1.).astype(int)
jfrac = (j % nk) / (nk - 1.0)

def gnorm(unnorm_model, nk):
  return unnorm_model[k] + (unnorm_model[l] - unnorm_model[k]) * jfrac

def normalize_model(unnorm_model, nk):
  return unnorm_model / gnorm(unnorm_model, nk)

###  F = 0.5*chi^2 - alpha*S (from Bridle et al. 1998)
###    Set alpha such that F_min = 1.5*npts
def entropy_map(map_pixels, data, weights, R, alpha=1, minval=1e-3, maxlim=10):
  map_pixels[map_pixels <= 0] = minval
  norm_pixels = map_pixels / np.sum(map_pixels)
  entropy = np.sum(norm_pixels * np.log(norm_pixels))
  model = np.dot(map_pixels.ravel(), R)
  chisq = pc._chi2.chi2(model, data, weights)
  if (map_pixels > maxlim).any():  chisq *= 1e4
  return 0.5*chisq - alpha*entropy

alpha=400
thisfit = [guess]
fits = []

guess = np.ones((nlat*nlon), dtype=float32)/np.pi
entropy_map(guess, observation_cnorm[0:nk*nobs], w_observation, Rmatrix, alpha=alpha)

while True:
   thisfit = an.fmin(entropy_map, thisfit[0], args=(observation_cnorm[0:nk*nobs], w_observation, Rmatrix), kw=dict(alpha=alpha), full_output=True, maxiter=200000, maxfun=200000)
   fits.append(thisfit)
   print thisfit[1], thisfit[0].min(), thisfit[0].max()

model_observation = np.dot(fits[-1][0], Rmatrix)
chisq_2 = (w_observation*(observation_cnorm[0:nk*nobs] - model_observation)**2).sum()
w_observation *= 1.0 * npts / chisq_2
this_entropy = entropy_map(fits[-1][0], observation_cnorm[0:nk*nobs], w_observation, Rmatrix, alpha=alpha)
### weights and alpha until chisq=npts, and this_entropy=1.5*npts


sav = dict(map=map, bestfit=fits[-1], observation_cnorm=observation_cnorm, nk=nk, nobs=nobs, w_observation=w_observation, kerns=kerns, Rmatrix=Rmatrix, mjd=mjd, per=per, alpha=alpha, inc=inc, radius=radius,vsini=vsini)
savfn = ('bdmap_crudeinput_inc=%1.1f_r=%1.1frj_alpha=%1.5e__nlat=%i_nlon=%i' % \
      (inc, radius/an.rjup, alpha, nlat, nlon)).replace('.',',') + '.pickle'
tools.savepickle(sav, savfn)


latlon_corners = np.array([c.corners_latlon for c in map.cells])
phi_corner = latlon_corners[:,1,:]
theta_corner = latlon_corners[:,0,:]

hiphi = np.linspace(0, 2*np.pi, nlat*5)
hitheta = np.linspace(0, np.pi, nlon*5)
himap = griddata(np.concatenate((phi_corner-2*np.pi, phi_corner, phi_corner+2*np.pi)).ravel(), np.concatenate((theta_corner, theta_corner, theta_corner)).ravel(), np.concatenate((np.tile(fits[-1][0], (4,1)).T, np.tile(fits[-1][0], (4,1)).T, np.tile(fits[-1][0], (4,1)).T)).ravel(), hiphi, hitheta)


figure(tools.nextfig(), [18, 7])
for ii in range(7):
  subplot(2, 7, ii+1)
  jj = ii*2
  hiphi2, hitheta2 = np.meshgrid(hiphi, hitheta)
  pmap = Basemap(projection='ortho', lat_0=-inc*180./np.pi, lon_0=phase[jj]*180/np.pi)
  x1, y1 = pmap(360. - hiphi2*180./np.pi, -hitheta2*180./np.pi +  90)
  CS = pmap.contourf(x1, y1, himap, np.linspace(fits[-1][0].min(),fits[-1][0].max(), 21)) #, cmap=cm.YlOrBr)
  pmap.drawmeridians(np.arange(0, 360, 30)); 
  pmap.drawparallels(np.arange(-90, 90, 30))
  title('%1.1f hr' % (phase[jj]*per/(3600.*2*np.pi)))
  subplot(2,7,ii+8)
  plot(dv, observation_cnorm[jj*nk:(jj+1)*nk], 'k')
  plot(dv, model_observation[jj*nk:(jj+1)*nk], 'r')
  xlabel('RV [km/s]', fontsize=16)
  if jj==0: ylabel('Intensity', fontsize=16)
  ylim([.958, 1.002]); minorticks_on()




testmap = np.ones(map.ncell)
testmap[map.ncell/2+nlon*3] = 0.
testmap[map.ncell/2+nlon*4-1] = 0.
hitestmap = griddata(np.concatenate((phi_corner-2*np.pi, phi_corner, phi_corner+2*np.pi)).ravel(), np.concatenate((theta_corner, theta_corner, theta_corner)).ravel(), np.concatenate((np.tile(testmap, (4,1)).T, np.tile(testmap, (4,1)).T, np.tile(testmap, (4,1)).T)).ravel(), hiphi, hitheta)
testmodel = np.dot(testmap, Rmatrix)

figure()
for ii in range(14):
  subplot(2, 14, ii+1)
  jj = ii
  hiphi2, hitheta2 = np.meshgrid(hiphi, hitheta)
  pmap = Basemap(projection='ortho', lat_0=-inc*180./np.pi, lon_0=phase[jj]*180/np.pi)
  x1, y1 = pmap(360. - hiphi2*180./np.pi, -hitheta2*180./np.pi +  90)
  CS = pmap.contourf(x1, y1, hitestmap, np.linspace(0,1, 21), cmap=cm.jet)
  pmap.drawmeridians(np.arange(0, 360, 30)); 
  pmap.drawparallels(np.arange(-90, 90, 30))
  title('%1.1f hr' % (phase[jj]*per/(3600.*2*np.pi)))
  subplot(2,14,ii+15)
  plot(dv, testmodel[jj*nk:(jj+1)*nk] - np.median(testmodel.reshape(14, nk),0), 'r')
  plot([0,0], ylim(), '--k')
  xlabel('RV [km/s]', fontsize=16)
  minorticks_on()


##################################################
2013-08-06 17:28 IJMC: 
##################################################

from scipy import optimize
import dime  # Doppler Imaging & Maximum Entropy
import phasecurves as pc
dime.setup(observation.size, nk)



flatguess = 10*np.ones(nlat*nlon)
flatmodel = dime.normalize_model(np.dot(flatguess, Rmatrix), nk)
minfm = flatmodel.min()
cutoffval = 1. - (1. - minfm) / 22.
w_observation = (flatmodel < cutoffval).astype(float) / np.median(kerns.mean(1).std(0))**2


#### Try some line-profile vsini &  LLD fitting:

norm_meanprofile = mean_profile - (mean_profile[0] + (mean_profile[-1] - mean_profile[0])/(nk - 1.) * np.arange(nk))
rguess = [28, 0.5, 1., 1.0]
test = spec.modelline(rguess, intrinsic_profile, dv)
up = [(0, 200), (0, 3), (-100, 100), (0, 100)]
mcmcargs = (spec.modelline, intrinsic_profile, dv, norm_meanprofile, w_observation[0:nk], dict(uniformprior=up))
fit1 = an.fmin(pc.errfunc, rguess, args=mcmcargs)
test1 = spec.modelline(fit1, intrinsic_profile, dv)


import emcee
nwalkers = 100
ndim = len(rguess)
sampler = emcee.EnsembleSampler(nwalkers, ndim, pc.lnprobfunc, args=mcmcargs, threads=1)
p0 = [np.random.normal(fit1, np.abs(fit1)/1000.) for i in xrange(nwalkers)]
p1, prob, state = sampler.run_mcmc(p0, 500) # Burn-in
sampler.reset()
p2, prob, state = sampler.run_mcmc(p1, 1000) # Burn-in
bestparams = sampler.flatchain[nonzero(sampler.lnprobability.ravel()==sampler.lnprobability.ravel().max())[0][0]]
bestparams = an.fmin(pc.errfunc, bestparams, args=mcmcargs)
test3 = spec.modelline(bestparams, intrinsic_profile, dv)

## NOW GO AND PUT THIS vsini AND lld AND COMPUTE A NEW RMATRIX
################




# Scale the observations to match the model's equivalent width:
out, eout = an.lsq((observation_norm, np.ones(nobs*nk)), flatmodel, w=w_observation)
sc_observation_norm = observation_norm * out[0] + out[1]

alpha = 4000
nx = nlat*nlon
fitargs = (sc_observation_norm, w_observation, Rmatrix, alpha)



x0 = flatguess.copy()
s0 = dime.entropy_map_norm_sp(x0, *fitargs)
grad0 = dime.getgrad_norm_sp(x0, *fitargs)
x1 = x0 - 0.0001*grad0
s1 = dime.entropy_map_norm_sp(x1, *fitargs)
grad1 = dime.getgrad_norm_sp(x1, *fitargs)
x2 = x1 - 0.0001*grad1
s2 = dime.entropy_map_norm_sp(x2, *fitargs)
grad2 = dime.getgrad_norm_sp(x2, *fitargs)
x3 = x2 - 0.0001*grad2
s3 = dime.entropy_map_norm_sp(x3, *fitargs)

xx = flatguess.copy()
for ii in range(50):
  ss = dime.entropy_map_norm_sp(xx, *fitargs)
  gradgrad = dime.getgrad_norm_sp(xx, *fitargs)
  xx = xx - 0.001*gradgrad


bestparams = xx.copy()
factor = 1.
ent = 9e99
dent = 9e9

maxiter = 2000
iter = 0
recalcGrad = True
while iter<maxiter and np.abs(dent)>0.0003:
  if iter==0: ent = dime.entropy_map_norm_sp(bestparams, *fitargs)
  iter += 1
  if recalcGrad: grad = dime.getgrad_norm_sp(bestparams, *fitargs)
  newparam = bestparams - factor * grad
  newent = dime.entropy_map_norm_sp(newparam, *fitargs)
  if newent < ent:
    dent = newent - ent
    ent = newent
    if recalcGrad is True: factor *= 1.5  # we updated twice in a row!
    bestparams = newparam.copy()
    recalcGrad = True
  else:
    factor /= 2
    recalcGrad = False
    print iter, ent, factor, dent
  

pp = dime.normalize_model(np.dot(bestparams, Rmatrix), nk)
mm = dime.normalize_model(np.dot(xx, Rmatrix), nk)
m0 = dime.normalize_model(np.dot(x0, Rmatrix), nk)
m1 = dime.normalize_model(np.dot(x1, Rmatrix), nk)
m2 = dime.normalize_model(np.dot(x2, Rmatrix), nk)
m3 = dime.normalize_model(np.dot(x3, Rmatrix), nk)


sfit1 = optimize.fmin_tnc(dime.entropy_map_norm_sp, flatguess, fprime=dime.getgrad_norm_sp, args=fitargs, maxfun=10000, bounds=[(1e-6, 20)]*nx, ftol=1)
bfit1 = optimize.fmin_l_bfgs_b(dime.entropy_map_norm_sp, flatguess, fprime=dime.getgrad_norm_sp, args=fitargs, maxfun=10000, bounds=[(1e-6, 20)]*nx, factr=1e12, iprint=5)
cfit1 = optimize.fmin_slsqp(dime.entropy_map_norm_sp, flatguess, fprime=dime.getgrad_norm_sp, args=fitargs, iter=100000, bounds=[(1e-6, 20)]*nx, acc=1, full_output=True, iprint=1)


ss  = dime.normalize_model(np.dot(sfit1[0], Rmatrix), nk)
bb  = dime.normalize_model(np.dot(bfit1[0], Rmatrix), nk)
cc  = dime.normalize_model(np.dot(cfit1[0], Rmatrix), nk)

junk = np.concatenate((sfit1[0], bfit1[0], cfit1[0], bestparams))
cl = [junk.min(), junk.max()]

figure()
subplot(2,2,1)
imshow(bestparams.reshape(nlat, nlon)); clim(cl)
title('bestparams')
subplot(2,2,2)
imshow(sfit1[0].reshape(nlat, nlon)); clim(cl)
title('sfit')
subplot(2,2,3)
imshow(bfit1[0].reshape(nlat, nlon)); clim(cl)
title('bfit')
subplot(2,2,4)
imshow(np.array(cfit1[0]).reshape(nlat, nlon)); clim(cl)
title('cfit')

figure()
subplot(2,2,1)
imshow(1. - (mm.reshape(14,nk) - pp.reshape(14,nk).mean(0)))
title('bestparams')
subplot(2,2,2)
imshow(1. - (ss.reshape(14,nk) - ss.reshape(14,nk).mean(0)))
title('sfit')
subplot(2,2,3)
imshow(1. - (bb.reshape(14,nk) - bb.reshape(14,nk).mean(0)))
title('bfit')
subplot(2,2,4)
imshow(1. - (cc.reshape(14,nk) - cc.reshape(14,nk).mean(0)))
title('cfit')






##################################################
2013-08-08 16:13 IJMC: 
##################################################

from pylab import *
from mpl_toolkits.basemap import Basemap
import os
import tools
import dia, dao
from scipy import interpolate
from scipy import signal
import nsdata as ns
import analysis as an
import maps
import pyfits
import spec
from scipy import optimize
import dime  # Doppler Imaging & Maximum Entropy
import phasecurves as pc


# LOad the corrected spectra:
home = os.path.expanduser('~')
_dat = home + '/proj/pcsa/data/raw/20130505/'
dat = tools.loadpickle(_dat + 'fainterspectral-fits_6.pickle')
nobs = 14
for key in dat.keys():
  exec('%s = dat["%s"]' % (key, key))

# Set up for Least Squares Deconvolution
_model = './'
pmod = 'output'
npix = wobs.shape[1]
_mod = home + '/proj/pcsa/data/model/'
modelfn = modelfn.replace('/Users/ianc', os.path.expanduser('~'))
model = pyfits.getdata(modelfn)
tind = (model[0]>2.27) * (model[0] < 2.355)
lam_template = model[0,tind]
template = model[1,tind]
template /= np.median(template)
pmod = os.path.split(modelfn)[1].replace('BT-Settl.spec.7.fits', '_dao').replace('.',',')

f_linelist   = _mod + pmod + '_edited.clineslsd'
f_pspec      = _mod + pmod + '.fits'
f_pspec_cont = _mod + pmod + 'C.fits'

# Load the LSD files:
(lineloc, lineew, linespec) = dao.getlines(f_linelist)  
pspec_cont = pyfits.getdata(f_pspec_cont) 
hdr_pspec_cont = pyfits.getheader(f_pspec_cont) 
wspec = hdr_pspec_cont['crval1'] + np.arange(pspec_cont.size)*hdr_pspec_cont['cdelt1']
spline = interpolate.UnivariateSpline(wspec, pspec_cont, s=0.0, k=1.0)
deltaspec0 = ns.linespec(lineloc, lineew, wspec, verbose=False, cont=pspec_cont)

# Compute LSD:
nk = 123
kerns = np.zeros((nobs, 4, nk), dtype=float)
modkerns = np.zeros((nobs, 4, nk), dtype=float)
### In this scheme, run LSD separately for each frame's wavelength solution:
for jj in range(4):
  for kk in range(14):
    deltaspec = ns.linespec(lineloc*(1.+9e-5), lineew, chiplams[kk,jj]*1e4, verbose=False, cont=spline(chiplams[kk,jj]*1e4))
    m,kerns[kk,jj],b,c = dia.dsa(deltaspec, obs1[kk,jj]/chipcors[kk,jj], nk)    
    m,modkerns[kk,jj],b,c = dia.dsa(deltaspec, chipmodnobroad[kk,jj]/chipcors[kk,jj], nk)    

# Compute LSD velocity grid:
dbeta = np.diff(chiplams).mean()/chiplams.mean()
dx = -dbeta * arange(np.floor(-nk/2.+.5), np.floor(nk/2.+.5))
dv = an.c*dx / 1e3 # km/s

kernstack = (kerns.mean(1) - kerns.mean(1).mean(0))
fkernstack = array([signal.medfilt(kk, 3) for kk in kernstack])
intrinsic_profile = modkerns.mean(0).mean(0)
intrinsic_profile -= intrinsic_profile[0] + (intrinsic_profile[-1] - intrinsic_profile[0])/(nk - 1.) * np.arange(nk)

systematic_rv_offset = (intrinsic_profile==intrinsic_profile.max()).nonzero()[0] - (dv==0).nonzero()[0]
intrinsic_profile = np.interp(np.arange(nk), np.arange(nk) - systematic_rv_offset, intrinsic_profile)
kerns = np.array([[np.interp(np.arange(nk), np.arange(nk) - systematic_rv_offset, kerns[jj, kk]) for kk in range(4)] for jj in range(nobs)])


# Prepare the Doppler Imaging "response matrix"

# Load the timestamps:
fnsLi = ([['specL_chip%i_frame%i.fits' % (ii,kk) for ii in range(1,5)] for kk in range(14)])
mjd = np.array([pyfits.getval(fn[0], 'MJD-OBS') for fn in fnsLi])
#airmass = np.array([pyfits.getval(fn[0], 'HIERARCH ESO TEL AIRM END') for fn in fnsLi])

per = 4.87*3600 # seconds
phase = (mjd-mjd[0])*86400 * 2*np.pi/ per
vsini = 30e3 # m/s

modIP = 1. - np.concatenate((np.zeros(300), intrinsic_profile, np.zeros(300)))
modDV = - np.arange(np.floor(-modIP.size/2.+.5), np.floor(modIP.size/2.+.5)) * dbeta * an.c / 1e3
flineSpline2 = interpolate.UnivariateSpline(modDV[::-1], modIP[::-1], k=1., s=0.) 


observation = 1. - kerns.mean(1).ravel()
observation_norm = observation.copy()
for ii in range(nobs):
    i0, i1 = ii*nk, (ii+1)*nk
    inds = np.concatenate((np.arange(i0, i0+7), np.arange(i1-7, i1)))
    continuumfit = np.polyfit(inds, observation_norm[inds], 1)
    observation_norm[i0:i1] /= np.polyval(continuumfit, np.arange(i0, i1))



# Start the Imaging Analysis
nlat, nlon = 30, 60
alpha = 1000

nx = nlat*nlon
dime.setup(observation_norm.size, nk)
flatguess = 10*np.ones(nlat*nlon)
testincs = np.arccos(np.linspace(1, 0.4, 9))

allfits = []
for inc in testincs:
  map = maps.map(nlat=nlat, nlon=nlon, i=inc)
  ncell = map.ncell
  latlon_corners = np.array([c.corners_latlon for c in map.cells])
  phi_corner = latlon_corners[:,1,:]
  theta_corner = latlon_corners[:,0,:]

  LLD = 1.0
  # Compute the "Response Matrix" of Vogt 1987:
  Rmatrix = np.zeros((ncell, nobs*dv.size), dtype=float32)
  for kk, rot in enumerate(phase):
      speccube = np.zeros((ncell, dv.size), dtype=float32)
      this_map = maps.map(nlat=nlat, nlon=nlon, i=inc, deltaphi=-rot)
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
    perfect_fit = optimize.fmin_l_bfgs_b(dime.entropy_map_norm_sp, flatguess, fprime=dime.getgrad_norm_sp, args=fitargs, maxfun=10000, bounds=[(1e-6, 20)]*nx, factr=1e7, iprint=5)
    perfect_model = dime.normalize_model(np.dot(perfect_fit[0], Rmatrix), nk)
    w_observation /=  w_observation.max() * (sc_observation_norm - perfect_model)[w_observation>0].std()**2

  fitargs = (sc_observation_norm, w_observation, Rmatrix, alpha)
  bfit = optimize.fmin_l_bfgs_b(dime.entropy_map_norm_sp, flatguess, fprime=dime.getgrad_norm_sp, args=fitargs, maxfun=10000, bounds=[(1e-6, 20)]*nx, factr=1e7, iprint=5)
  allfits.append(bfit)
  bestparams = bfit[0]
  model_observation = dime.normalize_model(np.dot(bestparams, Rmatrix), nk)
  metric, chisq, entropy = dime.entropy_map_norm_sp(bestparams, *fitargs, retvals=True)

  hiphi = np.linspace(0, 2*np.pi, nlat*5)
  hitheta = np.linspace(0, np.pi, nlon*5)
  himap = griddata(np.concatenate((phi_corner-2*np.pi, phi_corner, phi_corner+2*np.pi)).ravel(), np.concatenate((theta_corner, theta_corner, theta_corner)).ravel(), np.concatenate((np.tile(bestparams, (4,1)).T, np.tile(bestparams, (4,1)).T, np.tile(bestparams, (4,1)).T)).ravel(), hiphi, hitheta)

  fig=figure(tools.nextfig(), [18, 7])
  fig.text(0.5, 0.9, 'MaxEnt Results: $i=%i deg, a=%i, \chi^2=%1.1f, S=%1.3f, Q=%1.2f $' % (inc*180/np.pi, alpha, chisq, entropy, metric), horizontalalignment='center')
  for ii in range(7):
    subplot(2, 7, ii+1)
    jj = ii*2
    hiphi2, hitheta2 = np.meshgrid(hiphi, hitheta)
    pmap = Basemap(projection='ortho', lat_0=-inc*180./np.pi, lon_0=phase[jj]*180/np.pi)
    x1, y1 = pmap(360. - hiphi2*180./np.pi, -hitheta2*180./np.pi +  90)
    CS = pmap.contourf(x1, y1, himap, np.linspace(bestparams.min(),bestparams.max(), 21)) #, cmap=cm.YlOrBr)
    pmap.drawmeridians(np.arange(0, 360, 30)); 
    pmap.drawparallels(np.arange(-90, 90, 30))
    title('%1.1f hr' % (phase[jj]*per/(3600.*2*np.pi)))
    subplot(2,7,ii+8)
    plot(dv, sc_observation_norm[jj*nk:(jj+1)*nk], 'k')
    plot(dv, model_observation[jj*nk:(jj+1)*nk], 'r')
    xlabel('RV [km/s]', fontsize=16)
    if jj==0: ylabel('Intensity', fontsize=16)
    ylim([.95, 1.003]); minorticks_on()


tools.printfigs('firstresults_nlat=%i_nlon=%i_alpha=%i_vsini=%1.1f.pdf' % (nlat, nlon, alpha, vsini/1e3))


##################################################
2013-08-11 13:33 IJMC: 
##################################################

Made considerable progress in the mapping routines; the script
'dopplermap.py' now does a credible job.  However, maps produced for
both A & B components exhibit zonal banding (like Jupiter), which
suggests a mismatch between the assumed and true instrinsic line
profile.

Let's try simulating some spectroscopic observations, run LSD on it,
and then invert the profiles.


##################################################
# Generate a simulated, rotationally-broadened spectrum:
##################################################

jj=1
modelfn = '/Users/ianc/proj/pcsa/data/model/lte015-5.0-0.0a+0.0.BT-Settl.spec.7.fits'
specmodelfn = '/Users/ianc/proj/pcsa/data/model/lte014-4.5-0.0a+0.0.BT-Settl.spec.7.fits'
specmodelfn = '/Users/ianc/proj/pcsa/data/model/lte013.5-5.5-0.0a+0.0.BT-Settl.spec.7.fits'
specmodelfn = '/Users/ianc/proj/pcsa/data/model/lte014-5.0-0.0a+0.0.BT-Settl.spec.7.fits'

model = pyfits.getdata(specmodelfn)
vsini = 30e3
inc = 0.2
nlat, nlon = 20, 40
lolim = chiplams[:,jj].min() - 0.003
hilim = chiplams[:,jj].max() + 0.003
lolim2, hilim2 = lolim+0.001, hilim -0.001
tind = (model[0]>lolim) * (model[0] < hilim)
lam_template = model[0,tind]
template = model[1,tind]
template /= np.median(template)
specSpline = interpolate.UnivariateSpline(lam_template, template, k=1., s=0.)

LLD=0.6
specRmatrix = np.zeros((nlat*nlon, nobs*npix), dtype=float32)
for kk, rot in enumerate(phase):
    this_map = maps.map(nlat=nlat, nlon=nlon, i=inc, deltaphi=-rot)
    speccube = np.zeros((this_map.ncell, npix), dtype=float32)
    this_doppler = 1. + vsini*this_map.visible_rvcorners.mean(1)/an.c/np.cos(inc)
    good = (this_map.projected_area>0) * np.isfinite(this_doppler)
    for ii in good.nonzero()[0]:
        speccube[ii,:] = specSpline(chiplams[kk,jj] * this_doppler[ii])
    limbdarkening = (1. - LLD) + LLD * this_map.mu
    specRblock = speccube * ((limbdarkening*this_map.projected_area).reshape(this_map.ncell, 1)*np.pi/this_map.projected_area.sum())
    specRmatrix[:,npix*kk:npix*(kk+1)] = specRblock

# Generate some simulated data:
fakemap = ones((nlat, nlon))
xx, yy = np.meshgrid(np.arange(nlon), np.arange(nlat))
for ll in range(3):
  fakemap[np.sqrt((xx+6-12*(ll+1))**2 + (yy-10)**2) <=3.5] = 0.9 


fakemap = fakemap.ravel()
simdata = np.dot(fakemap, specRmatrix).reshape(nobs, npix)
ss = np.argsort(simdata, axis=1)
simdata = (simdata/array([simdata[ll,ss[ll,-64:]] for ll in range(nobs)]).mean(1).reshape(nobs,1))
simdata += np.random.normal(size=simdata.shape)/100.

####################
# Run LSD on the fake data:
####################

# Prepare for Least-Squares Deconvolution:
pmod = os.path.split(modelfn)[1].replace('BT-Settl.spec.7.fits', '_dao').replace('.',',')
f_linelist   = _mod + pmod + '_edited.clineslsd'
f_pspec      = _mod + pmod + '.fits'
f_pspec_cont = _mod + pmod + 'C.fits'
(lineloc, lineew, linespec) = dao.getlines(f_linelist)  
pspec_cont = pyfits.getdata(f_pspec_cont) 
hdr_pspec_cont = pyfits.getheader(f_pspec_cont) 
wspec = hdr_pspec_cont['crval1'] + np.arange(pspec_cont.size)*hdr_pspec_cont['cdelt1']
spline = interpolate.UnivariateSpline(wspec, pspec_cont, s=0.0, k=1.0)
deltaspec0 = ns.linespec(lineloc, lineew, wspec, verbose=False, cont=pspec_cont)
tind2 = (wspec/1e4>lolim) * (wspec/1e4 < hilim)
deltaspec0 = deltaspec0[tind2]
wspec = wspec[tind2]


# Run LSD:
nk = 123
kerns = np.zeros((nobs, 4, nk), dtype=float)
modkerns = np.zeros((nobs, 4, nk), dtype=float)
for kk in range(nobs):
  deltaspec = ns.linespec(lineloc*(1.+9e-5), lineew, chiplams[kk,jj]*1e4, verbose=False, cont=spline(chiplams[kk,jj]*1e4))
  m,kerns[kk,jj],b,c = dia.dsa(deltaspec, simdata[kk], nk)    
  m,modkerns[kk,jj],b,c = dia.dsa(deltaspec, chipmodnobroad[kk,jj]/chipcors[kk,jj], nk)    

dbeta = np.diff(chiplams).mean()/chiplams.mean()
dx = -dbeta * arange(np.floor(-nk/2.+.5), np.floor(nk/2.+.5))
dv = an.c*dx / 1e3 # km/s

intrinsic_profile = modkerns[:,chips].mean(0).mean(0)
intrinsic_profile -= intrinsic_profile[0] + (intrinsic_profile[-1] - intrinsic_profile[0])/(nk - 1.) * np.arange(nk)
systematic_rv_offset = (intrinsic_profile==intrinsic_profile.max()).nonzero()[0] - (dv==0).nonzero()[0]
intrinsic_profile = np.interp(np.arange(nk), np.arange(nk) - systematic_rv_offset, intrinsic_profile)
sro2 = 17.25
sro2 = (kerns[:,jj].mean(0)==kerns[:,jj].mean(0).max()).nonzero()[0] - (dv==0).nonzero()[0]
kerns = np.array([[np.interp(np.arange(nk), np.arange(nk) - sro2, kerns[ll, kk]) for kk in range(4)] for ll in range(nobs)])

# Prepare for Doppler Imaging:
fnsLi = ([['specL_chip%i_frame%i.fits' % (ii,kk) for ii in range(1,5)] for kk in range(14)])
mjd = np.array([pyfits.getval(fn[0], 'MJD-OBS') for fn in fnsLi])
phase = (mjd-mjd[0])*86400 * 2*np.pi/ per

modIP = 1. - np.concatenate((np.zeros(300), intrinsic_profile, np.zeros(300)))
modDV = - np.arange(np.floor(-modIP.size/2.+.5), np.floor(modIP.size/2.+.5)) * dbeta * an.c / 1e3
flineSpline2 = interpolate.UnivariateSpline(modDV[::-1], modIP[::-1], k=1., s=0.) 

observation = 1. - kerns[:,jj].ravel()
observation_norm = observation.copy()
for ii in range(nobs):
    i0, i1 = ii*nk, (ii+1)*nk
    inds = np.concatenate((np.arange(i0, i0+7), np.arange(i1-7, i1)))
    continuumfit = np.polyfit(inds, observation_norm[inds], 1)
    observation_norm[i0:i1] /= np.polyval(continuumfit, np.arange(i0, i1))

# Now invert the line profiles:
dime.setup(observation_norm.size, nk)
flatguess = 100*np.ones(nx)
bounds = [(1e-6, 300)]*nx

map = maps.map(nlat=nlat, nlon=nlon, i=inc)
ncell = map.ncell
latlon_corners = np.array([c.corners_latlon for c in map.cells])
phi_corner = latlon_corners[:,1,:]
theta_corner = latlon_corners[:,0,:]

# Compute the "Response Matrix" of Vogt 1987:
Rmatrix = np.zeros((ncell, nobs*dv.size), dtype=float32)
for kk, rot in enumerate(phase):
    speccube = np.zeros((ncell, dv.size), dtype=float32)
    this_map = maps.map(nlat=nlat, nlon=nlon, i=inc, deltaphi=-rot)
    this_doppler = 1. + vsini*this_map.visible_rvcorners.mean(1)/an.c/np.cos(inc)
    good = (this_map.projected_area>0) * np.isfinite(this_doppler)
    for ii in good.nonzero()[0]:
        speccube[ii,:] = flineSpline2(dv + (this_doppler[ii]-1)*an.c/1000.)
        #speccube[ii,:] = flineSpline(lamdv / this_doppler[ii])
    limbdarkening = (1. - LLD) + LLD * this_map.mu
    Rblock = speccube * ((limbdarkening*this_map.projected_area).reshape(this_map.ncell, 1)*np.pi/this_map.projected_area.sum())
    Rmatrix[:,dv.size*kk:dv.size*(kk+1)] = Rblock

flatmodel = dime.normalize_model(np.dot(flatguess, Rmatrix), nk)
minfm = flatmodel.min()
cutoffval = 1. - (1. - minfm) / 22.
w_observation = (flatmodel < cutoffval).astype(float) / np.median(kerns[:,jj].std(0))**2
# Scale the observations to match the model's equivalent width:
out, eout = an.lsq((observation_norm, np.ones(nobs*nk)), flatmodel, w=w_observation)
sc_observation_norm = observation_norm * out[0] + out[1]

fitargs = (sc_observation_norm, w_observation, Rmatrix, 20000)
perfect_fit = an.gfit(dime.entropy_map_norm_sp, flatguess, fprime=dime.getgrad_norm_sp, args=fitargs, ftol=ftol, disp=1, maxiter=1e4, bounds=bounds)
pp = dime.normalize_model(np.dot(perfect_fit[0], Rmatrix), nk)



##################################################
2013-08-11 16:08 IJMC: 

Tests show that using near-saturated lines does not cause LSD and
subsequent Doppler Imaging to fail: gross features can be recovered in
simulated data (using BT-Settl model spectra) using the same
techniques as applied to real data.  These simulations also show mean
zonal banding of only +/- 1%: near-zero at the poles, 1% darkening at
equator, and ~1% brightening at midlatidues.

These simuations also show that for ~1000 pixels, alpha-values of
1000-10000 should be employed. If single-pixel features are being
resolved: don't trust them!

Finally, it's not clear that even such model mismatches cause the
banding structure I observe in the objects' recovered maps.

The largest-scale features are recovered even for atmospheric model
parameters ~200K or 0.5 dex (in logg) from the best-matching values.


##################################################
2013-08-12 09:34 IJMC: 

Now try simulating data using Derek Homeier's fancy model:

jj=1
modelfn = '/Users/ianc/proj/pcsa/data/model/lte015-5.0-0.0a+0.0.BT-Settl.spec.7.fits'

_dat = os.path.expanduser('~') + '/proj/pcsa/data/raw/20130505/'
dat = tools.loadpickle(_dat + 'fainterspectral-fits_6.pickle')
nobs = 14
for key in dat.keys():
  exec('%s = dat["%s"]' % (key, key))


allfits = array([[ii[0] for ii in kk] for kk in individual_fits])



specmodelfn = '/Users/ianc/proj/pcsa/data/model/lte1400-5.00-0.0a+0.0.BT-settl-giant-2012.cifist.tm1.0-1.1.R0-g12p0co.irf.h5'
specmodelfn = '/Users/ianc/proj/pcsa/data/model/lte1400-5.00-0.0a+0.0.BT-settl-giant-2012.cifist.tm1.0-1.1.R0-g12p0k7.irf.h5'
fp = tables.open_file(specmodelfn, 'r')
modelMu = fp.root.Atmosphere.murf.read()
modelAtmo = fp.root.Atmosphere.spectrum.read_where('(22800<=wl)&(wl<=23500)')
fp.close()
lam_template = modelAtmo['wl']/1e4


def getMuSpec(mu):
  muInds, fracs = tools.findFrac(modelMu, mu, retinds=True)
  spec = np.sum([frac*modelAtmo['rfmu'][:,ind] for frac, ind in zip(fracs, muInds)], axis=0)
  return spec

vsini = 30e3
inc = 0.2
nlat, nlon = 20, 40

#lolim = chiplams[:,jj].min() - 0.003
#hilim = chiplams[:,jj].max() + 0.003
#lolim2, hilim2 = lolim+0.001, hilim -0.001
#tind = (model[0]>lolim) * (model[0] < hilim)
#lam_template = model[0,tind]
#template = model[1,tind]
#template /= np.median(template)
#specSpline = interpolate.UnivariateSpline(lam_template, template, k=1., s=0.)

specRmatrix = np.zeros((nlat*nlon, nobs*npix), dtype=float32)
for kk, rot in enumerate(phase):
    this_map = maps.map(nlat=nlat, nlon=nlon, i=inc, deltaphi=-rot)
    speccube = np.zeros((this_map.ncell, npix), dtype=float32)
    this_doppler = 1. + vsini*this_map.visible_rvcorners.mean(1)/an.c/np.cos(inc) #- allfits[:,:,2].mean()
    good = (this_map.projected_area>0) * np.isfinite(this_doppler)
    for ii in good.nonzero()[0]:
        speccube[ii,:] = getMuSpec(this_map.mu[ii]) #np.interp(chiplams[kk,jj] * this_doppler[ii], lam_template, getMuSpec(this_map.mu[ii]))
    limbdarkening = 1. #########(1. - LLD) + LLD * this_map.mu
    specRblock = speccube * ((limbdarkening*this_map.projected_area).reshape(this_map.ncell, 1)*np.pi/this_map.projected_area.sum())
    specRmatrix[:,npix*kk:npix*(kk+1)] = specRblock

# Generate some simulated data:
fakemap = ones((nlat, nlon))
#xx, yy = np.meshgrid(np.arange(nlon), np.arange(nlat))
#for ll in range(3):
#  fakemap[np.sqrt((xx+6-12*(ll+1))**2 + (yy-10)**2) <=3.5] = 0.5 


fakemap = fakemap.ravel()
simdata = np.dot(fakemap, specRmatrix).reshape(nobs, npix)
ss = np.argsort(simdata, axis=1)
simdata = (simdata/array([simdata[ll,ss[ll,-64:]] for ll in range(nobs)]).mean(1).reshape(nobs,1))
simdata += np.random.normal(size=simdata.shape)/100.


##################################################
2013-08-19 16:33 IJMC: 
##################################################

Look at the results of MCMC fitting for the different spectral orders,
and for using data when the spot is (a) visible and (b) not visible.

Here's what you do:

# --------------------------------------------------
# Run the beginning of mcmc_spectral_fits. Then:

fns = ['fainterspectral-fits_6_mcmcfits_guess--teff=1340_logg=4.20.pickle', 
       'fainterspectral-fits_6_mcmcfits_guess--teff=1450_logg=4.50.pickle', 
       'fainterspectral-fits_6_mcmcfits_guess--teff=1340_logg=5.40.pickle', 
       'fainterspectral-fits_6_mcmcfits_guess--teff=1540_logg=5.40.pickle',
       'fainterspectral-fits_6_mcmcfits_guess--teff=1320_logg=4.20_inds=45678.pickle',
       'fainterspectral-fits_6_mcmcfits_guess--teff=1450_logg=4.50_inds=45678.pickle',
       'fainterspectral-fits_6_mcmcfits_guess--teff=1340_logg=5.40_inds=45678.pickle',
       'fainterspectral-fits_6_mcmcfits_guess--teff=1540_logg=5.40_inds=45678.pickle',
       'fainterspectral-fits_6_mcmcfits_guess--teff=1340_logg=4.20_inds=01210111213.pickle',
       'fainterspectral-fits_6_mcmcfits_guess--teff=1450_logg=4.50_inds=01210111213.pickle',
       'fainterspectral-fits_6_mcmcfits_guess--teff=1340_logg=5.40_inds=01210111213.pickle',
       'fainterspectral-fits_6_mcmcfits_guess--teff=1540_logg=5.40_inds=01210111213.pickle']


nfns = len(fns)
dat = map(tools.loadpickle, fns)

rmses = np.zeros((nfns, 4), dtype=float)
teffs = np.zeros((nfns, 4), dtype=float)
loggs = np.zeros((nfns, 4), dtype=float)
models = np.zeros((nfns, 4, 1024), dtype=float)
for ii in range(nfns):
  for jj in range(4):
    d = dat[ii]
    args = (interpModel,) + d['args'][jj]
    obs = args[-3]
    wobs = args[-2]
    mask = wobs > 1
    models[ii,jj], bestw = interpModel(d['bestfits'][jj], *args[1:-3], retlam=True)
    rmses[ii, jj] =     (obs - models[ii,jj])[mask].std()
    teffs[ii,jj] = d['bestfits'][jj][0]
    loggs[ii,jj] = d['bestfits'][jj][1]

# --------------------------------------------------


##################################################
2013-08-25 17:06 IJMC: 
##################################################

Look at the 'rmses' and see which initial conditions lead to the
best-fitting models.  This seems to be the initialization with high
surface gravity.

Then, look at how the best-fit Teff and Logg change for the two sets!

# Run the beginning of mcmc_spectral_fits. Then:

import os
import tools

lis = os.listdir('.')
fns = []
for fn in lis:
  if fn.find('=45678.')>0: fns.append(fn)

  if fn.find('=01210111213.')>0: fns.append(fn)

  if fn.find('=012345')>0: fns.append(fn)


nfns = len(fns)


rmses = np.zeros((nfns, 4), dtype=float)
teffs = np.zeros((nfns, 4), dtype=float)
loggs = np.zeros((nfns, 4), dtype=float)
models = np.zeros((nfns, 4, 1024), dtype=float)
for ii in range(nfns):
  d = tools.loadpickle(fns[ii])
  for jj in range(4):
    args = (interpModel,) + d['args'][jj]
    obs = args[-3]
    wobs = args[-2]
    mask = wobs > 1
    models[ii,jj], bestw = interpModel(d['bestfits'][jj], *args[1:-3], retlam=True)
    rmses[ii, jj] =     (obs - models[ii,jj])[mask].std()
    teffs[ii,jj] = d['bestfits'][jj][0]
    loggs[ii,jj] = d['bestfits'][jj][1]


bestinds = [(rmses[:,ii]==rmses[:,ii].min()).nonzero()[0][0] for ii in range(4)]
besttemps  = np.zeros(4)
ebesttemps = np.zeros(4)
bestloggs  = np.zeros(4)
ebestloggs = np.zeros(4)
for ii in range(4):
  d = tools.loadpickle(fns[bestinds[ii]])
  besttemps[ii], bestloggs[ii] = d['bestfits'][ii][0:2]
  ebesttemps[ii] = an.dumbconf(d['chains'][ii][:,0], .683, mid=besttemps[ii])[0]
  ebestloggs[ii] = an.dumbconf(d['chains'][ii][:,1], .683, mid=bestloggs[ii])[0]


########################################
Preliminarily, this gives results for (clear minus spot) of:

T_clear - T_spot = -3.6 +/- 2.7 K
g_clear - t_spot = +0.020 +/- 0.013 (logg, cgs)

A subsequent analysis gives:
T_clear - T_spot = -4.0 +/- 1.6 K
g_clear - t_spot = +0.014 +/- 0.012 (logg, cgs)


------------------------------
2013-09-05 15:39 IJMC: 
------------------------------

TO examine smudge-plot possibilities:

cls = [-0.0025, .0025], [-0.0035, .0035], [-0.002, .003], [-0.0024, .003]
cminds = array([23, 27, 43, 77, 81, 87, 104, 106, 115, 119, 122, 128, 129, 138])-1
cminds = array([ 26,  42,  76,  86, 118, 127, 128, 137])
cminds = [86]
cms = cm.datad.keys()
ioff()
for cmind in cminds:
  fig = figure()
  fig.text(0.5, 0.95, cms[cmind]) 
  for jj,cl0 in enumerate(cls):
    subplot(2,2,jj+1)
    imshow(kernstack, x=dv, y=hr, cmap=cms[cmind])
    ylabel('Elapsed Time [hr]', fontsize=fs*0.8)
    xlabel('Velocity [km/s]', fontsize=fs*0.8)
    minorticks_on()
    xlim(xl)
    yl = ylim()
    plot([26.1]*2, ylim(), '--k', linewidth=3)
    plot([-26.1]*2, ylim(), '--k', linewidth=3)
    ylim(yl)
    clim(cl0)
    colorbar()

tools.printfigs('colormap_choices3.pdf')
close('all')


cl0 = [-0.0023999999999999998, 0.0030000000000000001]
map = 'gist_earth'
figure(3, [6,4])
ax=subplot(111, position=[.1, .15, .8, .8])
imshow(kernstack, x=dv, y=hr, cmap=map)
ylabel('Elapsed Time [hr]', fontsize=fs)
xlabel('Velocity [km/s]', fontsize=fs)
minorticks_on()
xlim(xl)
yl = ylim()
plot([26.1]*2, ylim(), '--', linewidth=3, color='k')
plot([-26.1]*2, ylim(), '--', linewidth=3, color='k')
ylim(yl)
clim(cl0)
colorbar()



--------------------------------------------------
2013-11-12 18:43 IJMC: 
--------------------------------------------------
# Notes to make a Gnomic projection-on-a-cube of Luhman 16B brown dwarf:

hiphi = np.linspace(0, 2*np.pi, nlat*50)
hitheta = np.linspace(0, np.pi, nlon*50)
hiphi2, hitheta2 = np.meshgrid(hiphi, hitheta)
himap = griddata(np.concatenate((phi_corner-2*np.pi, phi_corner, phi_corner+2*np.pi)).ravel(), np.concatenate((theta_corner, theta_corner, theta_corner)).ravel(), np.concatenate((np.tile(bestparams, (4,1)).T, np.tile(bestparams, (4,1)).T, np.tile(bestparams, (4,1)).T)).ravel(), hiphi, hitheta)


lat0s = [90, 0, 0, 0, 0, -90]
lon0s = [5, 5, 95, 185, 275, 5]
corners = [[.75, .3333, .25, .3333], 
	   [.25, .6666, .25, .3333], 
	   [.50, .3333, .25, .3333], 
	   [.25, .0,    .25, .3333], 
	   [.0,  .3333, .25, .3333], 
	   [.25, .3333, .25, .3333]]

figure(tools.nextfig(), figsize=[10, 7.5])
for hh in range(6):
  ax = subplot(111, position=corners[hh])
  fax = ax.get_frame()
  fax.set_linewidth(8)
  pmap = Basemap(width=12.8e6,height=12.8e6,projection='gnom',lat_0=lat0s[hh],lon_0=lon0s[hh])
  x1, y1 = pmap(360. - hiphi2*180./np.pi, -hitheta2*180./np.pi +  90)
  CS = pmap.contourf(x1, y1, himap, np.linspace(bestparams.min(),bestparams.max(), 21), cmap=cm.gist_heat) #np.linspace(0, 300, 21)) # #, cmap=cm.YlOrBr)   
  pmap.drawmeridians(np.arange(0, 360, 30)); 
  pmap.drawparallels(np.arange(-90, 90, 30))

savefig('cubeplots_vec_1fig_l8.png', transparent=True, dpi=350)
savefig('cubeplots_vec_1fig_l8.pdf', transparent=True)



