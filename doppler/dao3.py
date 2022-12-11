"""
# Original author: Ian Crossfield (Python 2.7)

#2008-12-01 17:52 IJC:  Created

"""
######################################################

# 18-02-2020
# Emma Bubb - change to run on Python 3

######################################################

import os, shutil, pdb

def runspec(w, model, out, wlim, fwhm, lineids=None, clobber=False, plotalot=False, daoexec='/Users/ianc/python/daospec_src/daospec', order=15):
    """Prepare and run DAOSpec on the planet spectra provided by Travis
       Barman, and previously prepared with lineid.pro and
       :func:`getspec`.

       INPUT: 
          w -- 1D array:  wavelength grid for input spectrum (Angstroms)
          model -- 1D array:  input spectrum 
          out -- str: prefix name of the output files
          wlim -- float list; the desired wavelength limits for DAOSpec, in
                         Angstroms.
          fwhm -- float; the FWHM of features, in pixels

       OPTIONAL INPUTS:
          lineids  -- str: *.7.gz filename containing Barman model
          daoexec -- str: path to 'daospec' executable.
          order -- int: polynomial order for continuum fit.

       EXAMPLE: 

         try:
             from astropy.io import fits as pyfits
         except:
             import pyfits
         , dao
         _path = '/Users/ianc/proj/pcsa/data/model/'
         Lfn = _path + 'lte0125-3.5.rainout.HD75732e.redist=0.50.hires_id.7.gz'
         Mfn = _path + 'lte0125-3.5.rainout.HD75732e.redist_0.50.hires.7.fits'
         M = pyfits.getdata(Mfn)
         out = '55cnce_test_output'
         dao.runspec(M[0], M[1], out, [2e4, 2.4e4], 9.5, lineids=Lfn, clobber=True)

       :NOTE:
          You need DAOSPEC to run this routine; it is available here:
          http://www.bo.astro.it/~pancino/projects/daospec.html
          """
    #2010-11-09 12:06 IJC:  Added note as to location of DAOSPEC.  Added 'daoexec' option.
    #2010-11-10 IJC: Rather changed script structure; syntax was corrupted/outdated.
    #2014-01-23 15:47 IJMC: Increased flexibility of searching for
    #                       'daoout' file. Not sure why this was
    #                       suddenly necessary!
    
    ### Dependencies:  DAOspec, lineid.pro, rdidgf.pro, getmodellines.pro


    try:
        from astropy.io import fits as pyfits
    except:
        import pyfits
    
    from numpy import array, diff, float32

#lineids  = 'lte0125-3.5.rainout.HD75732e.redist=0.50.hires_id.7.gz'
#model    = 'lte0125-3.5.rainout.HD75732e.redist_0.50.hires.7.fits'
#model = '55cnce_barman_model.fits'
#linelist = '55cnce_lines.txt'  # Model linelist, extracted by lineid.pro

#wlim = [2.03e4, 2.4e4]
#fwhm = 9.5  # FWHM, in pixels

    linelistfn = out+'.txt'
    modelfn = out+'.fits'
    linesout = out+'.clines'

    # Rewrite the input spectrum in DAOSpec-friendly 1D format, with
    # correct header keywords:
    c1 = pyfits.Card('CRPIX1', 1, 'Reference pixel')
    c2 = pyfits.Card('CRVAL1', w[0], 'Coordinates at reference pixel')
    c3 = pyfits.Card('CDELT1', diff(w).mean(), 'Coordinate increment per pixel')
    c4 = pyfits.Card('CTYPE1', 'LINEAR', 'Units of coordinate')
    print( w[0], diff(w).mean() ) # EB updated print statement
    h = pyfits.Header([c1,c2,c3,c4])
    hdu = pyfits.PrimaryHDU(model.astype(float32), h)
    hdu.writeto(modelfn, clobber=clobber)


 
# Use lineid.pro to get the line identifications
    if lineids != None: #EB updated not equal to sign <>
        if os.path.isfile('getmodellines.pro'):
            os.remove('getmodellines.pro')
        f = open('getmodellines.pro', 'w')
        f.write('.r lineid\n')
        f.write('init_line_id\n')
        f.write("rdidgf,'" + lineids + "',w,f,b\n")
        f.write('mol_strg=15\n')
        f.write('!x.range=[' + str(wlim[0]) + ',' + str(wlim[1]) + ']\n')
        f.write('!y.range=[0,0]\n')
        f.write('!x.style=1\n')
        f.write('plot,w,f\n')
        f.write("pmolid,'" + linelistfn + "', 0.65, 0.7, 1, mol_strg\n")
        f.write('exit\n')
        f.close()

        idlexec = os.popen('which idl').read().strip()
        os.system(idlexec + ' getmodellines')

        shutil.copy(linelistfn, 'laboratory.dat')
    else:
        if os.path.isfile('laboratory.dat'):
            os.remove('laboratory.dat')

    logname = modelfn + '.log'

#Set up and run DAOSpec
    f = open('daospec.sh', 'w')
    f.write('# !' + os.popen('which csh').read().strip() + '\n')
    f.write('%s << DONE > %s \n' % (daoexec, logname)  )
    f.write('or=%i \nfw=%1.2f\nfi=0 \n' % (order, fwhm) )
    f.write('sh='+str(wlim[0])+' \nlo=' + str(wlim[1])+'\n')
    f.write('re=0 \nmi=-0.5 \nma=0.5 \nve=3 \ncr=1\n')
    f.write('wa=0 \nsm=0.5 \nre=20 \nsc=1\n')
    f.write('\n' + modelfn + '\n \nDONE \n')
    f.close()

    #f = open('daospec.sh', 'w')
    #f.write('# !' + os.popen('which csh').read().strip() + '\n')
    ##if os.path.isfile(daoout):
    ##    f.write('rm ' + '\n')
    #f.write('%s << DONE \n' % daoexec  )
    #f.write('\n' + modelfn + '\n \nDONE \n')
    #f.close()
    
    os.system('chmod u+x daospec.sh')
    os.system('./daospec.sh')

    daoout_found = False
    p = os.popen('ls %s.K*' % out)
    filenames = [out+'.daospec'] + p.readlines()
    p.close()
    file_iter = 0
    while not daoout_found:
        this_filename = filenames[file_iter].strip()
        if os.path.isfile(this_filename):
            daoout = this_filename
            daoout_found = True
        file_iter += 1

#    if os.path.isfile(    daoout.replace('.daospec', 'lte013,5-4,5-0,0a+0,0,_dao.KEYBOARD INPU



# Confirm and label lines
    f = open(daoout, 'r')
    f.readline(); f.readline()
    raw = f.readlines()
    f.close()
    conf = []
    ew   = []
    eew  = []
    spec = []
    for row in raw:
        rowel = row.split()
        conf.append((rowel[0]))
        ew.append((rowel[2]))
        eew.append((rowel[3]))
        if len(rowel)>5:
            spec.append(rowel[-1])
        else:
            spec.append('unidentified')
    f = open(linesout, 'w')
    for ii in range(len(conf)):
        f.write(str(conf[ii]) + '\t' + str(ew[ii]) + '\t' + 
                str(eew[ii]) + '\t' + spec[ii] + '\n')
    f.close()

# Re-measure the strengths of the located lines:
    pspec = pyfits.getdata(modelfn)
    pspec_cont = pyfits.getdata(modelfn.replace('.fits', 'C.fits'))
    lines = array(conf).astype(float)
    
    if len(lines)>0:
        try:
            (lineloc, linestrength)=findlinedepths(w, lines, pspec, pspec_cont, 4.4*diff(w).mean(), dstep=2346., plotalot=plotalot)
        except:
            (lineloc, linestrength)=findlinedepths(w, lines, pspec, pspec_cont, 4.4*diff(w).mean(), dstep=1794., plotalot=plotalot)
        f = open(linesout+'lsd', 'w')
        for ii in range(len(lineloc)):
            f.write(str(lineloc[ii]) + '\t' + str(linestrength[ii]*1000.) + '\t' + 
                    '-1.0\tIJC_gaussian_linefit\n')
        f.close()
    else:
        print( "Did not find any lines!!  Exiting.") # EB updated print statement

    return


def getlines(f_linelist):
    """Read the line locations and equivalent widths from a DAOSPEC output file.

    Example:
      f_linelist = 'model_spec.clines'
      (lineloc, lineew, linespec) = getlines(f_linelist)
      """
    #2009-02-22 10:15 IJC: Initiated

    from numpy import zeros

    # Get the line locations and EWs:
    f = open(f_linelist, 'r')
    raw = f.readlines()
    f.close()
    
    # dat = array([map(float, line.split()[0:2]) for line in raw]) #EB the map function does not work the same in python 3
    
    #EB - this loop gets the form needed
    dat = zeros([len(raw), 2], dtype=float)
    i =0                                                    
    for line in raw:                                         
        dat[i,:]= list(map(float, line.split()[0:2]))
        i=i+1

    lineloc = dat[:,0]
    lineew = dat[:,1]/1e3
    linespec = [line.split()[-1] for line in raw]

    return (lineloc, lineew, linespec)



def findlinedepths(w, line, spec, cont, sigma, dstep=2000., plotalot=False):
    """Fit line strengths, if you know the location.

    Returns the line locations and equivalent widths (area of Gaussian
    profiles) for the given continuum level.

    This attempts to account for line blending using matrix inversion;
    if W is longer than several thousand elements, the matrix will be
    prohibitively large, and discontinuities may result.  A preferable
    solution is to break up your spectrum into a number of smaller
    chunks.

    
    :INPUTS:
      w - wavelength grid

      line - list of line locations (in units of w)

      spec - spectrum for which line widths will be crudely measured.

      cont - assumed continuum level used, appropriate for 'spec'
      
      sigma - Assumed sigma ( = FWHM/2.35) of each (assumed Gaussian)
              line profile, in units of wavelength grid!

      dstep - Number of pixels in each segment, when breaking up a
              spectrum into smaller sections.

    :Returns:
      (locations, coefficients)

    :TO_DO:
      Allow different line profile shapes.
      """

    # 2009-03-05 19:37 IJC: Created @ UCLA
    # 2013-05-27 08:51 IJMC: Added useful documentation, in midair over Brazil.

    from fit_atmo import gaussian
    from scipy import interpolate
    from numpy import array, arange, zeros, float, dot, concatenate
    from pylab import pinv, figure, plot, legend

    line = array(line).copy()
    nlam = len(w)
    gparam = [1.0, sigma, 0.0]

    coeflist = array([])
    linelist = array([])
    newspec = array([])
    nstep = int(nlam/dstep)+1
    #print( len(line), nlam, nstep) # EB updated print statement
    for ii in range(nstep):
        windex = arange(ii*dstep,min((ii+1)*dstep,len(w)), dtype=int)
        ind1 = (line>w[windex].min())
        ind2 =  (line<w[windex].max())
#        if 
        lineind = (line>w[windex].min()) * (line<w[windex].max())

        theselines = line[lineind]
        nline = len(theselines)
        nlam = len(w[windex])
        kernelmatrix = zeros((nlam, nline), float)
        print( ii, w[windex].min(), w[windex].max(), lineind.sum(), kernelmatrix.shape) # EB updated print statement
        print( theselines) # EB updated print statement
        for jj in range(nline):   
            gparam[2]=float(theselines[jj]);  
            #print( ii,jj, w[windex].shape, gparam) # EB updated print statement
            kernelmatrix[:,jj]=gaussian(gparam, w[windex])

        svec = -(spec[windex]-cont[windex])
        if nline>0:
            invkern = pinv(kernelmatrix)
            coef = dot(invkern, svec)
            linelist = concatenate((linelist, theselines))
            coeflist = concatenate((coeflist, coef))
            newspec = concatenate((newspec, dot(kernelmatrix,coef)))
            if plotalot:
                figure()
                plot(w[windex], spec[windex], w[windex], cont[windex], \
                         w[windex], cont[windex]-dot(kernelmatrix,coef))
                legend(['spec', 'cont', 'cont-coef*kern'])

    ### Correct for continuum level to calculate equivalent width:
    cspline = interpolate.UnivariateSpline(w, cont, s=0.0, k=3.)
    coefcorrection = cspline(linelist)
    coeflist = coeflist/coefcorrection

    if plotalot:
        from dia import rconvolve1d, dsa
        from nsdata import linespec
        gparam2 = [1.0, 4.4, 0]
        deltaspec = linespec(linelist, coeflist, w, verbose=False, cont=cont)
        kern = gaussian(gparam2, arange(-100.,100.)+0.5)
        testspec = rconvolve1d(deltaspec, kern, 'same')
        figure()
        plot(w, spec); plot(w, cont); plot(w, deltaspec, '--'); plot(w, testspec)
        legend(['spec', 'cont', 'deltaspec', 'testspec'])
        tempindex = (cont>0)
        dsa(deltaspec[tempindex], spec[tempindex], 200, verbose=True)
        import pickle
        f = open('pickletemp', 'w')
        pickle.dump([linelist, coeflist, w[tempindex], spec[tempindex], \
                         cont[tempindex], kern], f)
        f.close()
    return (linelist, coeflist)
    

def getspec(input, output=None, clobber=False):
    """Simple script to put the necessary DAOSPEC FITS header keys in a
    spectrum.  Wavelength scale must be LINEAR!

    :INPUTS:
      input : str or 2xN NumPy array
        Spectral model. The spectrum should be of shape 2xN, where the
        first row is the wavelength and the second is the spectral
        data. If a string, the filename of a FITS file containing such data.

      output : str or None
        Output filename. If None, file will be overwritten (so set clobber=True)

      clobber : bool
        Whether or not an existing file should be overwritten.
        """
    #2008-12-01 17:22 IJC: Created
    # 2013-05-11 10:31 IJMC: COnverted from a crude script to a function.


    try:
        from astropy.io import fits as pyfits
    except:
        import pyfits
    
    import numpy as np

    if isinstance(input, str):
        dat = pyfits.getdata(input)
        hdr = pyfits.getheader(input)
        if output is None:
            output = input
    else:
        dat = input
        hdr = pyfits.header.Header()


    if dat.ndim==2:
        w = dat[0,:]
        spec = np.array(dat[1,:], dtype=np.float)

        crpix = 1
        crval = w[0]
        cdelt = np.diff(w).mean()

        hdr.update('CRPIX1', crpix, comment='Reference pixel')
        hdr.update('CRVAL1', crval, comment='Coordinate at reference pixel')
        hdr.update('CDELT1', cdelt, comment='Coordinate increment per pixel')
        hdr.update('CTYPE1', 'LINEAR', comment='Units of coordinate')

        pyfits.writeto(output, spec, header=hdr, clobber=clobber)
    else:
        print( "FITS array was not of size 2xN. No file written!") # EB updated print statement

    return

