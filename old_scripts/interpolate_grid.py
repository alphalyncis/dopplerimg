import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from scipy.interpolate import griddata
from scipy import interpolate
import matplotlib.patches as patches
from uncertainties import ufloat
from uncertainties.umath import *
from astropy.io import fits

t = Table.read('hybrid_solar_age.dat', format='ascii')

# make interpolated grid of Lbol as a function of mass and age
minmass = 0.003
maxmass = 0.02
#minmass = 0.0087
#maxmass = 0.0187
age = np.linspace(0.005,0.03,10000)
mass = np.linspace(minmass, maxmass, 10000)
AGE, MASS = np.meshgrid(age,mass)
lbolgrid = griddata((t['age(Gyr)'], t['M/Msun']), t['logL/Lsun'], (AGE, MASS), method='linear')
teffgrid = griddata((t['age(Gyr)'], t['M/Msun']), t['Teff(K)'], (AGE, MASS), method='linear')
logggrid = griddata((t['age(Gyr)'], t['M/Msun']), t['logg'], (AGE, MASS), method='linear')
radiusgrid = griddata((t['age(Gyr)'], t['M/Msun']), t['R/Rsun'], (AGE, MASS), method='linear')
