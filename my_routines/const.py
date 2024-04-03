from astropy.constants import *
import numpy as np
import astropy.units as u

# jupiter's semimajor axis
a_jup = 1e6*778.479*u.km # https://nssdc.gsfc.nasa.gov/planetary/factsheet/jupiterfact.html

a_io = 1e3*422*u.km # https://nssdc.gsfc.nasa.gov/planetary/factsheet/galileanfact_table.html

# Io's orbital inclination w.r.t. Jupiter's axis of rotation
io_i = 0.04*u.deg # https://nssdc.gsfc.nasa.gov/planetary/factsheet/galileanfact_table.html

# Jupiter's obliquity w.r.t. its orbital inclination
beta_jup = 3.13*u.deg # https://nssdc.gsfc.nasa.gov/planetary/factsheet/jupiterfact.html

# jupiter's orbital inclination w.r.t. Earth's
jup_i = 1.304*u.deg # https://nssdc.gsfc.nasa.gov/planetary/factsheet/jupiterfact.html

h = 6.62607015e-34
hbar = h/2/np.pi # SI
c = 2.99792458e8 # m/s
fgamma = 7.0e-4
mamu = 1.66053906660e-27
mubar = 0.6 # for the Sun
G = 6.67430e-11 # SI
Msun = 1.9884e30 # kg
Rsun = 6.95e8 # m
h = 0.7
H0 = 100*h
mpc2m = 3.0856775814914e16*1e6

# in units of Hz
H0sec = H0*10.**(3)/(mpc2m)

eps0 = 3*c**2/8./np.pi/G*(H0*10.**(3)/(mpc2m))**2 # critical density
n0 = 5 # particle density in m^-3
