#!/usr/bin/python
# the path in the line above should point to the location of the python
# executable; having this as the first line allows the script to be
# run directly provided you make it executable by chmod +x integrate.py
"""
integrate.py -- main program for running integration functions
  The integrand is defined in the function integrand() below.
  The integration functions themselves are in the file integrate_sub.py.
  Can be called from shell as integrate.py a b	  [see comment above]
    where a and b are the limits of integration 
  Can be called from ipython as '%run integrate.py a b'
"""

import math
import numpy as np		# may be used in integrand
import sys			# needed for command line reads below
from subprocess import call	# to allow a shell call to rename file

# import the integration functions, two different cases
from integrate_sub_starter import *

# define several control variables
nstepmax=2e8    	# maximum number of allowed integration steps
tolerance=1.e-6		# require convergence to this fractional error
verbose=1		# write iterations to output files

def integrand(x):
    '''
    x^{-0.5}
    '''
    return(x**(-1.5))

def int1(x):
    '''
    \frac{1}{x^{3/2}(1 + x^{3/2})}
    '''
    return (x**(1.5)*(1.+x**1.5))**(-1.)
limits1 = (1., 2.)

def int2(x):
    '''
    \frac{sin}{x}
    '''
    return (np.sin(x)/x)
limits2 = (1., 100.)

def int3(x):
    '''
    \frac{\sin^2}{x^2}
    '''
    return (np.sin(x)/x)**2
limits3 = (1., 1000.)

def int4(x):
    '''
    \frac{1}{x + 1/x}
    '''
    return (x + 1./x)**(-1.)
limits4 = (1., 1000.)

def int5(x):
    '''
    \frac{1}{1 + \exp{-2x}}
    '''
    return (1. + np.exp(-2.*x))**(-1.)
limits5 = (0., np.log(1000.))

def int6(x):
    '''
    x^{-0.5}
    '''
    return x**(-0.5)
limits6 = (0., 4.)

def int7(x):
    '''
    \frac{\sin}{x}
    i.e., the sine integral
    '''
    return np.sin(x)/x
limits7 = (0., np.pi)

def int8(x):
    '''
    \frac{1}{x^2 + x^3}
    '''
    return (1.+ x**(-1.))**(-1.)
limits8 = (0., 1.)

def int9(x):
    '''
    \frac{\sin^2x}{x^2}
    '''
    return np.sin(1./x)**2
limits9 = (0., 1.)

ints = [int1, int2, int3, int4, int5]
limits = [limits1, limits2, limits3, limits4, limits5]

################################################################
# 5)
################################################################
# ints = [int6, int7]
# limits = [limits6, limits7]

################################################################
# 6)
################################################################

# ints = [int8, int9]
# limits = [limits8, limits9]

if __name__=='__main__':
    # read the integration limits from the command line
    if len(sys.argv) > 1:
        limits = [(float(sys.argv[1]), float(sys.argv[2]))]
        ints = [integrand]

    for int, limit in zip(ints, limits):
        [value, nc]=integrate_driver(int,euler,limit[0],limit[1],tolerance,nstepmax,verbose)

        print 'Euler Integration for {} Converged to {} in {} steps'.format(int.__name__, value, nc)
        if (verbose):
            call(["mv","iterations.out","euler.out"])

        [value, nc]=integrate_driver(int,trapzd,limit[0],limit[1],tolerance,nstepmax,verbose)
        print 'Trapezoidal Integration for {} Converged to {} in {} steps'.format(int.__name__, value, nc)
        if (verbose):
            call(["mv","iterations.out","trapzd.out"])
        if True:
            [value, nc]=integrate_driver(int,midpoint,limit[0], limit[1],tolerance,nstepmax,verbose)
            print 'Midpoint Rule Integration Converged to ',value,' in ',nc,' steps'
            if (verbose):
                call(["mv","iterations.out","midpoint.out"])

        [value, nc]=simpson_driver(int,limit[0],limit[1],tolerance,nstepmax,verbose)
        print 'Simpson Integration for {} Converged to {} in {} steps'.format(int.__name__, value, nc)
        if (verbose):
            call(["mv","iterations.out","simpson.out"])
 
