"""
integrate_sub.py -- subroutine functions for numerical integration
  integrate_driver does adaptive stepping and can be called with euler,
    trapzd, or midpoint as the underlying method
"""

import numpy as np

def integrate_driver(func,integrator,a,b,tolerance,nstepmax,verbose):
    """
    Integrate a function func() using the specified integrator routine
    integrator = euler, euler_loop, trapzd, or midpoint
    a = lower limit of integration
    b = upper limit of integration
    tolerance = fractional convergence required for integral
    nstepmax = maximum number of steps allowed
    verbose = 1 -> write individual iterations to "iterations.out"

    Number of steps starts at 4 and doubles until convergence or nstep>nstepmax
    """

    if (verbose):
        f=open("iterations.out","w")
    nstep=4
    oldint=0.0    
    integral=integrator(func,a,b,nstep)

    while ((np.fabs(oldint/integral-1.0) > tolerance) and (2*nstep<nstepmax)):
        oldint=integral
        nstep*=2
	integral=integrator(func,a,b,nstep)
        if (verbose):
	    hstep=(b-a)/nstep
            outstring="%8d %.8g %.8g\n" % (nstep,hstep,integral)
            f.write(outstring)

    if (verbose):
        f.close()
    if (np.fabs(oldint/integral-1.0) > tolerance):
        print "Warning, fractional convergence is only ", \
	  np.fabs(oldint/integral-1.0)
    return [integral, nstep]

def simpson_driver(func,a,b,tolerance,nstepmax,verbose):
    """
    Integrate a function func() using the simpson rule
    a = lower limit of integration
    b = upper limit of integration
    tolerance = fractional convergence required for integral
    nstepmax = maximum number of steps allowed
    verbose = 1 -> write individual iterations to "iterations.out"
    Number of steps starts at 4 and doubles until convergence or nstep>nstepmax
    """

    if (verbose):
        f=open("iterations.out","w")
    nstep=4
    oldint=0.0    
    integral=(4.*trapzd(func,a,b,nstep) - trapzd(func,a,b,nstep//2))/3.

    while ((np.fabs(oldint/integral-1.0) > tolerance) and (2*nstep<nstepmax)):
        oldint=integral
        nstep*=2
	integral = (4.*trapzd(func,a,b,nstep) - trapzd(func,a,b,nstep//2))/3.
        if (verbose):
	    hstep=(b-a)/nstep
            outstring="%8d %.8g %.8g\n" % (nstep,hstep,integral)
            f.write(outstring)

    if (verbose):
        f.close()
    if (np.fabs(oldint/integral-1.0) > tolerance):
        print "Warning, fractional convergence is only ", \
	  np.fabs(oldint/integral-1.0)
    return [integral, nstep]

def euler_loop(func,a,b,nstep):
    """
    Evaluate [\int_a^b func(x) dx] using Euler rule with nstep steps
    Use loop analogous to C or fortran
    """
    hstep=(b-a)/nstep
    y=a                                
    integral=hstep*func(y)
    for i in xrange(nstep-1):
        y+=hstep
        integral+=func(y)*hstep
    return(integral)

def euler(func,a,b,nstep):
    """ 
    Evaluate [\int_a^b func(x) dx] using Euler rule with nstep steps
    Use numpy array operations
    """
    hstep=(b-a)/nstep
    x=np.linspace(a,b-hstep,nstep)
    y=func(x)*hstep
    return (np.sum(y))


def trapzd(func,a,b,nstep):
    '''
    Evaluate [\int_a^b func(x) dx] using trapezoidal rule with nstep steps
    Use numpy array operations
    '''
    hstep = (b - a)/nstep
    y = 0.5*func(a)*hstep
    y += euler(func,a+hstep,b,nstep-1)
    y += 0.5*func(b)*hstep
    return y

def midpoint(func,a,b,nstep):
    """
    Evaluate [\int_a^b func(x) dx] using midpoint rule with <nstep> steps
    Use numpy array operations
    """
    hstep = (b-a)/nstep
    x = np.linspace(a, b-hstep, nstep) + hstep/2.
    y = func(x)*hstep
    return np.sum(y)
