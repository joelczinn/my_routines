
__version__ = "0.5"
'''
0.1
-- added imports() to check what and what has not been imported
-- then use this in new function sigma_to_percent which converts z-scores into percentages.
0.2
-- find_nearest_val and find_nearest_ind were switched
0.3
-- added pass statements so that when condorized, will behave
0.4
-- added smooth function
0.5
-- added fullprint from http://stackoverflow.com/questions/1987694/print-the-full-numpy-array
'''

import numpy as np
import os
#condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized from ipdb import set_trace as pause

class empty(object):
    def __init__(self):
        pass


import types
def fullprint(*args, **kwargs):
  from pprint import pprint
  import numpy
  opt = numpy.get_printoptions()
  numpy.set_printoptions(threshold='nan')
  pprint(*args, **kwargs)
  numpy.set_printoptions(**opt)

def incremental_std(arr):
    # incrementally calculated std by way of keeping track of 
    # meansq : \sum_i \langle x_i \rangle ^2
    # sqmean : \sum_i \langle x_i^2 \rangle
    # and then use
    # \sigma = \sqrt{sqmean - meansq}

    
    std = []
    meansq = 0.
    sqmean = 0.
    c = 0.
    
    for j in range(n_bins):
        # c += 1.
        # meansq[j+offset] = meansq[j+offset]*(c-1.)/c + np.sum(amps[stride*j:-1:onewrap])/c
        # sqmean[j+offset] = sqmean[j+offset]*(c-1.)/c + np.sum(amps[stride*j:-1:onewrap]**2)/c
            
    std = np.sqrt(np.array(sqmean - meansq))
    return 9999

def imports():
    '''
    returns imported modules
    '''
    for name, val in globals().items():
        if isinstance(val, types.ModuleType):
            yield val.__name__

# def strip_ext(file):
#     '''
#     Inputs
#     file : str
#      path that you want to strip of the extension (e.g., /this/is/cool.notcool)
#     Outputs
#     stripped filed : str
#      e.g., /this/is/cool
#     Notes
#     uses os.path to do this
#     does not work on multiple files!
#     '''
#     return os.path.join(os.
def ang2p(ang):
    return 2.*np.pi/ang

def p2ang(p):
    return 2.*np.pi/p

def f2ang(f):
    return f*np.pi*2.

def ang2f(ang):
    return ang/2./np.pi

def f2p(f):
    return 1./f


def p2f(p):
    return 1./p

def isinfnan(x):

    return np.logical_or((np.isinf(x)), (np.isnan(x)))
def delete_inf_nan(x, y=None, ind=False, fill=None, full_ind=False):
    '''
    Returns an array without the values that are NaN or Inf. Optionally can replace these entries with a fill factor, <fill>.
    full ind : bool
     if true, returns a boolean array of the same shape as input array, to be contrasted with the ind option, which returns an int ndarray of
     length equal to or less than input array.
    '''
    if fill is not None:
        if y:
            a = np.where((np.isinf(y)) | (np.isnan(y)))[0]
        else:
            a = np.where((np.isinf(x)) | (np.isnan(x)))[0]
        if ind:
            return a
        x[a] = fill
        
        if y:
            y[a] = fill
            return (x, y)
        else:
            return x
        

    if y:
        a = np.where((~np.isinf(y)) & (~np.isnan(y)))[0]
    else:
        a = np.where((~np.isinf(x)) & (~np.isnan(x)))[0]
    if ind:
        return a


    if full_ind:
        return np.isinf(x) | np.isnan(x)

    if y:
        return (x[a], y[a])
    else:
        return x[a]

def n_elements(array):
    '''
    Will count the number of elements along the zeroth dimension of a list or numpy array.
    '''
    if type(array) == type([]):
        return len(array)
    if type(array) == type(np.array([])):
        return array.shape[0]
    return 0

def plot(x, y=None):
    '''
    Automatically shows a plot, with interactive window on
    '''
    #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized import pylab as plt
    if y is not None:
        #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized plt.plot(x,y); plt.ion(); plt.show()
#condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized pause()
        pass
    else:
        #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized plt.plot(x); plt.ion(); plt.show()
        pass
#condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized #condorized pause()

def find_nearest_ind(array, value):
    '''
    taken from unutbu http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
    '''

    ind = (np.abs(array-value)).argmin()
    return ind

def find_nearest_val(array, value):
    '''
    taken from unutbu http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
    '''

    return array[find_nearest_ind(array, value)]

def smooth(raw, w, type='boxcar'):
    """
    Inputs
     w : float or int
      triangle:
       w is the FWHM of the triangle
      boxcar:
       w is the width of the boxcar - 2
     
    """

    N = raw.shape[0]

    sm = raw.copy()

    if type == 'boxcar':
        print w % 2 > 0
        if not w % 2 > 0:
            w = w - 1
            print 'changed w from {} to {}.'.format(w+1, w)
        inds = np.linspace(w, N-w, num=N-2*w+1)
        w = int(w)-2

        w_i = np.ones(w+2)

    if type == 'triangle':
        # base of the triangle in units of indices
        inds = np.linspace(w, N-w-1, num=N-2*w-1)
        b = w*2.
        w = int(b) - 1
        center = int((w+1)/2.)
        w_i = np.zeros(shape=(w+2))
        for i in range(w+1+1):

            w_i[i] = (b - 2*np.abs(center - i))/b

        
    for ind in inds:

        sm[ind] = np.sum(w_i*raw[ind - (w+1)/2 : ind + (w+1)/2+1])/np.sum(w_i)

        

    return sm


def sigma_to_percent(sigma):
    '''
    uses z-score to convert to percent
    '''
    if 'norm' not in imports():
        import scipy.stats.norm as norm
    return norm(sigma)*2. -1.
