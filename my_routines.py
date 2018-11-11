from __future__ import print_function
__version__ = "0.7"
'''
(C) Joel C. Zinn 2018
zinn.44@osu.edu
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
0.6
-- added plot_grid()
0.7
-- added radec2ellb(), added fits2pd
'''

import pyfits
import pandas as pd
import numpy as np
import os
from scipy.stats import ks_2samp

# JCZ 111118
# commented this out
# from SYDOSU.smooth_c import smooth_c

import matplotlib as mpl
import pylab as plt
# JCZ 311018
# commented this out.
# from SYDOSU.util import process_one, read_ts, freq_cut
import glob

class empty(object):
    def __init__(self):
        pass

def fits2pd(fitsfile):
    '''
    given an inputs fits file, will return pandas dataframe of that data.
    Inputs
    fitsfile : str
     fits file with at least a column called <match_on> (default 'epic'). Assumes only has on HUD
    Outputs
    fits_df : pd DataFrame
     pd DataFrame
    '''
    fits = pyfits.open(fitsfile)[1].data
    dict = {}
    for key in fits.names:
        # JCZ 070817
        # should not affect normal operations, but will successfully add data with more than one dimension.
        entry = fits.field(key).byteswap().newbyteorder()
        entry_dim = len(entry.shape)
        n_dim = entry.shape[np.argsort(entry.shape)[0]]
        dim_suffixes = np.arange(n_dim).astype(str)
        # JCZ 251017
        # for some reason now that i have added lengths to string columns, there are some that have no dimension? so need to skip over those...
        try:
            dim_suffixes[0] = ''
        except:            
            continue

        # add the dimensions one-by-one, along the short dimension -- for the case of two dimensions. for more dimensions, the array is simply flattened
        if entry_dim > 1 and entry_dim < 3:

            for dim, suffix in zip(range(n_dim), dim_suffixes):
                print('could not add {} to the DataFrame because it had multiple dimensions. breaking it up in to {} entries, instead...'.format(key, np.min(entry.shape)))
                print('adding entry {}'.format(key+suffix))
                if np.argsort(entry.shape)[0] == 0:
                    dict[key+suffix] = entry[dim, :]
                else:
                    dict[key+suffix] = entry[:, dim]
        elif entry_dim == 1:
            dict[key] = entry
    fits_df = pd.DataFrame(dict)
    return fits_df

def radecellb(ra, dec):
    '''
    converts ra and dec into ell and b
    '''
    from astropy.coordinates import ICRS, Galactic, FK4, FK5
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    c = SkyCoord(frame=ICRS, ra=ra*u.degree, dec=dec*u.degree)
    c = c.transform_to(frame=Galactic)
    return c.l.value, c.b.value

def radeclonlat(ra, dec):
    '''
    converts ra and dec into ecliptic lon and lat
    '''
    from astropy.coordinates import ICRS, Galactic, FK4, FK5, HeliocentricTrueEcliptic
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    c = SkyCoord(frame=ICRS, ra=ra*u.degree, dec=dec*u.degree, distance=1000*u.kpc)
    c = c.transform_to(frame=HeliocentricTrueEcliptic)
    return c.lon.value, c.lat.value

def plot_grid(funcs, args, kwds, labels, xlabel, ylabel, n_rows=3, n_cols=3):
    '''
    creates a grid of plots
    '''
    try:
        mpl.style.use('jcz_paper_latex')
    except:
        print('could not load jcz_paper_latex style file. will not use it.')

    j = 0    
    i = 0
    row_index = 0
    count = 0
    for _i in xrange(len(funcs)):
        for j in xrange(n_cols):
            if j == 0 and row_index == 0:
                fig, axes = mpl.pyplot.subplots(ncols=n_cols, nrows=n_rows, figsize=(12, 13), sharey='row', sharex='col')
                fig.subplots_adjust(hspace=0.0, wspace=0.0)        
                axes = axes.flatten()
                
            if j == n_cols//2:
                _xlabel = True
                _ylabel = False
            else:
                _xlabel = False
                _ylabel = False
            if j == 0:
                _xlabel = False
                _ylabel = False
            if row_index == 1 and j == 0:
                _ylabel = True
                
            if _ylabel:
                kwd['ylabel'] = ylabel
            else:
                kwd['ylabel'] = ''
            if _xlabel:
                kwd['xlabel'] = xlabel
            else:
                kwd['xlabel'] = ''
            kwd['ax'] = axes[i]
            funcs[_i](*args[_i], **kwds[_i])
            axes[i].text(0.6, 0.8, labels[_i], fontsize=15, weight='bold', transform=axes[i].transAxes)
            
            if row_index > 0 and row_index < n_rows:
                yticks = axes[i].yaxis.get_major_ticks()
                yticks[-1].label1.set_visible(False)
            else:
                yticks = axes[i].yaxis.get_major_ticks()
                yticks[0].label1.set_visible(False)

            i += 1
            if (row_index == n_rows-1 and j == n_rows-1) or (_i == len(funcs)-1 and j == n_cols-1):
                row_index=0
                i = 0
                j = 0
                plt.savefig('three_spec_'+str(count)+'.png', format='png')
                plt.clf()
                count += 1

            elif j == n_cols-1:
                j = 0
                row_index += 1
            else:
                j += 1




        




    
def ks_test(a,b):
    '''
    gives per cent confidence at which the two distributions, a and b, are inconsistent. Lower means more different and higher means more similar
    '''
    return 100. - (ks_2samp(a, b))[1]*100.

def str_ind(array, st, get_bool=False):
    '''
    returns the indices in <array> that contains <st>                                                                     
    Inputs                                                                                                                
    array : ndarray                                                                                                       
    st : str                         
    [ get_bool : bool ]
     if True, returns a boolean array of the same length as <array>, instead of an index array with a length not necessarily equal to the length of <array>. Default False.
    Outputs                                                                                                               
    idx : ndarray                                                                                                         
    containing indices corresponding to <array> that have <st> in them                                                    
    Notes
     will interpret elements in the passed array that are not strings as strings.
    '''

    if get_bool:
        indices = [True if st in str(x) else False for i,x in enumerate(array)]
    else:
        indices = [i for i,x in enumerate(array) if st in str(x)]
    # JCZ 160118
    # changed this to be a numpy array instead of list
    return np.array(indices)

def add_header(header, file):
    '''
    adds a header to an existing csv / ASCII file.
    Inputs
    header : str
     header that DOES NOT END IN \n -- this will be added automatically by add_header().
    file : str
     file that will be re-written
    
    '''
    with open(file, "r+") as f:
        #Read complete data of CSV file
        old = f.read()
        #Get cursor to start of file
        f.seek(0)
        #Write header and old data to file.
        f.write(header+ "\n" + old)
    return

import types
def fullprint(*args, **kwargs):
  from pprint import pprint
  import numpy
  opt = numpy.get_printoptions()
  numpy.set_printoptions(threshold='nan')
  pprint(*args, **kwargs)
  numpy.set_printoptions(**opt)

def add_errors(combined, combined_base, use_g_to_v=False, xs=None, es=None):
    '''
    add errors to the quantities according to their formal errors, assuming Gaussian errors
    Inputs
    combined : array that is modified in place by adding errors to columns specified in <xs> and <es>, based on mean values from combined_base
    combined_base : pd df
     for the keys in <xs>, <combined> columns will be replace by combined[xs] = combined_base[xs] + np.random.normal(loc=0.0, scale= combined_base[es])
    Outputs
    combined
    '''
    # add errors  according to the reported errors for MC analysis
    # things to add errors to
    # JCZ 201016
    # perturb JHK, as well...
    if use_g_to_v and (xs is None and es is None):
        xs = ['BTmag', 'VTmag', 'gmag', 'rmag', 'imag', 'numax', 'dnu', 'TEFF_COR', 'A_V', 'FE_H_ADOP_COR', 'jmag', 'hmag', 'kmag']
        # their error columns
        es = ['e_BTmag', 'e_VTmag', 'e_gmag', 'e_rmag', 'e_imag', 'numax_sig', 'dnu_sig', 'TEFF_COR_ERR', 'e_A_V', 'FE_H_ADOP_COR_ERR', 'e_jmag', 'e_hmag', 'e_kmag']
    elif (xs is None and es is None):
        xs = ['BTmag', 'VTmag', 'numax', 'dnu', 'TEFF_COR', 'A_V', 'jmag', 'hmag', 'kmag', 'FE_H_ADOP_COR', 'jmag', 'hmag', 'kmag']
        # their error columns
        es = ['e_BTmag', 'e_VTmag', 'numax_sig', 'dnu_sig', 'TEFF_COR_ERR', 'e_A_V', 'e_jmag', 'e_hmag', 'e_kmag', 'FE_H_ADOP_COR_ERR', 'e_jmag', 'e_hmag', 'e_kmag']
    
    N = len(combined)
    # JCZ 290317
    # don't need to copy it
    # _combined = combined.copy()
    for x,e in zip(xs, es):
        
        try:
            combined[x] = combined_base[x] + np.random.normal(loc=0., scale=combined_base[e], size=N)
        except ValueError:
            print('encountered ValueError for {} item.'.format(x))
    return combined

  
# def incremental_std(arr):
#     # incrementally calculated std by way of keeping track of 
#     # meansq : \sum_i \langle x_i \rangle ^2
#     # sqmean : \sum_i \langle x_i^2 \rangle
#     # and then use
#     # \sigma = \sqrt{sqmean - meansq}

    
#     std = []
#     meansq = 0.
#     sqmean = 0.
#     c = 0.
    
#     for j in range(n_bins):
#         # c += 1.
#         # meansq[j+offset] = meansq[j+offset]*(c-1.)/c + np.sum(amps[stride*j:-1:onewrap])/c
#         # sqmean[j+offset] = sqmean[j+offset]*(c-1.)/c + np.sum(amps[stride*j:-1:onewrap]**2)/c
#         pass    
#     std = np.sqrt(np.array(sqmean - meansq))
#     return 9999

def DM2d(DM):
    '''
    given a distance modulus, returns distance in pc
    '''

    return 10.**(1. + DM/5.0)

def DM_err2d_err(DM, DM_err):
    '''
    returns error on distance, given DM and DM_err
    Inputs
    DM : float
     distance modulus
    DM_err : float
     error on distance modulus
    Outputs
    (d, d_err) : tuple of float
     distance and error on distance (pc)
    '''
    return np.log(10.0)/5.0*DM_err*DM2d(DM)
def d2DM(d):
    '''
    given a distance in pc, returns distance modulus
    '''

    return 5.*np.log10(d) - 5.
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

def delete_func(x, y=None, ind=False, fill=None, full_ind=False, f1=lambda x: False, f2=lambda x: False):
    '''
    Returns an array without the values that satisfy either of the functions passed, f1 and f2. Optionally can replace these entries with a fill factor, <fill>.
    Full ind : bool
     if true, returns a boolean array of the same shape as input array, to be contrasted with the ind option, which returns an int ndarray of
     length equal to or less than input array.
    f1 : func
     One of two functions that can be specified. If an element of x or y is True, then it will be deleted (or its index will not be in the result when ind or full_ind is True). Default lambda x: False (i.e., return everything).
    f2 : func
     One of two functions that can be specified. If an element of x or y is True, then it will be deleted (or its index will not be in the result when ind or full_ind is True). Default lambda x: False (i.e., return everything).

    Notes
    General version of delete_inf_nan with f1=np.isinf and f2=np.isnan
    '''
    if fill is not None:
        if y is not None:
            a = np.where((f1(y)) | (f2(y)))[0]
        else:
            a = np.where((f1(x)) | (f2(x)))[0]
        if ind:
            return a
        x[a] = fill
        
        if y is not None:
            y[a] = fill
            return (x, y)
        else:
            return x
        

    if y is not None:
        a = np.where((~f1(y)) & (~f2(y)))[0]
    else:
        a = np.where((~f1(x)) & (~f2(x)))[0]
    if ind:
        return a


    if full_ind:
        return f1(x) | f2(x)

    if y is not None:
        return (x[a], y[a])
    else:
        return x[a]


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
        print(w % 2 > 0)
        if not w % 2 > 0:
            w = w - 1
            print('changed w from {} to {}.'.format(w+1, w))
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

def perc2sigma(perc):
    '''                                                                                                               
    convert percent to sigma, assuming normality
    '''
    from scipy import stats
    return stats.norm.interval(perc/100., loc=0.0, scale=1.0)[1]

def sigma_to_percent(sigma):
    '''
    uses z-score to convert to fraction
    '''
    if 'norm' not in imports():
        import scipy.stats
        
    return (scipy.stats.norm.cdf(sigma)*2.) - 1.0

def sig_diff(x, y, xerr, yerr):
    '''
    returns significance of difference in units of sigma for two measurements with errors
    <xerr> and <yerr>, assuming Gaussianity.
    Inputs
    x : float
    y : float
    xerr : float
    yerr : float
    Outputs
    diff : float
     in units of \sigma = \sqrt{xerr^2 + yerr^2}
    '''
    diff = np.abs(x-y)/np.sqrt(xerr**2 + yerr**2)
    return diff
