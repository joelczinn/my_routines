from __future__ import print_function
__version__ = "0.7"
'''
(C) Joel C. Zinn 2018
zinn.44@osu.edu
0.8
--- added a routine that will give coordinates of a plate from leavitt's catalogue.
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
from scipy import stats, integrate
import numpy as np

def jhk2kepmag(J, H, K):
    '''
    from Huber+ 2016 ultimately from Howell+ 2012
    Inputs
     J : float or ndarray
     H : float or ndarray
     K : float or ndarray
    Outputs
     kepmag : color-based kepmag
    for things outside the vlaid color range, returns np.nan

    '''

    if not ((np.isscalar(J) and np.isscalar(H) and np.isscalar(K)) or
            (not np.isscalar(J) and not np.isscalar(H) and not np.isscalar(K))):
        raise Exception('J, H, & K must all either be scalars or all arrays')
    if np.isscalar(J):
        
        J = np.array([J])
        H = np.array([H])
        K = np.array([K])
    if not (len(J) == len(H) == len(K)):
        raise Exception('J H and K must all be the same size')
    x = J-K    
    dwarf = ~((J - H > 0.75)*(H - K > 0.1))
    
    kepmag = 0.42443603 + 3.7937617*x -  2.3267277*x**2 + 1.460255*x**3 + K
    kepmag[dwarf] = 0.314377 + 3.85667*x + 3.176111*x**2  -  25.3126*x**3 + 40.7221*x**4  -  19.2112*x**5 + K

    
    
    return kepmag



def g2tmag(G, BP, RP):
    '''
    from Stassun+ 2019
    Inputs
     G : float or ndarray
     BP : float or ndarray
     RP : float or ndarray
    Outputs
     tmag : color-based tmag
    for things outside the vlaid color range, returns np.nan
    Note that it is strictly valid only for deredenned photometry. Also there seems to be no distinctino between dwarf and giant, so this is not a very very valid transformation...
    '''

    if not ((np.isscalar(G) and np.isscalar(BP) and np.isscalar(RP)) or
            (not np.isscalar(G) and not np.isscalar(BP) and not np.isscalar(RP))):
        raise Exception('G, BP, & RP must all either be scalars or all arrays')
    if np.isscalar(G):
        
        G = np.array([G])
        BP = np.array([BP])
        RP = np.array([RP])
    if not (len(G) == len(BP) == len(RP)):
        raise Exception('G BP and RP must all be the same size')
    x = BP-RP
    good = (x > -0.2)*(x < 3.5)
    tmag = BP*np.nan
    tmag[good] = G - 0.00522555*x**3 + 0.0891337*x**2 - 0.633923*x + 0.0324473
    
    return tmag


def plate():
    '''
    this code is to get the positions of the Cepheids from the table in Leavitt 1912
    '''
    import astropy.units as u
    from astropy.coordinates import ICRS, SkyCoord, FK4
    coords = ["0 50 0.9 -73 07 00"]
    unit=(u.hourangle, u.deg)
    gc = SkyCoord(coords[0], frame=FK4, equinox='B1900', unit=unit)
    offsetra = 1606 - 12752
    offsetdec = 3218 - 10393


    offsetra = 13640 - 12752
    offsetdec = 10519 - 10393

    gc = gc.spherical_offsets_by(offsetra*u.arcsec, offsetdec*u.arcsec)
    print(gc)
    print(gc.transform_to(FK4(equinox='J2000')).ra.to_string(unit=u.hour))
    print(gc.transform_to(FK4(equinox='J2000')).dec.to_string(unit=u.degree, sep=('deg', 'm', 's')))


def dist(x, pdf, n):
    '''
    returns a random draw from the pdf that is passed. note that this involves discretizing the function along x, so 
    if it is necessary to get a certain resolution, make sure x is reflective of that. will generate random draws from the 
    center of the bins of x, in other words.
    Inputs
    x : ndarray
    pdf : ndarray
     the probability density, but it doesn't have to be normalized. this is done in dist()
    n : int
     how many random draws to return
     
    Outputs
    rand : ndarray(n)
     random draws from the pdf.
     
    '''
    if len(x) <= 1:
        raise Exception("x needs to have at least 2 entries")
    if len(x) != len(pdf):
        raise Exception("x and pdf have to have the same length")

    norm = integrate.cumulative_trapezoid(pdf, x)
    
    p = np.hstack((norm[0],np.diff(norm)))/norm[-1]
    
    custm = stats.rv_discrete(name='custm', values=(np.arange(len(x)-1), p))
    # return random draws from the center of the bin
    inds = custm.rvs(size=int(n))
    rand = x[inds+1] + np.diff(x)[inds]/2.
    return rand
        # for testing
        # import pylab as plt
        # masses = np.array([1.1,1.2,1.3,1.4,1.5])
        # pdf = np.array([0.1,0.2,1.0,2.0,0.1])
        # ndown = 20
        # nup = 10000
        # m_min = 0.9
        # m_max = 3.0
        # masses = np.hstack((np.linspace(0,m_min, ndown), np.linspace(m_min, m_max, nup)))
        # pdf = np.hstack((np.zeros(ndown), np.ones(nup)))
        # rand = dist(masses, pdf, 1000000)
        # plt.clf()
        # plt.hist(rand, bins=50)
        # plt.show()

def process_one(file, verbose=False, kepler=False, drop_duplicates=True, match_on='epic', sep=None, drop_no_epic=True, make_index=False):
    '''
    read in the table and rename the file column to make it 'EPIC' so that this can be later matched on.
    Inputs
    [ make_index : bool ]
     If True, will make epic to be the index of the dataframe.Default False. NOT COMPATIBLE YE TWTH COMBINE_TABLES, SO IT'S CURRENTLY NOT EVEN POSSIBLE FOR THIS TO BE DONE FROM A CALL TO COMBINE_TABLES. WOULD HAVE TO MANUALLY CHANGE THE COMBINE_TABLES CODE TO MAKE THAT POSOSIBLE. !!!
    [ drop_no_epic : bool ]
     if True, will not return rows for which no EPIC could be found.
    
    '''
    # put the whitespace separator last in case there are issues with some columns not being filled in
    # JCZ 221017
    # added nested try/except structure
    if sep is None:
        try:
            one = pd.read_table(file, comment='#', sep='\s+', header=0)
            if len(one.keys()) > 1:
                if verbose:
                    print('{} read successfully'.format(file))
            elif len(one.keys()) == 1:  
                one = pd.read_table(file, comment='#', sep=',', header=0)
            elif len(one.keys()) == 1:
                one = pd.read_table(file, comment='#', sep='\|', header=0)
            else:
                print('{} could not be read using , | or whitespace delimiters'.format(file))

        except:
            try:
                one = pd.read_table(file, comment='#', sep=',', header=0)
                if len(one.keys()) > 1:
                    if verbose:
                        print('{} read successfully'.format(file))
                elif len(one.keys()) == 1:
                    one = pd.read_table(file, comment='#', sep='\|', header=0)
                else:
                    print('{} could not be read using , | or whitespace delimiters'.format(file))
            except:
                one = pd.read_table(file, comment='#', sep='\|', header=0)
                if len(one.keys()) > 1:
                    if verbose:
                        print('{} read successfully'.format(file))
                else:
                    print('{} could not be read using , | or whitespace delimiters'.format(file))

            
                
    else:
        try:
            one = pd.read_table(file, comment='#', sep=sep, header=0)
        except:
            print('{} could not be read using the \'{}\' delimiter'.format(file,sep))
    # JCZ 070317
    # added this as an option to be kepler-friendly or not. jie's KICs don't prepend 00, so there are issues...
    
    if kepler:
        func = lambda x : (re.search('([0-9]{6,10})', x).group(1))
    else:
        func = lambda x : (re.search('([0-9]{6,10})', x).group(1))

    if verbose:
        print('found the following keys:')
        print(one.keys())
    if match_on == 'epic' or match_on == 'epic_orig':
        found = False
    else:
        found = True
    # JCZ 120317
    # aded kepid to deal with MAST files without having to rename the KIC col
    # JCZ 071118
    # 'epic' HAS TO BE FIRST !!!

    for k in ['epic', 'file', 'filename', 'EPIC', 'KICID', 'ID', 'KIC', 'KID', 'files', 'kepid', 'epic_orig','Kepler_ID']:
        if k in one.keys() and found is False:
            # now cut on the pipeline_rating:
            one['__epic_orig'] = one[k].copy()
            try:
                one[k] = one[k].astype(int)
            except:
                pass
            if ((one[k]).dtype != 'int64' and (one[k]).dtype != 'float64') and found is False:
                
                one[k] = one[k].apply(func).astype(int)
            # merge one and pipeline output table (flags_df) so that if they are different sizes, one = one[one['pipeline_rating'] > 0] will work
            # JCZ 071118
            # added k != 'epic', because otherwise there will be two epic columns
            if found is False and k != 'epic':
                one = one.rename(columns={k:'epic'})
            # JCZ 040317
            # added because will have two epic_orig cols if epic and epic_orig are already present in table
            if 'epic_orig' in one.keys():
                one = one.rename(columns={'epic_orig':'epic_orig_orig'})
            one = one.rename(columns={'__epic_orig':'epic_orig'})
            print(one.keys())

            found = True
    if not found:
        if verbose:
            print('couldnt find a column named file epic or filename so assuming the first unnamed column is the correct one. assuming no header in the input file')
        one = pd.read_table(file, comment='#', sep='\s+|,', engine='python', header=None)
        one = one.convert_objects(convert_numeric=True)
        print(one)

        one['epic_orig'] = one[0].copy()
        if not((one[0]).dtype == 'int64' or one[0].dtype == 'float64'):
            one[0] = one[0].astype(str).apply(func).astype(int)
        if one[0].dtype == 'float64':
            one[0] = one[0].astype(int)
        one = one.rename(columns={0:'epic'})
        one['epic'] = one['epic'].__array__().astype(int)
    # except:
        
    #     if 'epic' not in one.keys():
        
    #         print 'could not find file or EPIC column for {}...'.format(file)
    #     # if epic does exist as a field then it must be converted to an int
    #     else:
    #         print 'found epic field -- converting to int'
    #         one['epic'] = one['epic'].astype(int)
    # JCZ 290916
    # drop any duplicates.

    if drop_duplicates:
        if verbose:
            print('dropping duplicates')
        one = one.drop_duplicates(subset=[match_on])
    # JCZ 151117
    # drop any objects that don't have an EPIC
    if drop_no_epic:
        one = one.dropna(subset=[match_on])

        # print one[match_on]
    if one[match_on].dtype == 'float64':
        one[match_on] = one[match_on].astype(int)
    if make_index:
        one.set_index(match_on, inplace=True, drop=False)
    return one

    
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

def lonlatellb(lon, lat):
    '''
    converts into ecliptic lon and lat into galactic ell and b
    '''
    from astropy.coordinates import ICRS, Galactic, FK4, FK5, HeliocentricTrueEcliptic
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    c = SkyCoord(frame=HeliocentricTrueEcliptic, lon=lon*u.degree, lat=lat*u.degree, distance=1000*u.kpc)
    c = c.transform_to(frame=Galactic)
    return c.l.value, c.b.value

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
