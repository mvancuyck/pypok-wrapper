import numpy as np
import scipy.constants as cst
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u
from IPython import embed
import matplotlib.cm as cm
from scipy.optimize import curve_fit
from scipy import interpolate


def gaussian_random_field(l, clt, ny, nx, res, ny_exit = None, nx_exit = None, l_cutoff=None,  sigma = 0 * u.rad, cl_map = None , force = True):
    """
    Create a map of a Gaussian random field, given its angular power spectrum or from a power law with a given index 
    
    l (astropy.Quantity): the multipol

    clt (astropy.Quantity): the given power spectrum 

    ny, nx: the size of the map to generate

    res (astropy.Quantity): the resolution of the map

    ny_exit, nx_exit (optional): <ny, nx, size of the map to return, to avoid periodic boundaries

    l_cutoff (optional, astropy.Quantity): the maximum multipol to taken into account in the map generation

    sigma (optional, astropy.Quantity): the size of the instrumental Gaussian beam profile in sigma unit to apply, if any. 

    ampl (optional): the amplitude of the power law to compute the angular power spectrum

    index (optional): the index of the power law to compute the angular power spectrum

    cl_map (optional, astropy.Quantity): map containing the power spectrum values for each l values of the map. 
                                         If provided, not recomputed. 

    force (bool): set the negative values in the power spectrum map to zero. 

    return: 
    
    real_space_map (astropy.Quantity): the generated map in real space.

    cl_map (astropy.Quantity): the angular power spectrum map
    """


    lmap_rad_y = ny*res.to(u.rad)
    lmap_rad_x = nx*res.to(u.rad)
    #Generate gaussian amplitudes
    norm = 1/res.to(u.rad) #np.sqrt( (nx/lmap_rad_x)*(ny/lmap_rad_y)) 

    np.random.seed()

    noise = np.random.normal(loc=0, scale=1, size=(ny,nx))
    dmn1  = np.fft.fft2( noise )

    #Interpolate input power spectrum
    if(cl_map is None):

        l_map = give_map_spatial_freq(res.to(u.rad), ny, nx)

        if(l_cutoff is not None): lmax = np.minimum(l_cutoff.to(u.rad**-1).value, l.max().to(u.rad**-1).value) * u.rad**-1
        else: lmax = np.minimum(l.to(u.rad**-1).max().value, l_map.max().to(u.rad**-1).value) * u.rad**-1

        if(sigma.value != 0): bl = np.exp(-(l.to(u.rad**-1).value**2)*(sigma.to(u.rad).value**2)) #l(l+1) pour CMB
        else: bl = np.ones(l.shape)
        
        cl_map = np.zeros(l_map.shape)
        w = np.where((l_map>l.min()) & (l_map<=lmax))
        if(not w[0].any()): print("wrong k range")
        else:
            #Power law spectrum
            print("interpolate")
            f = interpolate.interp1d( l.to(u.rad**-1).value, clt.value*bl,  kind='linear')
            cl_map[w] = f(l_map[w])
            w1 = np.where( cl_map <= 0)
            if(w1[0].shape[0] != 0 and force): cl_map[w1] = 0
            cl_map = cl_map * clt.unit

    #Fill amn_t
    amn_t = dmn1 * norm * np.sqrt( cl_map )
    
    #Output map
    real_space_map = np.real(np.fft.ifft2( amn_t )).to(u.Jy/u.sr)
    
    if(ny_exit is not None and nx_exit is not None):
        scale_y = ny / ny_exit
        iy0 = int(( (scale_y - 1)*ny_exit/2))
        scale_x = nx / nx_exit
        ix0 = int(( (scale_x - 1)*nx_exit/2))
        return real_space_map[iy0:int(iy0+ny_exit), ix0:int(ix0+nx_exit)], cl_map

    else: return real_space_map, cl_map

def give_map_spatial_freq(res, ny, nx, output_unit = "l"):
    
    #res must be in rad!
    lmap_y = ny*res #rad
    lmap_x = nx*res #rad
    map_ky = np.float64(np.zeros((ny, nx)))
    map_kx = np.float64(np.zeros((ny, nx)))
    map_k  = np.float64(np.zeros((ny, nx)))
    for m in range(0,nx):
        if(m <= nx/2): m1 = m
        else: m1 = m - nx
        for n in range(0,ny):
            if(n <= ny/2): n1 = n
            else: n1 = n - ny
            kx = np.float64(m1/lmap_x) 
            ky = np.float64(n1/lmap_y)
            
            map_kx[n,m] = kx
            map_ky[n,m] = ky
            map_k[n,m] = np.float64(np.sqrt( kx**2 + ky**2))

    map_l = map_k * 2 * np.pi #* (res**-1).unit
    map_k = map_k 
    if(output_unit == "l"): return map_l
    elif(output_unit == "k"): return map_k

    
def set_l_infos(ny, nx, nmin_map, res, beta=1, delta_l_over_l = 0):
    #compute the maximum value of l
    l_nyquist = np.pi / res #rad**-1
    #compute the map of the radial l multipols
    l_map  =  give_map_spatial_freq(res, ny, nx) #rad-1
    #Measuring from the map of l the miminum value of l
    l_min = np.min(l_map[np.nonzero(l_map)]) #rad-1
    n = np.min((ny,nx))
    #Computing the minimum size of a l bin
    dl_min = 2 * 2 * np.pi / (nmin_map * res) #rad-1
    l_range = [l_min, l_nyquist] #rad**-1
    #Setting the binning of the l map
    l_bin_tab, l_bin_width = make_bintab(l_range, dl_min, delta_l_over_l) #rad-1 
    #Compute bin addresses
    l_out, edges = np.histogram(l_map, bins = l_bin_tab, weights = l_map)
    histo, edegs = np.histogram(l_map, bins = l_bin_tab)
    l_out = l_out / histo
    #compute the radial multipols to the power beta 
    map_l_power_beta = l_map.copy()
    map_l_power_beta[l_map !=0] = (l_map[l_map != 0]**beta)
    map_l_power_beta *= l_map**beta
    
    return l_nyquist, l_min, dl_min, l_bin_tab, l_out, l_map, map_l_power_beta


def make_bintab(l, delta_l_min, dll = 0, delta_l_max=0, ):
    lmax      = l[1]
    lmin      = l[0]
    if(dll == 0): 
        bintab = np.arange(lmin, lmax, delta_l_min )
        bin_width = np.ones(len(bintab)-1) * delta_l_min
    else:

        l1 = lmin
        delta_l = 0 
        bintab = []
        bintab.append(lmin)
        bin_width = []
        while(l1 + delta_l <= lmax):
    
            delta_l = np.minimum( np.maximum(l1*dll, delta_l_min) , (lmax - l1) )

            if(delta_l_max != 0): 
                delta_l = np.minimum(delta_l, delta_l_max)
            l1 = l1 +  delta_l
            bintab.append(l1)
            bin_width.append(delta_l)
        bintab = np.asarray(bintab)

    if( bintab.max() <= lmax):
        bintab[-1]    = lmax
        bin_width[-1] = bintab[-1]-bintab[-2]
    bintab[0] = bintab[0]*0.99
    bintab    = np.insert(bintab,0,0)

    return bintab , bin_width 
