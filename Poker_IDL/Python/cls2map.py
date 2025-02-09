import astropy.units as u 
import numpy as np
from astropy.stats import gaussian_fwhm_to_sigma


def cls2map(l, clt, nx, ny, res_arcmin, map_t, cu_t, k_map, k_mapx, k_mapy, zero_index=None, index=None, l_cutoff=0, no_k_map_reset = 0, no_cu_t_rest = 0, ampl = 0, seed = 0, force = 0 , fwhm_arcmin = 0):

    lmap_rad_x = nx*res_arcmin.to(u.rad)
    lmap_rad_y = ny*res_arcmin.to(u.rad)
    
    #Generate gaussian amplitudes
    noise    = dblarr( nx, ny)
    n2       = nx*ny

    noise = np.randomn( seed, n2)
    dmn1     = fft( noise, /double)

    noise[*] = randomn( seed, n2, /double)
    dmn2     = fft( noise, /double)

    # Init amn fields
    # nx and ny appear in the definition of norm because the amplitudes in Fourier space
    # are generated by the fft of a white noise map in real space
    norm = sqrt( (nx/lmap_rad_x)*(ny/lmap_rad_y))

    #Define k modes
    if(no_k_map_reset == None):
        give_map_k, res_arcminrcmin2rad, dblarr(nx,ny), k_map, k_mapx, k_mapy
        k_map  = k_map  * 2.0d0*!dpi
        k_mapx = k_mapx * 2.0d0*!dpi
        k_mapy = k_mapy * 2.0d0*!dpi

        if(l_cutoff != 0):
            lmax = l_cutoff
        else:
            if(index != None or zero_index != None):
                lmax =  k_map.max()
            else:
                lmax = l.max()

        if(fwhm_arcmin is not None):
            sigma = fwhm_arcmin.to(u.rad) * gaussian_fwhm_to_sigma
            dl    = double(l)
            bl    = np.exp(-dl*(dl+1)*sigma^2)
        else:
            bl = 1.0d0


    #Interpolate input power spectrum
    if ( no_cu_t_reset is None ):
        cu_t = k_map*0.0d0
        w = np.where( k_map gt 0 and k_map le lmax, nw)
        if(nw==0):
            print("wrong k range")
            stop
        else:
            #Power law spectrum
            if (ampl is None): ampl = 1.0d0
            if (index is not None or zero_index is not None):
            if(fwhm_arcmin is not None):
                sigma = fwhm_arcmin.to(u.rad) * gaussian_fwhm_to_sigma
                bk    = exp(-k_map*(k_map+1)*sigma^2)
            else:
                bk = 1.0d0

            if(zero_index is None): index = 0

            cu_t[w] = ampl * k_map[w]^index * bk[w]

        else:
            cu_t[w] = ampl * interpol( clt*bl, l, k_map[w])
            w1 = where( cu_t lt 0, nw1)
            if(nw1!= 0):
            if(force is not None):
                cu_t[w1] = 0.0d0
            else: 
               print( "cu_t < 0 for "+strtrim( nw1,2)+" values")
               stop
   
    #Fill amn_t
    amn_t = dmn1 * norm * sqrt( cu_t)
    
    #Output map
    map_t = double( fft( amn_t, /inverse, /double))

    return



