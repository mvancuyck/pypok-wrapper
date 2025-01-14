import astropy.units as u 
import numpy as np
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.io import fits



def cmn2cb( map_b, map_cmn, cb, sigma_b, sigma_b_1):

    w  = where( map_b != None)
    m1 = map_b[w]
    m2 = map_cmn[w]

    h         = histogram( m1, bin=1.0d0, reverse_ind=R)
    cb        = double(h*0.0d0)
    sigma_b   = double(h*0.0d0)
    sigma_b_1 = double(h*0.0d0)
    for i in range(0, n_elements(h)-1):

        if(R[i+1]-R[i] == 2):
            cb[i]        = avg(    m2[ R[ R[i] : R[i+1]-1]])
            sigma_b_1[i] = stddev( m2[ R[ R[i] : R[i+1]-1]])
            sigma_b[i]   = stddev( m2[ R[ R[i] : R[i+1]-1]])/sqrt((R[i+1]-R[i]))

        if(R[i+1]-R[i] == 1):
            cb[i]        = m2[ R[ R[i]]]
            sigma_b_1[i] = !values.f_nan
            sigma_b[i]   = !values.f_nan
            
    return




def cls2map(l, clt, nx, ny, res, map_t, cu_t, k_map, k_mapx, k_mapy, l_cutoff=0,  no_k_map_reset = 0, no_cu_t_rest = 0, ampl = 0, seed = 0, force = 0 , fwhm_arcmin = 0, zero_index=None, index=None):

    lmap_rad_x = nx*res.to(u.rad)
    lmap_rad_y = ny*res.to(u.rad)
    
    #Generate gaussian amplitudes
    noise    = np.zeros(( nx, ny)) 
    n2       = nx*ny

    noise = np.random.randint(0, 1, n2)  
    dmn1     = np.fft.fft2( noise)

    noise = np.random.randint(0, 1, n2)
    dmn2     = np.fft.fft2( noise)

    norm = np.sqrt( (nx/lmap_rad_x)*(ny/lmap_rad_y))

    if( no_k_map_reset is False):
        k_map, k_mapx, k_mapy = give_mape_k(res, np.zeros((nx,ny))) 
        
        k_map  = k_map  * 2*np.pi
        k_mapx = k_mapx * 2*np.pi
        k_mapy = k_mapy * 2*np.pi


    if(l_cutoff is not None): lmax = l_cutoff
    else:
        if(index is not None or zero_index is not None):lmax = k_map.max()
        else: lmax = l.max()


    if(fwhm_arcmin is not None):
        sigma = fwhm_arcmin.to(u.rad) * gaussian_fwhm_to_sigma
        bl    = exp(-l*(l+1)*sigma**2)
    else:
        bl = 1

    return k_map, k_mapx, k_mapy



def mc_reduce(tab_res, quick = False):
    
    nmc = tab_res.shape[0]
    nbins = tab_res.shape[1]
    tab_avg = np.mean(tab_res, axis = 0)
    
    sigma_avg = []
    for i in range(0, nmc):
        sigma_avg.append(np.stddev(tab_res[i,:])/np.sqrt(nmc))
    sigma_avg = np.asarray(sigma_avg)

    if(quick = False):
        cov_mat = np.zeros((nbins,nbins))
        xcor = np.zeros((nbins, nbins))
        #Cmpute covariance matrix
        for b in range(0, nbins):
            for b1 in range(0, nbins):
                cov_mat[b,b1] = np.mean((tab_res[b,:] - tab_avg[b]) * (tab_res[b1,:]-tab_avg[b1]))
                xcorr[b,b1] = cov_mat[b,b1]/np.sqrt(cov_mat[b,b]*cov_mat[b1,b1])
        #Cross correlation
        for b in range(0, nbins):
            for b1 in range(0, nbins):
                xcorr[b,b1] = cov_mat[b,b1]/np.sqrt(cov_mat[b,b]*cov_mat[b1,b1])
        return tab_avg, sigma_avg, cov_mat, xcorr
    else: return tab_avg, sigma_avg
