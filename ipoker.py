import astropy.units as u 
import numpy as np
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.io import fits
from IPython import embed
from scipy import signal
from scipy import interpolate
import datetime
import scipy.constants as cst
import subprocess
from scipy.ndimage import gaussian_filter
import powspec
from scipy.ndimage import gaussian_filter1d
from scipy.ndimage import generic_filter
import astropy.convolution as conv
import matplotlib.pyplot as plt
import argparse 
from set_multipols import *
from progress.bar import Bar
from mbb import compute_mbb


def ipoker(data_map_large, res, beta, l_out, l_bintab, map_l_power_beta, l_map, M_bb_m1, noise_pseudo_pk = None, map_2_large =None):

    """
    Evaluate the pseudo power spectrum and the angular power spectrum estimate of the map given a mixing matrix and an estimation of the noise pseudo power spectrum. 

    pars: the parameters dictionnary
    data_map_large (astropy.Quantity) : the enlarged and masked map containing the data from which we want the power spectrum 
    l_out (astropy.Quantity): the l values at which we estimate the power spectrum
    map_l_power_beta (astropy.Quantity): map of l values to the power beta
    bins_mask: map indicating on which bins a pixel belong
    M_bb_m1: the inverse of the mixing matrix computed with beams taken into account.
    M_bb_nobeam_m1: the inverse of the mixing matrix computed with ams not taken into account.
    noise_pseudo_pk (optional, astropy.Quantity): the estimation of the noise's pseudo power spectrum
    map_2_large(optional, astropy.Quantity): the second  enlarged and masked map for the cross correlation

    return: 

    pseudo_pk_mes (astropy.Quantity): the pseudo power spectrum of the map
    pk_out (astropy.Quantity): the estimated angular power spectrum
    """


    
    measured_cl_map = compute_p2( data_map_large, res, map_2_large)  
   
    pseudo_pk_mes, edges = np.histogram(l_map, bins = l_bintab, weights = map_l_power_beta * measured_cl_map)
    histo, edges = np.histogram(l_map, bins = l_bintab)
    pseudo_pk_mes /= histo

    if(noise_pseudo_pk is None): noise_pseudo_pk = np.zeros(len(pseudo_pk_mes))

    #Binned k^beta*P(k) output power spectrum
    pb_out = np.real(np.dot(M_bb_m1, pseudo_pk_mes - noise_pseudo_pk))
    #Effective power spectrum at bin average assuming the process is described by a
    #power law of index=beta
    pk_out = pb_out * l_out**(-beta)
    return pseudo_pk_mes, pk_out

def apodizing_a_map(data_map, mask, patch = None, iy0=None, ix0=None,  keep_avg = False):
    """
    Remove the mean of the Map in the non-masked regions

    Put the data map in the zero-padded window and apply the smoothed mask if any.
    
    Map (astropy.Quantity): the data map

    mask: the mask of the data

    patch: the padding array

    keep_avg (bool): Remove the mean of the Map in the non-masked. 
                     Default is False  
    
    return: 

    Map_large (astropy.Quantity): the enlarged data map, apodized by the window, without its mean if required
    """
    Map = data_map.copy()
    #subtract the mean of the map
    if(not keep_avg): Map -= Map[np.where(mask != 0)].mean()

    if(patch.shape == mask.shape): Map_large = Map * mask
    else: 
        Map_large = np.zeros(patch.shape)
        #Put the map in the large patch and apply the apodized mask
        Map_large[iy0:iy0+Map.shape[0], ix0:ix0+Map.shape[1]] = Map
        Map_large = Map_large * patch 
        
    return Map_large
        
def make_patch(mask, scale, ny, nx, ny_large, nx_large, apod_length = 0, kernel_type = "sinc_cos"):

    #Coordinates of references of the Map inside the large patch
    iy0 = int(( (scale - 1)*ny/2))
    ix0 = int(( (scale - 1)*nx/2))

    patch = np.zeros((ny_large, nx_large))
    #Put the mask in a large map.
    patch[iy0:iy0+ny, ix0:ix0+nx] = mask
    #If apodization is required by the user:
    if(apod_length > 0 ): patch = patch_apod(patch, apod_length, kernel_type = kernel_type, map_shape = (ny,nx), coord_0 = (iy0,ix0) )

    return patch

def patch_apod(window, apod_length, kernel_type = "sin_cos", map_shape = (100,100), coord_0 = (25,25) ):

    #Shape of the map inside the window
    ny = map_shape[0]
    nx = map_shape[1]

    #Coordinates of reference of the map inside the window
    iy0 = coord_0[0]
    ix0 = coord_0[1]

    #Set an odd size for the kernel
    kernelSize = int(2*apod_length+1)
    kernel = np.ones((kernelSize, kernelSize))
    
    y = np.arange(kernelSize)
    y_axis = np.ones(kernelSize)

    if(kernel_type == "sin_cos"):

        #Apodizing both the egdes and holes
        #with the function from  the POKER paper Ponthieu et al. 2018
        #minizing the amplitude in Fourier Space

        w = (y <= apod_length)
        y_axis[w] = y[w]/apod_length - 1/(2*np.pi)*np.sin(2*np.pi*y[w]/apod_length)
        w = ( (kernelSize-1-y) <= apod_length)
        y_axis[w] = (kernelSize-1-y[w])/apod_length - 1/(2*np.pi)*np.sin(2*np.pi*(kernelSize-1-y[w])/apod_length)

        for y in range(0, kernelSize):
            kernel[y,:] = kernel[y,:] * y_axis
        for x in range(0, kernelSize):
            kernel[:,x] = kernel[:,x] * y_axis

    if(kernel_type == "Han"):

        #Apodizing both the edges and the holes 
        #with a Han function

        for y in range(0, kernelSize):
            kernel[y,:] = kernel[y,:] * signal.hann(kernelSize)
        for x in range(0, kernelSize):
            kernel[:,x] = kernel[:,x] * signal.hann(kernelSize)

    if(kernel_type == "sincos_edges_only"):
        
        mask = window[iy0:iy0+ny,ix0:ix0+nx].copy()
        
        w = (y <= apod_length)
        y_axis[w] = y[w]/apod_length - 1/(2*np.pi)*np.sin(2*np.pi*y[w]/apod_length)
        w = ( (kernelSize-1-y) <= apod_length)
        y_axis[w] = (kernelSize-1-y[w])/apod_length - 1/(2*np.pi)*np.sin(2*np.pi*(kernelSize-1-y[w])/apod_length)

        for y in range(0, kernelSize):
            kernel[y,:] = kernel[y,:] * y_axis
        for x in range(0, kernelSize):
            kernel[:,x] = kernel[:,x] * y_axis

        window[iy0:iy0+ny,ix0:ix0+nx] = 1

    w1  = conv.convolve_fft(window, kernel )
    
    w1[:iy0, 0:] = 0
    w1[iy0+ny:, 0:] = 0
    w1[0:, :ix0] = 0
    w1[0:, ix0+nx:] = 0
    
    wp = np.where(w1[iy0:iy0+ny,ix0:ix0+nx] <= 0.9)
    
    w2 = w1.copy()
    w2[iy0:iy0+ny,ix0:ix0+nx][wp] = 0
    
    w3 = conv.convolve_fft(w2, kernel )
    
    w3[:iy0, 0:] = 0
    w3[iy0+ny:, 0:] = 0
    w3[0:, :ix0] = 0
    w3[0:, ix0+nx:] = 0
    
    w3[iy0:iy0+ny,ix0:ix0+nx] = w3[iy0:iy0+ny,ix0:ix0+nx] * window[iy0:iy0+ny,ix0:ix0+nx]
    
    final_window = np.abs(w3)

    if(kernel_type == "sincos_edges_only"): final_window[iy0:iy0+ny,ix0:ix0+nx] *= mask

    return final_window

def make_mask(ny, nx, nholes, radius_ratio):  

    mask = np.ones((ny, nx))

    radius = np.minimum(nx,ny) / radius_ratio
    list_yc = np.random.randint(0, ny, nholes)
    list_xc = np.random.randint(0, nx, nholes)
    for y in range(0,ny):
        for x in range(0,nx):
            for h in range(0,nholes):
                if( np.sqrt( (x-list_xc[h])**2 + (y-list_yc[h])**2 ) <= radius):
                    #periodic boundary conditions
                    #----
                    if( x <= nx-1 ): x=x
                    else: x=x-nx
                    if( y <= ny-1 ): y=y
                    else: y=y-ny
                    #----                    
                    mask[y,x] = 0
    return mask 

def tf_beam_map( l_map, sigma):
    return np.exp(-1 *l_map**2 * (sigma.to(u.rad).value**2)/2)

def pk_beam_map( l_map, sigma):
    return np.exp(-1 *l_map**2 * (sigma**2))

def set_pars(P, load_mbb_directly = False, bypass_mtt_bb = False, nnodes=1, partition = '' ):

    #Getting the resolution from the data's header
    hdr = fits.getheader(f"{P['map_path']}/{P['map_name']}")
    res = (hdr["CDELT1"] * u.arcmin).to(u.rad).value

    nx, ny = hdr["NAXIS1"], hdr["NAXIS2"]
    nx_large, ny_large = int(nx*P['scale']), int(ny*P['scale'])

    l_nyquist, l_min, dl_min, l_bins_tab, l_out, l_map, map_l_power_beta = set_l_infos(ny_large, nx_large, ny_large, res, beta=P['beta'], delta_l_over_l= P['dll'])

    if(P['mask_name'] is None):
        if(P['make_a_mask'] is False): mask = np.ones((ny,nx))
        else: mask = make_mask(ny, nx, P['n_holes'],P['hole_ratio'])
    else: mask = fits.getdata(f"{P['mask_path']}/{P['mask_name']}")

    iy0 = int(( (P['scale'] - 1)*ny/2))
    ix0 = int(( (P['scale'] - 1)*nx/2))
    
    if(P['scale'] == 1 and P['apod_length'] == 0): patch = mask
    else: patch = make_patch(mask, P['scale'], ny, nx, ny_large, nx_large, P['apod_length'], P['kernel_type'])

    if(P['sigma_beam'] >0):
        sigma = P['sigma_beam']
        a_beam = tf_beam_map(l_map, sigma*u.rad)**2
        pkg_a = (np.exp(-1 *l_map**2 *(sigma**2)/2))**2
        b_beam = np.ones(l_map.shape)
    else:
        a_beam = np.ones(l_map.shape)
        pkg_a = np.ones(l_out.shape)
        b_beam = np.ones(l_map.shape)

    #If no mask nor beam are provided, save computation time
    
    if( (P['sigma_beam'] <= 0) and (P['mask_name'] is None) and (P['apod_length'] == 0) and (P['scale'] == 1) ): P['bypass_mtt_bb'] = True
    
    P['mask_array'] = mask
    P["patch"] = patch
    P["beamA"] = a_beam
    P["beamB"] = b_beam
    P["date"] = str(datetime.datetime.now())
    P["res"] = res
    P["l_max"] = l_nyquist
    P["l_min"] = l_min
    P["dl_min"] = dl_min
    P["l"] = l_out
    P["l_bins"] = l_bins_tab
    P['map_l_power_beta'] = map_l_power_beta
    P["nx"] = nx
    P["ny"] = ny
    P["nx_large"] = nx_large
    P["ny_large"] = ny_large
    P["iy0"] = iy0
    P["ix0"] = ix0
    P["l_map"] = l_map

    #Get mode mixing matrix
    if(P['bypass_mtt_bb']):
        print("bypassé")
        M_bb = np.identity(len(l_bins_tab)-1)
    elif(load_mbb_directly == False): 
        #M_bb = compute_mbb(P) 
        #Save the needed arrays for the mixing matrix computation:
            f = fits.PrimaryHDU(patch)
            hdu = fits.HDUList([f])
            hdu.writeto(P["mask_file"],  overwrite = True)

            f = fits.PrimaryHDU(a_beam)
            hdu = fits.HDUList([f])
            hdu.writeto(P["beam_a_file"],  overwrite = True)

            f = fits.PrimaryHDU(b_beam)
            hdu = fits.HDUList([f])
            hdu.writeto(P["beam_b_file"],  overwrite = True)

            pars = set_pars(res, 
                            file_mtt_bb_x = P["file_mtt_bb_x"] ,  
                            scale = P["scale"], 
                            delta_l_over_l = P['dll'],
                            beta = P["beta"],
                            log_binning = int(np.float(np.bool(P['dll']))), 
                            nx = nx, ny = ny,
                            nx_large =  int(ny * P["scale"]), ny_large =  int(nx * P["scale"]), 
                            remove_1st_bin = int(np.float(P["remove_1st_bin"])),
                            keep_avg = int(np.float(0)),  
                            apod_length = P["apod_length"] )    
            params2ascii( pars )
        
            #Save the pars dictionnary, in that precise order, needed for the computation of the mixing matrix.
            print("begin poker_mbb")

            ####################################### 
            #Compute the mode mixing matrix.
            if(salloc=='salloc'): subprocess.run( ["salloc", "--ntasks-per-node=24", f"--nodes={nnodes}", "-p", f"{partition}" , "mpirun", "poker_mbb_mpi", P["outparfile"]])
            if(salloc=='sbash' ): subprocess.run( [ "mpirun", "poker_mbb_mpi", P["outparfile"]])

            print( "MBB COMPUTED. ")
            #######################################
            
            #Loads results:
            M_bb = fits.getdata( P["file_mtt_bb_x"] )


    else:  M_bb = fits.getdata(f"{P['mbb_path']}/{P['mbb_name']}")
            
    #Discard DC bin to improve M_bb's conditioning
    if(P['remove_1st_bin']):
        #Inverse matrix of M_bb' in equation 17
        M_bb = M_bb[1:,1:]
        M_bb_m1 = np.linalg.inv( M_bb )
        P["l"] = P["l"][1:]
        P["l_bins"] = P["l_bins"][1:]
    else: 
        M_bb_m1 = np.linalg.inv(M_bb)
    
    P["m_bb"] =  M_bb
    P["m_bb_m1"] = M_bb_m1

    return P

def compute_p2(data_map, res,  map2 = None):
    """
    Measure the power spectrum  in a map

    data_map (astropy.Quantity): the map in which to measure the power spectum
    res (astropy.Quantity): the resolution of the map
    map2 (optional, astropy.Quantity): the second map in case of cross correlation

    return:
    
    p2 (astropy.Quantity): the measured power spectrum
    """

    npix = data_map.shape[0] * data_map.shape[1]
    norm = res**2 / npix
    if(map2 is not None):
        ft  = np.fft.fft2( data_map)
        ft1 = np.fft.fft2( map2)
        p2  = 0.5*( (ft*np.conj(ft1)) + (np.conj(ft)*ft1) )  * norm
    else: 
        p2 = np.abs( (np.fft.fft2( data_map))**2) * norm

    return p2

def mc_reduce(tab_res, quick = False):
    """
    reduce the results of n Monte Carlo simulations.

    tab_res (astropy.Quantity): table containing the power spectra of n Monte Carlo simulation

    quick (bool): to return only mean and standard deviation.
    
    return:

    tab_avg (astropy.Quantity): the average

    sigma_avg (astropy.Quantity): the standard deviation 

    cov_mat (optional, astropy.Quantity): the covariance matrix

    xcor (optional, astropy.Quantity): the bin-bin correlation matrix
    """
    
    nmc = tab_res.shape[0]

    nbins = tab_res.shape[1]

    tab_avg = np.mean(tab_res, axis = 0)

    sigma_avg = np.std(tab_res, axis = 0) #/ np.sqrt(nmc)

    if(not quick):

        cov_mat = np.zeros((nbins,nbins))

        xcor = np.zeros((nbins, nbins))

        #Compute covariance matrix
        for b in range(0, nbins):
            for b1 in range(0, nbins):
                cov_mat[b,b1] = np.mean((tab_res[:,b] - tab_avg[b]) * (tab_res[:,b1]-tab_avg[b1]))

        #Cross correlation
        for b in range(0, nbins):
            for b1 in range(0, nbins):
                xcor[b,b1] = cov_mat[b,b1]/np.sqrt(cov_mat[b,b]*cov_mat[b1,b1])

        return tab_avg, sigma_avg, cov_mat, np.abs(xcor) 

    else: return tab_avg, sigma_avg

def load_params(path):
    #Load the parameter file as a dictionnary

    file = open(path)

    params = {}
    for line in file:
        line = line.strip()
        if not line.startswith("#"):
            no_comment = line.split('#')[0]
            key_value = no_comment.split("=")
            if len(key_value) == 2:
                params[key_value[0].strip()] = key_value[1].strip()

    for key in params.keys():
        params[key] = eval(params[key])

    return params

'''
def worker_nmc(Nmc):
    global _args
    l, cl_signal, ny, nx, iy0, ix0, res, cl_map_signal, l, map_l_power_beta, l_map, l_bintab, l_nyquist, patch, M_bb_m1, noise_pseudo_pk, keep_avg, remove_1st_bin = _args 
    Pk_list = []
    Pseudo_pk_list = []
    
    for N in Nmc:
        
        map_t, _ = gaussian_random_field(l, cl_signal, ny_large, nx_large, res, l_cutoff = l_max, ny_exit = ny, nx_exit = nx, sigma = sigma_beam, cl_map = cl_map_signal)
        map_t = apodizing_a_map(map_t, mask, patch = patch, iy0=iy0, ix0=ix0,  keep_avg = keep_avg )
        pseudo_pk, pk_out = ipoker(map_t, res, beta, l, map_l_power_beta, l_map, l_bins,M_bb_m1, remove_1st_bin = remove_1st_bin, noise_pseudo_pk = noise_pseudo_pk)
        
        Pk_list.append(pk_out.to(u.Jy**2/u.sr))
        Pseudo_pk_list.append(pseudo_pk)
        
    return (np.asarray(Pk_list), np.asarray(Pseudo_pk_list))


def Nmc_pll(nmc, l, cl_signal, ny, nx, iy0, ix0, res, beta, l, map_l_power_beta, l_map, l_bintab, l_nyquist, l_max, patch, M_bb_m1, sigma_beam = 0*u.rad, noise_pseudo_pk = None, keep_avg = False, remove_1st_bin = True, ncpus = 24):
    """
    POKER error bar computation using MC simulation
    """
    #beam_area = np.sum(kernel_channel.array)
    Nmc= np.arange(1,nmc+1,1)

    #Initiates cl_map_signal
    print("Initiate cl map")
    _, cl_map_signal = gaussian_random_field(l, cl_signal, ny,nx,res,beta,  l_cutoff= l_max,  sigma = sigma_beam)

    with Pool(ncpus, initializer=worker_init, initargs=(l, cl_signal, ny, nx, iy0, ix0, res, beta, cl_map_signal, l, map_l_power_beta, l_map, l_bintab, l_max, patch, M_bb_m1, noise_pseudo_pk, keep_avg, remove_1st_bin )) as p:
        # Transform full cube (nchan, npix, npix) as (npix*npix, nchan)
        pk_list = p.map(worker_nmc, np.array_split(Nmc, ncpus) )
    
    pk_out_list = np.vstack(pk_list) * u.Jy**2/u.sr
    pk_final, sigma_pk_final, cov_mat, corr_mat = mc_reduce(pk_out_list)
    return pk_final, sigma_pk_final 
#--- Worker for MC simulations
'''

def params2ascii( pars, name = None ):
    """
    Write the pars file
    """
    if(name is None):
        f = open(pars["outparfile"], "w")
    else: f = open(name, "w")
    for key in zip(pars.items()):
        f.write(str(key[0][0])+" = "+  str(key[0][1]) +"\n")
    f.close()
    return True
