import numpy as np
from astropy.io import fits
import astropy.units as u
from pythonized_ipoker import *
from IPython import embed 
from astropy import wcs

#Charging the mask to apply on the data map
file_mask = "mask_raw.fits" #mask_apod.fits"
mask = fits.getdata(file_mask) 

#Charging the data map
map_file = "dust_map.fits"
data_map = fits.getdata(map_file)

#Getting the resolution from the data's header
res =  2 * u.arcmin #ATTENTION 5 * u.arcsec #w.wcs.cdelt[0] * u.deg

#Dust power spectrum
index     = -3
beta      = 3 #binning index

#Charging the noise map
measured_noise = fits.getdata("/home/nponthieu/Soft/Data/dust_noise.fits") #noise of the data

#Charging the multipols l and the corresponding values of input dust power spectrum and of the input noise power spectrum
l         = fits.getdata("/home/nponthieu/Soft/Data/dust_l.fits")
clt       = fits.getdata("/home/nponthieu/Soft/Data/dust_clt.fits")
clnoise   = fits.getdata("/home/nponthieu/Soft/Data/dust_clnoise.fits")

embed()
delta_l_over_l = 0.1 #For a logarithmic binning
remove_1st_bin = True #to improve Mbb's conditioning
include_noise  = True #to add noise to the simulated datanorm
#keep average of the data is not set by default.
#However, for MC simulation of the noise, the average of the noise map is kept.
keep_avg = False
nmc = 500 #Number of random realizations used to determine the noise pseudo-power spectrum and the covariance matrix
parfile = "poker.par" #ASCII params file for F90
scale = 1.5 #2 #Size of the zero padded large patch that contains the oberved patch
#If the data maps (two in case of cross correlation) are smoothed by the beam, 
#Then the beam information shoud be pass to the F90 to be included in the mixing matrix computation.
#In case of cross correlation, the mask is the combination of the masks of each map to be cross correlated
a_beam = np.ones(( int(data_map.shape[1] * scale) ,int(data_map.shape[1] * scale)))
b_beam = np.ones(( int(data_map.shape[1] * scale) ,int(data_map.shape[1] * scale)))


#Initialization
mask_large,  k_nyquist, kmin,dk_min, k_map, k_bin_tab, k_bin_width, map_k_binning, map_binned, xbintab, x_mtt_bb_m1, k_out = poker_initialization(parfile, data_map, mask, res, scale, beta, a_beam, b_beam, delta_l_over_l, 3, keep_avg = keep_avg, remove_1st_bin = remove_1st_bin  )

data_map_large = enlarges_masks_apodizes(data_map, mask, scale, window = mask_large, keep_avg = keep_avg)

input_data_pseudo_pk, input_data_pk_out, input_sigma_pk_out = ipoker(data_map_large, res, k_out, beta, map_binned, map_k_binning, x_mtt_bb_m1,  remove_1st_bin = remove_1st_bin )

#------------ Noise ----------------
if(include_noise):
    #Compute the noise pseudo-spectrum <N^>
    llcb_res = []

    #init cu_t_noise
    noise_map_large, cu_t_noise = cls2map( l, clnoise,  data_map_large.shape[1],  data_map_large.shape[0], res,  l_cutoff = k_nyquist)  



    for imc in range(0, nmc):
        print(imc)
        noise_map_large, cu_t_noise = cls2map( l, clnoise,  data_map_large.shape[1],  data_map_large.shape[0],  res, cu_t = cu_t_noise,  l_cutoff = k_nyquist)  

        noise_map_large = noise_map_large * mask_large

        pseudo_pk, pk_out, sigma_pk_out = ipoker(noise_map_large, res, k_out, beta, map_binned, map_k_binning, x_mtt_bb_m1,  remove_1st_bin = remove_1st_bin )

        llcb_res.append( pseudo_pk )

    llcb_res = np.asarray(llcb_res)
    pseudo_noise_powspec, sigma_avg = mc_reduce(llcb_res, quick = True)


else:
    clnoise           = np.zeros(clnoise.shape)
    cu_t_noise        = np.zeros(mask_large.shape)
    pseudo_noise_powspec = np.zeros((1))


#---------------------------------- Data power spectrum ------------------------------------------
if(include_noise): total_data_map = data_map + measured_noise
else: total_data_map = data_map


total_data_map_large = enlarges_masks_apodizes(total_data_map, mask, scale, window = mask_large, keep_avg = keep_avg)
data_pseudo_pk, data_pk_out, data_sigma_pk_out = ipoker(total_data_map_large, res, k_out, beta,  map_binned, map_k_binning, x_mtt_bb_m1, noise_pseudo_spec = pseudo_noise_powspec, remove_1st_bin = remove_1st_bin)


#-------- Error bars ----------------------

# Init signal cu_t
Map, cu_t = cls2map( l, clt,  data_map.shape[1] , data_map.shape[1], res,  l_cutoff = k_nyquist)

# Compute the theoretical binned input spectrum corresponding to cu_t
cb_t_in_th = cmn2cb( map_b, map_k_binning * cu_t )
if(remove_1st_bin): cb_t_in_th = cb_t_in_th[1:]

pk_res = []

for imc in range(0, nmc):
    #Signal
    map_t, cu_t = cls2map( l, clt, data_map.shape[1] , data_map.shape[1], res,  cu_t = cu_t, l_cutoff = k_nyquist)  
    map_t = enlarges_masks_apodizes(map_t, mask, scale, window = mask_large, keep_avg = keep_avg)

    #Add white noise
    noise, cu_t_noise =  cls2map( l, clnoise, data_map_large.shape[1] , data_map_large.shape[1], res_arcmin, cu_t = cu_t_noise,    l_cutoff = k_nyquist)  
    noise = noise * mask_large

    map_t = map_t + noise
    #Power spectrum
    pseudo_pk, pk_out, sigma_pk_out = ipoker(map_t_large, res, beta,  map_b, map_k_binning, x_mtt_bb_m1, noise_pseudo_spec = pseudo_noise_powspec, remove_1st_bin = remove_1st_bin)
    pk_res.append( pk_out )

pk_final, sigma_pk_final, cov_mat, xcorr = mc_reduce(pk_res)

#---Plot---

print("end")
