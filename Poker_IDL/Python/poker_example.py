import numpy as np
from astropy.io import fits
import astropy.units as u
from ipoker import *
from IPython import embed 
from astropy import wcs
#from functions import * 

#Mask

file_mask = "/home/nponthieu/Soft/Data/mask_raw.fits" #mask_apod.fits"
mask = fits.getdata(file_mask)


#Dust power spectrum
index     = -3
beta      = 3 #binning index
Map_file  = "/home/nponthieu/Soft/Data/dust_map.fits"
Map       = fits.getdata(Map_file)
noise     = fits.getdata("/home/nponthieu/Soft/Data/dust_noise.fits")
l         = fits.getdata("/home/nponthieu/Soft/Data/dust_l.fits")
clt       = fits.getdata("/home/nponthieu/Soft/Data/dust_clt.fits")
clnoise   = fits.getdata("/home/nponthieu/Soft/Data/dust_clnoise.fits")

delta_l_over_l = 0.1
log_binning = True

#to improve Mbb's conditioning
remove_1st_bin = True

#to add noise to the simulated data
include_noise  = True

#Number of random realizations used to determine the noise pseudo-power spectrum and the covariance matrix
nmc = 500

#Size of the zero padded large patch that contains the oberved patch
scale = 1.5 # 2

#keep_avg is not set by default to remove the average of the data from the map before power spectrum estimation
keep_avg = False

parfile = "poker.par" #ASCII params file for F90

#Map header
#Header = fits.getheader(Map_file) #probleme pour charger les wcs
#w = wcs.WCS(Header)

res =  2 * u.arcmin #ATTENTION FOIREUX 5 * u.arcsec #w.wcs.cdelt[0] * u.deg

nx_large = int( scale*Map.shape[0])
ny_large = int( scale*Map.shape[1])


k_out, pk_out, sigma_pk_out, k_nyquist= run_ipoker(parfile, Map, mask, res, scale, nx_large, ny_large, beta, delta_l_over_l , log_binning, remove_1st_bin, keep_avg ) 

"""

#------------ Noise ----------------
if(include_noise== True):

   #Add simulated noise to the simulated data
   data_map += noise

   #Compute the noise pseudo-spectrum
   llcb_res = []

   #Init cu_t_noise
   map_t, cu_t, k_mapp = cls2map( l, clnoise, nx, ny, res,  l_cutoff = k_nyquist)
   for imc in range(0, nmc):
      map_t, cu_t, k_mapp = cls2map( l, clnoise, nx, ny, res,  l_cutoff = k_nyquist)  
      k_out, pk_out, sigma_pk_out, k_nyquist= run_ipoker(data_map, res, mask,scale, beta, remove_1st_bin = remove_1st_bin, delta_l_over_l  )  #??
      
      llcb_res.append( pseudo_pk )

   pseudo_noise_powspec, sigma_avg = mc_reduce(llcb_res, quick = True)


else:
   clnoise           = np.zeros(clnoise.shape)
   cu_t_noise        = np.zeros(w8.shape)
   noise_pseudo_spec = 0

#-------- Data power spectrum -------------
k_out, pk_out, sigma_pk_out, k_nyquist= run_ipoker(data_map, res, mask,scale, beta, remove_1st_bin = remove_1st_bin, delta_l_over_l = delta_l_over_l  ) 


#-------- Error bars ----------------------
# Init signal cu_t
 map_t, cu_t, k_mapp = cls2map( l, clnoise, nx, ny, res,  l_cutoff = k_nyquist)

# Compute the theoretical binned input spectrum corresponding to cu_t
cmn2cb

pk_res = []

if(remove_1st_bin == True): cb_t_in_th = cb_t_in_th[1:]

for imc in range(0, nmc-1):
   #Signal
   cls2map, 

   #Add white noise
   cls2map,
   map_t = map_t + noise

   #Power spectrum
   run_ipoker, 

   pk_res.append( pk )


mc_reduce, pk_res, pk_final, sigma_pk_final, cov_mat, xcorr

#---Plot---
"""
