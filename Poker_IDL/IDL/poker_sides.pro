;; Copyright 2010, Nicolas Ponthieu, Julien Grain, Guilaine Lagache
;; 
;; This file is part of Poker.
;; 
;; Poker is free software: you can redistribute it and/or modify
;; it under the terms of the GNU General Public License as published by
;; the Free Software Foundation, either version 3 of the License, or
;; (at your option) any later version.
;; 
;; Poker is distributed in the hope that it will be useful,
;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;; GNU General Public License for more details.
;; 
;; You should have received a copy of the GNU General Public License
;; along with Poker.  If not, see <http://www.gnu.org/licenses/>.
;; 
;; For more information about Poker, see http://www.ias.u-psud.fr/poker
;; ======================================================================


; /home/nponthieu/Soft/Data

; Mask
maskfile = "/home/mvancuyck/pypoker/mask.fits"
;maskfile = "/home/nponthieu/Soft/Data/mask_apod.fits"

;; CMB

;index     = -3
;beta      = 3 ; binning index
;mapfile   = "/home/nponthieu/Soft/Data/cmb_map.fits"
;noisefile = "/home/nponthieu/Soft/Data/cmb_noise.fits"
;l         = readfits( "/home/nponthieu/Soft/Data/cmb_l.fits")
;clt       = readfits( "/home/nponthieu/Soft/Data/cmb_clt.fits")
;clnoise   = readfits( "/home/nponthieu/Soft/Data/cmb_clnoise.fits")


index     = 0
beta      = 0 ; binning index
mapfile   = "/home/mvancuyck/pypoker/map.fits" ;MJy/sr
;l         = readfits( "/home/nponthieu/Soft/Data/dust_l.fits")

delta_l_over_l = 0.1
log_binning = 1


 ; to improve Mbb's conditioning
remove_1st_bin = 1

; set to 1 to add noise to the simulated data
include_noise  = 0

; Number of random realizations used to determine the noise pseudo-power
; spectrum and the covariance matrix
nmc = 100

; Size of the zero padded large patch that contains the oberved patch
scale = 1.2 ; 2

; Get data
map   = readfits( mapfile) ;MJy/sr

res_arcmin = 0.08333333 ;sxpar( header, "cdelt1")
fwhm = 0.34353685 ;arcmin
sigma = fwhm*!fwhm2sigma*!arcmin2rad ;rad



; keep_avg is not set by default to remove the average of the data from the map
; before power spectrum estimation
keep_avg = 0
nx = n_elements( map[*,0])
ny = n_elements( map[0,*])
nx_large = long( scale*nx)
ny_large = long( scale*ny)

; Init various poker variables
parfile = "./poker.par"
mask= readfits( maskfile, header)
; mask = dblarr(nx,ny) +  1.0d0


give_map_k, res_arcmin*!arcmin2rad, dblarr(nx_large,ny_large), k_map ; k [rad-1]
k_map = k_map*2.0d0*!dpi ;l 
beam= exp( -k_map^2*sigma^2/2.0d0) ; a_lm convention (not Cl) for ipoker, fft(beam)



lmax = !dpi/(res_arcmin*!arcmin2rad) ; Nyquist
cu_t = dblarr(nx_large,ny_large)
wk = where( k_map gt 0 and k_map le lmax, nwk)
cu_t[wk] = k_map[wk]^index ;kmap * beta

ipoker, map, res_arcmin,k_out, pk_out, a_beam = beam, mask=mask, nx_large=nx_large, ny_large=ny_large, beta=beta, remove_1st_bin=remove_1st_bin, delta_l_over_l=delta_l_over_l, out_arrays=arrays, out_params=params


k_out   = arrays.k_out
nbins   = n_elements( k_out)
k_nyquist = !dpi/(res_arcmin*!arcmin2rad)


l = k_out
clt = pk_out

;;--------------------------------------- Noise ------------------------------------------------------
if include_noise eq 1 then begin

   ;; Add simulated noise to the simulated data
   noise = readfits( noisefile)
   map = map + noise

   ;; Compute the noise pseudo-spectrum
   llcb_res     = dblarr( nmc, nbins+1)

   ;; Init cu_t_noise
   cls2map, l, clnoise, nx_large, ny_large, res_arcmin, map_noise, cu_t_noise, k_map, $
            /force, l_cutoff=k_nyquist
   
   for imc=0L, nmc-1 do begin


      cls2map, l, clnoise, nx_large, ny_large, res_arcmin, noise, cu_t_noise, k_map, $
               /no_k_map_reset, /no_cu_t_reset, /force

      ipoker, noise[0:nx-1,0:ny-1], res_arcmin, k_out, dummy, $
              in_arrays=arrays, in_params=params, out_arrays=out_arrays
      
      llcb_res[imc, *] = out_arrays.pseudo_pk

   endfor

   print, "mc reduce noise"
   mc_reduce, llcb_res, noise_pseudo_spec



endif else begin
   clnoise           = 0.0d0
   ;cu_t_noise        = w8*0.0d0
   noise_pseudo_spec = 0.0d0
endelse

noise_pseudo_spec = 0.0d0
arrays.noise_pseudo_spec = noise_pseudo_spec

;;---------------------------------- Data power spectrum ------------------------------------------
ipoker, map, res_arcmin, k_out, pk_out, in_arrays=arrays, in_params=params

;;--------------------------------------- Error bars -----------------------------------------------
; Init signal cu_t
cls2map, l, clt, nx_large, ny_large, res_arcmin, map, cu_t, k_map, $
         /force, /no_k_map_reset, l_cutoff=k_nyquist

; Compute the theoretical binned input spectrum corresponding to cu_t
cmn2cb, arrays.map_b, arrays.map_k_binning*cu_t, cb_t_in_th

pk_res = dblarr( nmc, nbins)

if remove_1st_bin eq 1 then cb_t_in_th = cb_t_in_th[1:*]

for imc=0, nmc-1 do begin
   percent_status, imc, nmc, 5

   ;; Signal
   cls2map, l, clt, nx_large, ny_large, res_arcmin, map_t, cu_t, k_map, $
            l_cutoff = lmax, /no_k_map_reset, /no_cu_t_reset, /force

   ;; Add white noise
   ;cls2map, l, clt, nx_large, ny_large, res_arcmin, noise, cu_t_noise, k_map, /no_k_map_reset, /no_cu_t_reset, /force
   ;map_t = map_t + noise

   ;; Power spectrum
   ipoker, map_t[0:nx-1, 0:ny-1], res_arcmin, k_out, pk, $
           in_arrays=arrays, in_params=params

   pk_res[imc,*] = pk

endfor

mc_reduce, pk_res, pk_final, sigma_pk_final, cov_mat, xcorr

;poker_writefits, "pk_final_idl.fits", pk_final, /silent
;poker_writefits, "sigma_pk_final_idl.fits", sigma_pk_final, /silent

;;------------------------------------------ PLOTS ----------------------------------------------
xra = [100, 5e5]
if index eq !undef then begin
   fl = l*(l+1)/(2*!dpi)
   prefact = k_out*(k_out+1)/(2*!dpi)
   ytitle='k!u2!n P!dk!n/(2!7p!3) [!7l!3K!u2!n]'
   yra = [1e-2, 1e1]
endif else begin
   fl = l*0.0d0 + 1.0d0
   prefact = 1.0d0
   ytitle="P!d!k!3!n"
   yra = [1e-2, 1e1]
endelse


;; Estimator average performance
wind, 1, 1, /free, /large
xlog   = 1
ylog   = 1
xs     = 1
ys     = 1
xpos_0 = 0.15
ypos_0 = 0.25
xpos_1 = 0.97
ypos_1 = 0.95
!p.charsize=1.5
!p.position = [xpos_0, ypos_0, xpos_1, ypos_1]
!x.margin = [10,10]
!p.thick=2
!x.thick=2
!y.thick=2
;plot, xra, yra, ylog=ylog, xlog=xlog, $
;      ytitle=ytitle, xs=xs, ys=ys, /nodata, ycharsize=ycharsize, xcharsize=1e-5
oplot, l, fl*clt, line=2
;oplot, l, fl*(clt + clnoise)
oploterror, k_out, prefact * pk_final, prefact * sigma_pk_final * sqrt(nmc), psym=4;, col=70, errcol=70
;print, 'Power spectra and error bars'
;For ii=0, n_elements(k_out)-1 DO print, k_out(ii), prefact*pk_final(ii), prefact * sigma_pk_final(ii) * sqrt(nmc)  
legendastro, ['Input Data'], col=!p.color, line=0, /right, chars=chars, box=0
xyouts, 3900, 2.2e-7, "POKER", /data
oploterror, [7000], [2.5e-7], [1e-7], psym=4

!y.tickformat = ""
!p.position = [xpos_0, 0.07, xpos_1, ypos_0]
plot, k_out, pk_final/cb_t_in_th, xra=xra, xs=xs, $
      /nodata, yra=[0.5,1.5], /ys, /xlog, xtitle='k', ytitle='b', /noerase, ycharsize=ycharsize
oplot, [1, 2*max(k_out)], [1,1]
oploterror, k_out, k_out^beta*pk_final/cb_t_in_th, k_out^beta*sigma_pk_final/cb_t_in_th, psym=4;, col=70, errcol=70
legendastro, 'Estimator bias', box=0
!p.position = 0

;; Result
wind, 1, 1, /free
plot, xra, yra, ylog=ylog, xlog=xlog, $
      ytitle=ytitle, xs=xs, ys=ys, /nodata, ycharsize=ycharsize, xcharsize=1e-5, title='Data'
oplot, l, fl*clt, line=2
oplot, l, fl*(clt + clnoise)
oploterror, k_out, prefact * pk_out, prefact * sigma_pk_final * sqrt(nmc), psym=4, col=250, errcol=250


poker_writefits, "pk_final_idl.fits",prefact*pk_out, /silent
poker_writefits, "sigma_pk_final_idl.fits", prefact * sigma_pk_final * sqrt(nmc), /silent

;; Bin to bin correlation matrix
wind, 2, 2, /free
!y.tickformat=""
!p.charsize = 3
!x.margin = [5,5]
surface, abs(xcorr), /lego

stop

end
