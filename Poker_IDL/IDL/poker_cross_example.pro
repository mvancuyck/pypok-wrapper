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





nx = 100 ; 128
ny = 100 ; nx
res_arcmin = 1.0d0
nx_large = 128 ; nx
ny_large = 128 ; ny

;; Number of Monte Carlo simulations
nmc = 500 ; 100 ; 500

;;;; Input power spectra with different power laws
ampl_tt  = 10d0       & index_tt  = -1
ampl_gg  = 1000000d0  & index_gg  = -6
ampl_tg = 10          & index_TG = -3

;;;; Input power spectra with different power laws
;;ampl_tt = 10d0     & index_tt = -1
;;ampl_gg = 1000d0   & index_gg = -3
;;ampl_tg = 10000d0  & index_TG = -3

;; Input power spectra with different power laws
ampl_tt = 1e6      & index_tt = -3
ampl_gg = 1000d0   & index_gg = -3
ampl_tg = 10000d0  & index_TG = -3

; Beams
fwhm_t  = 5.0d0 ; arcmin
fwhm_g  = 3.0d0
bypass  = 0

; Define binning index
beta_tt = -index_tt
beta_gg = -index_gg
beta_tg = -index_tg

; Signal
kmax = 2*!dpi*sqrt(2)/2.0d0/(res_arcmin*!arcmin2rad)
kk   = dindgen( 10000)/9999.*( kmax-1) + 1
ctt = ampl_tt * kk^index_tt
cgg = ampl_gg * kk^index_gg
ctg = ampl_tg * kk^index_tg

; noise
n_tt = ampl_tt * 1000.0d0^index_tt
n_gg = ampl_gg * 1000.0d0^index_gg

col_t     = 250
col_g     = 70
col_tg    = 150

xra = [100, 1e4]
yra = [min([ctg^2/ctt, ctt,cgg,ctg])/10., 2]
wind, 1, 1, /free, /large
plot_oo, kk, ctt, yra=yra, /ys, xra=xra
oplot,   kk, ctt, col=col_t
oplot,   kk, cgg, col=col_g
oplot,   kk, ctg, col=col_tg
oplot, kk, ctt + kk*0.0d0 + n_tt, col=col_t, line=2
oplot, kk, cgg + kk*0.0d0 + n_gg, col=col_g, line=2
oplot, kk, ctg^2/ctt, col=100
legend, ['TT', 'GG', 'TG'], col=[col_t, col_g, col_tg], line=0

; to apodize boundaries
mask = dblarr(nx,ny) + 1.0d0
nholes = 0
radius_ratio = 30
apod_length = 30
poker_make_mask, nx, ny, scale, nholes, radius_ratio, mask_large, $
                 apod_length=apod_length, patch=patch, nx_large=nx_large, ny_large=ny_large

wind, 1, 1, /free
db, mask_large, title='mask_large'
wp = where( patch eq 1)
mask = dblarr(nx,ny)
mask[*] = mask_large[wp]

;wind, 2, 2, /free
;db, mask, title='mask'

; Format beams for ipoker
give_map_k, res_arcmin*!arcmin2rad, dblarr(nx_large,ny_large), k_map
k_map = k_map*2.0d0*!dpi
sigma_t = fwhm_t*!fwhm2sigma*!arcmin2rad
sigma_g = fwhm_g*!fwhm2sigma*!arcmin2rad
beam_t = exp( -k_map^2*sigma_t^2/2.0d0) ; a_lm convention (not Cl) for ipoker
beam_g = exp( -k_map^2*sigma_g^2/2.0d0)
remove_1st_bin = 1
delta_l_over_l = 0.1

; If you choose the same beta for each auto and cross spectrum, you can launch
; ipoker only once, otherwise do as follows:

; TT
ipoker, dblarr(nx,ny), res_arcmin, k, pk, a_beam=beam_t, beta=beta_tt, $
        mask=mask, remove_1st_bin=remove_1st_bin, nx_large=nx_large, ny_large=ny_large, $
        out_params=params_tt, out_arrays=arrays_tt, delta_l_over_l=delta_l_over_l

stop



; GG
ipoker, dblarr(nx,ny), res_arcmin, k, pk, a_beam=beam_g, beta=beta_gg, $
        mask=mask, remove_1st_bin=remove_1st_bin, nx_large=nx_large, ny_large=ny_large, $
        out_params=params_gg, out_arrays=arrays_gg, delta_l_over_l=delta_l_over_l

; GG
ipoker, dblarr(nx,ny), res_arcmin, k, pk, $
        map1=dblarr(nx,ny), a_beam=beam_t, b_beam=beam_g, beta=beta_tg, $
        mask=mask, remove_1st_bin=remove_1st_bin, nx_large=nx_large, ny_large=ny_large, $
        out_params=params_tg, out_arrays=arrays_tg, delta_l_over_l=delta_l_over_l

k_map      = arrays_tt.k_map
lmap_rad_x = double(nx_large*res_arcmin*!arcmin2rad)
lmap_rad_y = double(ny_large*res_arcmin*!arcmin2rad)
norm       = sqrt( (nx_large/lmap_rad_x)*(ny_large/lmap_rad_y))
ic         = dcomplex( 0.00d0, 1.00d0)

; Generate gaussian amplitudes
cu_tt = k_map*0.0d0
cu_gg = k_map*0.0d0
cu_tg = k_map*0.0d0

lmax = !dpi/(res_arcmin*!arcmin2rad) ; Nyquist
wk = where( k_map gt 0 and k_map le lmax, nwk)

cu_tt[wk] = ampl_tt * k_map[wk]^index_tt
cu_gg[wk] = ampl_gg * k_map[wk]^index_gg
cu_tg[wk] = ampl_tg * k_map[wk]^index_tg

;; Noise pseudo-spectra
nk = n_elements(k)
pk_tt_res = dblarr( nmc, nk+1)
pk_gg_res = dblarr( nmc, nk+1)
for imc=0L, nmc-1 do begin
   percent_status, imc, nmc, 10

   noise    = dblarr( nx_large, ny_large)
   n2       = long(nx_large)*long(ny_large)
   noise[*] = randomn( seed, n2, /double)
   noise    = noise - avg(noise)
   noise    = noise/stddev(noise)
   dmn1     = fft( noise, /double)

   noise[*] = randomn( seed, n2, /double)
   noise    = noise - avg(noise)
   noise    = noise/stddev(noise)
   dmn2     = fft( noise, /double)

   amn_t = dmn1*sqrt( n_tt)
   amn_g = dmn2*sqrt( n_gg)

   map_t = double( fft( amn_t, /inverse, /double))
   map_g = double( fft( amn_g, /inverse, /double))

   ;; Restrict to the lower left corner to simulate non periodic boundary conditions
   map_t = map_t[0:nx-1, 0:ny-1]
   map_g = map_g[0:nx-1, 0:ny-1]

   ipoker, map_t, res_arcmin, k, pk_tt, in_params=params_tt, in_arrays=arrays_tt, out_arrays=out_arrays_tt
   ipoker, map_g, res_arcmin, k, pk_gg, in_params=params_gg, in_arrays=arrays_gg, out_arrays=out_arrays_gg

   pk_tt_res[ imc,*] = out_arrays_tt.pseudo_pk
   pk_gg_res[ imc,*] = out_arrays_gg.pseudo_pk
endfor
mc_reduce, pk_tt_res, noise_pseudo_spec_tt
mc_reduce, pk_gg_res, noise_pseudo_spec_gg

;; Update input arrays with noise pseudo-spectra
arrays_tt.noise_pseudo_spec = noise_pseudo_spec_tt
arrays_gg.noise_pseudo_spec = noise_pseudo_spec_gg

;; Theoretical binned spectra
cmn2cb, arrays_tt.map_b, arrays_tt.map_k_binning * cu_tt, cb_tt_in_th
cmn2cb, arrays_tg.map_b, arrays_tg.map_k_binning * cu_tg, cb_tg_in_th
cmn2cb, arrays_gg.map_b, arrays_gg.map_k_binning * cu_gg, cb_gg_in_th

if remove_1st_bin eq 1 then begin
   cb_tt_in_th = cb_tt_in_th[1:*]
   cb_tg_in_th = cb_tg_in_th[1:*]
   cb_gg_in_th = cb_gg_in_th[1:*]
endif

;; Full error bars
nk = n_elements(k)
pk_tg_res = dblarr( nmc, nk)
pk_tt_res = dblarr( nmc, nk)
pk_gg_res = dblarr( nmc, nk)
pk_x_res  = dblarr( nmc, nk)

for imc=0L, nmc-1 do begin
   percent_status, imc, nmc, 10

   noise    = dblarr( nx_large, ny_large)
   n2       = long(nx_large)*long(ny_large)
   noise[*] = randomn( seed, n2, /double)
   noise    = noise - avg(noise)
   noise    = noise/stddev(noise)
   dmn1     = fft( noise, /double)

   noise[*] = randomn( seed, n2, /double)
   noise    = noise - avg(noise)
   noise    = noise/stddev(noise)
   dmn2     = fft( noise, /double)

   noise[*] = randomn( seed, n2, /double)
   noise    = noise - avg(noise)
   noise    = noise/stddev(noise)
   dmn3     = fft( noise, /double)

   noise[*] = randomn( seed, n2, /double)
   noise    = noise - avg(noise)
   noise    = noise/stddev(noise)
   dmn4     = fft( noise, /double)

   dmn1 = dmn1*norm
   dmn2 = dmn2*norm

; Init amn fields
; nx and ny appear in the definition of norm because the amplitudes in Fourier space
; are generated by the fft of a white noise map in real space
   amn_t   = dmn1*0.0d0
   amn_g   = dmn2*0.0d0
   amn_n_t = dmn3*sqrt( n_tt)
   amn_n_g = dmn4*sqrt( n_gg)

   ;; Fill amn
   amn_t = dmn1 * sqrt( cu_tt)
   amn_g[wk] = cu_tg[wk]/cu_tt[wk]*amn_t[wk] + sqrt( cu_gg[wk] - cu_tg[wk]*cu_tg[wk]/cu_tt[wk])*dmn2[wk]

   ;; Include beam
   amn_t = amn_t * beam_t
   amn_g = amn_g * beam_g

   ;; Add noise
   amn_t = amn_t + amn_n_t
   amn_g = amn_g + amn_n_g

   ;; Output map
   map_t = double( fft( amn_t, /inverse, /double))
   map_g = double( fft( amn_g, /inverse, /double))
   map_t = map_t[0:nx-1, 0:ny-1]
   map_g = map_g[0:nx-1, 0:ny-1]

   ipoker, map_t, res_arcmin, k, pk_tt,             in_params=params_tt, in_arrays=arrays_tt, out_arrays=out_arrays_tt
   ipoker, map_g, res_arcmin, k, pk_gg,             in_params=params_gg, in_arrays=arrays_gg, out_arrays=out_arrays_gg
   ipoker, map_t, res_arcmin, k, pk_tg, map1=map_g, in_params=params_tg, in_arrays=arrays_tg, out_arrays=out_arrays_tg

   pk_tt_res[ imc,*] = pk_tt
   pk_gg_res[ imc,*] = pk_gg
   pk_tg_res[ imc,*] = pk_tg
   
endfor
mc_reduce, pk_tg_res, pk_tg_avg, sigma_pk_tg, cov_mat_tg, xcorr_tg
mc_reduce, pk_tt_res, pk_tt_avg, sigma_pk_tt, cov_mat_tt, xcorr_tt
mc_reduce, pk_gg_res, pk_gg_avg, sigma_pk_gg, cov_mat_gg, xcorr_gg


yra=[ min( abs( [pk_tt_avg[1:*], pk_gg_avg[1:*], pk_tg_avg[1:*]]))/10., 2]
!p.multi=0
wind, 1, 1, /free, /large
plot_oo, k_map[wk], cu_tt[wk], yra=yra, /ys, /xs, /nodata
oplot, kk, ctt + kk*0.0d0 + n_tt, col=col_t, line=2
oplot, kk, cgg + kk*0.0d0 + n_gg, col=col_g, line=2
oplot,   k_map[wk], cu_tt[wk]*beam_t[wk]^2, line=2, col=col_t
oplot, k_map[wk], cu_tg[wk]^2/cu_tt[wk]
oplot, k_map[wk], cu_tt[wk], col=col_t
oploterror, k, pk_tt_avg, sigma_pk_tt, psym=4, col=col_t, errcol=col_t
oplot, k_map[wk], cu_gg[wk], col=col_g
oplot, k_map[wk], cu_gg[wk]*beam_g[wk]^2, col=col_g, line=2
oploterror, k, pk_gg_avg, sigma_pk_gg, psym=4, col=col_g, errcol=col_g
;oplot, [1e-10, 1e10], [1,1], line=2
;oplot, k_map[wk], cu_tg[wk]/sqrt(cu_tt[wk]*cu_gg[wk]), col=col_ratio
;oplot, k, pk_tg_avg/sqrt( pk_tt_avg*pk_gg_avg), psym=4, col=col_ratio
oplot, k_map[wk], cu_tg[wk], col=col_tg
oplot, k_map[wk], cu_tg[wk]*beam_t[wk]*beam_g[wk], line=2, col=col_tg
oploterror, k, pk_tg_avg, sigma_pk_tg, psym=4, col=col_tg, errcol=col_tg

oplot, k_map[wk], cu_tg[wk]^2/cu_tt[wk], line=2
;;legend, ['TT', 'GG', 'TG', 'TG/sqrt(TT*GG)', 'TG!u2!n/TT'], col=[col_t, col_g, col_tg, col_ratio, !p.color], line=0
;;legend, ['C!uTT!n = '+string( ampl_tt,format="(F5.2)")+' k!u'+String(index_tt,format="(F5.2)")+'!n', $
;;         'C!uGG!n = '+string( ampl_gg,format="(F5.2)")+' k!u'+String(index_gg,format="(F5.2)")+'!n', $
;;         'C!uTG!n = '+string( ampl_tg,format="(F5.2)")+' k!u'+String(index_tg,format="(F5.2)")+'!n'], /bottom, chars=1.5

legend, ['TT', 'GG', 'TG'], col=[col_t, col_g, col_tg], line=0
legend, ['C!uTT!n = '+string( ampl_tt,format="(F5.2)")+' k!u'+String(index_tt,format="(F5.2)")+'!n', $
         'C!uGG!n = '+string( ampl_gg,format="(F5.2)")+' k!u'+String(index_gg,format="(F5.2)")+'!n', $
         'C!uTG!n = '+string( ampl_tg,format="(F5.2)")+' k!u'+String(index_tg,format="(F5.2)")+'!n'], /bottom, chars=1.5


cos_var_tt = k*0.0d0
cos_var_gg = k*0.0d0
cos_var_tg = k*0.0d0
bintab = arrays_tt.bintab

;; Approximate cosmic variance
dk            = (bintab - shift( bintab, 1))[1:*]
w             = where( k gt 0)
;;fsky          = nx_large*ny_large*(res_arcmin*!arcmin2rad)^2/(4*!dpi)
fsky          = total(mask)*(res_arcmin*!arcmin2rad)^2/(4*!dpi)
cos_var_tt[w] = (sqrt(2./(fsky*(2*k+1)*dk)) * pk_tt_avg)[w]
cos_var_gg[w] = (sqrt(2./(fsky*(2*k+1)*dk)) * pk_gg_avg)[w]
cos_var_tg[w] = (sqrt(1./(fsky*(2*k+1)*dk)) * sqrt( pk_tg_avg^2 + pk_tt_avg*pk_gg_avg))[w]


;; Results
wind, 2, 2, /free, /large
!p.multi=[0,2,2]
yra = [min(cu_tt[wk])/10, max(cu_tt[wk])*10]
plot_oo, k_map[wk], cu_tt[wk], /xs, title='TT', yra=yra, /ys
oploterror, k, pk_tt_avg, sigma_pk_tt*sqrt(nmc), psym=4, col=col_t, errcol=col_t
oploterror, k*1.05, pk_tt_avg, cos_var_tt, psym=1, col=70, errcol=70
legend, ['MC', 'Analytic'], col=[col_t, 70], line=0, chars=1.5

yra = [min(cu_gg[wk])/10, max(cu_gg[wk])*10]
plot_oo, k_map[wk], cu_gg[wk], /xs, title='GG', yra=yra, /ys
oploterror, k, pk_gg_avg, sigma_pk_gg*sqrt(nmc), psym=4, col=col_t, errcol=col_t
oploterror, k*1.05, pk_gg_avg, cos_var_gg, psym=1, col=70, errcol=70
legend, ['MC', 'Analytic'], col=[col_t, 70], line=0, chars=1.5

yra = [min(cu_tg[wk])/10, max(cu_tg[wk])*10]
plot_oo, k_map[wk], cu_tg[wk], /xs, title='TG', yra=yra, /ys
oploterror, k, pk_tg_avg, sigma_pk_tg*sqrt(nmc), psym=4, col=col_t, errcol=col_t
oploterror, k*1.05, pk_tg_avg, cos_var_tg, psym=1, col=70, errcol=70
legend, ['MC', 'Analytic'], col=[col_t, 70], line=0, chars=1.5
!p.multi=0

;; Check for bias w.r.t. the binned input theoretical power spectra
wind, 1, 1, /free, /large
!p.multi=[0,1,3]
plot, k, k^beta_tt*pk_tt_avg/cb_tt_in_th, xra=xra, xs=xs, title='TT', $
      /nodata, yra=[0.5,1.5], /ys, /xlog, xtitle='k', ytitle='b', ycharsize=ycharsize
oplot, [1, 2*max(k)], [1,1]
oploterror, k, k^beta_tt*pk_tt_avg/cb_tt_in_th, k^beta_tt*sigma_pk_tt/cb_tt_in_th, psym=4
legend, 'Estimator bias', box=0

plot, k, k^beta_tg*pk_tg_avg/cb_tg_in_th, xra=xra, xs=xs, title='TG', $
      /nodata, yra=[0.5,1.5], /ys, /xlog, xtitle='k', ytitle='b', ycharsize=ycharsize
oplot, [1, 2*max(k)], [1,1]
oploterror, k, k^beta_tg*pk_tg_avg/cb_tg_in_th, k^beta_tg*sigma_pk_tg/cb_tg_in_th, psym=4
legend, 'Estimator bias', box=0

plot, k, k^beta_gg*pk_gg_avg/cb_gg_in_th, xra=xra, xs=xs, title='GG', $
      /nodata, yra=[0.5,1.5], /ys, /xlog, xtitle='k', ytitle='b', ycharsize=ycharsize
oplot, [1, 2*max(k)], [1,1]
oploterror, k, k^beta_gg*pk_gg_avg/cb_gg_in_th, k^beta_gg*sigma_pk_gg/cb_gg_in_th, psym=4
legend, 'Estimator bias', box=0
!p.multi=0

;; Look at cross-correlations between bins
wind, 2, 2, /free, /large
!p.multi=[0,2,2]
surface, abs(xcorr_tt), /lego, title='TT'
surface, abs(xcorr_tg), /lego, title='TG'
surface, abs(xcorr_gg), /lego, title='GG'
!p.multi=0

end
