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




;-------------------------------------------
;pro check_cross

; TG<<GG, ca marche mal
;; index_plot = 0
;; a1 = 1.0d0
;; a2 = 1.0d0
;; ampl_tt = 1.0d0  & index_tt  = -1
;; ampl_gg = 0.1d0  & index_gg  = -1
;; ampl_tg = 0.01d0  & index_TG = -1

; TG>>GG, ca marche bien
index_plot = 1
a1 = 1.0d0
a2 = 1.0d0
ampl_tt = 1.0d0  & index_tt  = -1
ampl_gg = 0.01d0  & index_gg  = -1
ampl_tg = 0.05d0  & index_TG = -1

;; ;;Meme pour des spectres tres differents
;; index_plot = 2
;; a1 = 1.0d0
;; a2 = 1.0d0
;; ampl_tt = 10d0 & index_tt  = -1
;; ampl_gg = 1000000d0  & index_gg  = -6
;; ampl_tg = 10 & index_TG = -3

;; ;; pentes inversees ?
;; index_plot = 3
;; a1 = 1.0d0 ; switch on/off the correlated part in map_g
;; a2 = 1.0d0 ; 0.0d0 ; 1.0d0 ; switch on/off the uncorrelated part in map_g
;; ampl_tt = 1.d0  & index_tt  = -1   ;; CIB
;; ampl_gg = 1.d-3 & index_gg  = 0.75 ;; Unnormalized Lensing
;; ampl_tg = 1.d-2 & index_TG = -2.2  ;; Unnormalized

index_plot = strtrim( long(index_plot),2)+'.'+strtrim(long(a1),2)+strtrim(long(a2),2)
print, index_plot

;;; Pentes inversees mais avec TG>>GG
;index_plot = 3
;a1 = 1.0d0
;a2 = 1.0d0
;ampl_tt = 1.d0  & index_tt  = -1  ;; CIB
;ampl_gg = 1.d-11 & index_gg  =  0.75 ;; Unnormalized Lensing
;ampl_tg = 100.0d0 & index_TG  = -2.2  ;; Unnormalized 

png = 0
nmc = 100
nx = 100 ; 128 ; 512
ny = nx
res_arcmin = 1.0d0
scale = 1.5
apod_length=20
nx_large = long( scale *  nx)
ny_large = long( scale *  ny)
npix = long(nx_large)*long(ny_large)

kmax = 2*!dpi*sqrt(2)/2.0d0/(res_arcmin*!arcmin2rad)
kk = dindgen( 10000)/9999.*( kmax-1) + 1

; Beams
fwhm_t = 0.0d0 ;5.0d0 ; arcmin
fwhm_g = 0.0d0
bypass  = 0 ; 

beta_tt = -index_tt
beta_gg = -index_gg
beta_tg = -index_tg

sigma_t = fwhm_t*!fwhm2sigma*!arcmin2rad
sigma_g = fwhm_g*!fwhm2sigma*!arcmin2rad

ctt = ampl_tt * kk^index_tt
cgg = ampl_gg * kk^index_gg
ctg = ampl_tg * kk^index_tg

n_tt = 0.0d0 ;0.01 * ampl_tt * 100.0d0^index_tt ; 0.0d0 ;
n_gg = 0.0d0 ;0.01 * ampl_gg * 100.0d0^index_gg ; 0.0d0 ;

col_t     = 250
col_g     = 70
col_tg    = 150
col_ratio = 100

xra = minmax(kk) ; [100, 1e4]
yra = [min([ctg^2/ctt, ctt,cgg,ctg])/10., 2]
wind, 1, 1, /free, /large
plot_oo, kk, ctt, yra=yra, /ys, xra=xra
oplot, [100,100], yra, line=1
oplot,   kk, ctt, col=col_t
oplot,   kk, cgg, col=col_g
oplot,   kk, ctg, col=col_tg
oplot,   kk, ctg^2/ctt, col=col_ratio
oplot, kk, ctt + kk*0.0d0 + n_tt, col=col_t, line=2
oplot, kk, cgg + kk*0.0d0 + n_gg, col=col_g, line=2
legendastro, ['TT', 'GG', 'TG', 'TG!u2!n/TT'], col=[col_t, col_g, col_tg, col_ratio], line=0

; to apodize boundaries
mask = dblarr(nx,ny) + 1.0d0
nholes = 0
radius_ratio = 30
poker_make_mask, nx, ny, scale, nholes, radius_ratio, mask_large, $
                 apod_length=apod_length, patch=patch, nx_large=nx_large, ny_large=ny_large
wind, 1, 1, /f
db, mask_large
wp = where( patch eq 1)
;mask[wp] = mask_large[wp]


; beams
give_map_k, res_arcmin*!arcmin2rad, dblarr(nx_large,ny_large), k_map
k_map = k_map*2.0d0*!dpi

beam_t = exp( -k_map^2*sigma_t^2/2.0d0) ; a_lm convention (not Cl) for ipoker
beam_g = exp( -k_map^2*sigma_g^2/2.0d0)
rem = 1

; If you choose the same beta for each auto and cross spectrum, you can launch
; ipoker only once, otherwise do as follows:

; TG
ipoker, dblarr(nx,ny), res_arcmin, k, pk, $
        map1=dblarr(nx,ny), a_beam=beam_t, b_beam=beam_g, beta=beta_tg, $
        mask=mask, rem=rem, nx_large=nx_large, ny_large=ny_large, $
        out_params=params_tg, out_arrays=arrays_tg, bypass=bypass


one_large = dblarr(nx_large,ny_large) + 1.0d0
; TT
ipoker, dblarr(nx,ny), res_arcmin, k, pk, $
        a_beam=beam_t, beta=beta_tt, $
        mask=mask, rem=rem, nx_large=nx_large, ny_large=ny_large, $
        out_params=params_tt, out_arrays=arrays_tt, bypass=bypass

ipoker, dblarr(nx,ny), res_arcmin, k, pk, $
        beta=beta_tt, $
        mask=mask, rem=rem, nx_large=nx_large, ny_large=ny_large, $
        out_params=pars, out_arrays=arrs, bypass=bypass


; GG
ipoker, dblarr(nx,ny), res_arcmin, k, pk, $
        a_beam=beam_g, beta=beta_gg, $
        mask=mask, rem=rem, nx_large=nx_large, ny_large=ny_large, $
        out_params=params_gg, out_arrays=arrays_gg, bypass=bypass



k_map  = arrays_tt.k_map
lmap_rad_x = double(nx_large*res_arcmin*!arcmin2rad)
lmap_rad_y = double(ny_large*res_arcmin*!arcmin2rad)
norm = sqrt( (nx_large/lmap_rad_x)*(ny_large/lmap_rad_y))
ic = dcomplex( 0.00d0, 1.00d0)

; Generate gaussian amplitudes
cu_tt  = k_map*0.0d0
cu_gg  = k_map*0.0d0
cu_tg = k_map*0.0d0

lmax = !dpi/(res_arcmin*!arcmin2rad) ; max( k_mapx) ; Nyquist
wk = where( k_map gt 0 and k_map le lmax, nwk)

cu_tt[wk] = ampl_tt * k_map[wk]^index_tt
cu_gg[wk] = ampl_gg * k_map[wk]^index_gg
cu_tg[wk] = ampl_tg * k_map[wk]^index_tg


;; Noise pseudo-spectra
if n_tt ne 0 or n_gg ne 0 then begin
   nk = n_elements(k)
   pk_tt_res = dblarr( nmc, nk+1)
   pk_gg_res = dblarr( nmc, nk+1)
   for imc=0L, nmc-1 do begin
      percent_status, imc, nmc, 10

      noise    = dblarr( nx_large, ny_large)
      n2       = long(nx_large)*long(ny_large)
      noise[*] = randomn( seed, n2, /double)
                                ;noise    = noise - avg(noise)
                                ;noise    = noise/stddev(noise)
      dmn1     = fft( noise, /double)
      dmn1     = dmn1 * norm

      noise[*] = randomn( seed, n2, /double)
                                ;noise    = noise - avg(noise)
                                ;noise    = noise/stddev(noise)
      dmn2     = fft( noise, /double)
      dmn2     = dmn2 * norm

      amn_t = dmn1*sqrt( n_tt)
      amn_g = dmn2*sqrt( n_gg)

      map_t = double( fft( amn_t, /inverse, /double))

      map_g = double( fft( amn_g, /inverse, /double))
      map_t = map_t[0:nx-1, 0:ny-1]
      map_g = map_g[0:nx-1, 0:ny-1]

      ipoker, map_t, res_arcmin, k, pk_tt, in_params=params_tt, in_arrays=arrays_tt, out_arrays=out_arrays_tt
      ipoker, map_g, res_arcmin, k, pk_gg, in_params=params_gg, in_arrays=arrays_gg, out_arrays=out_arrays_gg

      pk_tt_res[ imc,*] = out_arrays_tt.pseudo_pk
      pk_gg_res[ imc,*] = out_arrays_gg.pseudo_pk
   endfor
   mc_reduce, pk_tt_res, noise_pseudo_spec_tt
   mc_reduce, pk_gg_res, noise_pseudo_spec_gg

   

endif else begin
   noise_pseudo_spec_tt = k*0.0d0
   noise_pseudo_spec_gg = k*0.0d0
endelse
   
;; Update input arrays
arrays_tt.noise_pseudo_spec = noise_pseudo_spec_tt
arrays_gg.noise_pseudo_spec = noise_pseudo_spec_gg

;; Full error bars
nk = n_elements(k)
pk_tg_res = dblarr( nmc, nk)
pk_tt_res = dblarr( nmc, nk)
pk_gg_res = dblarr( nmc, nk)
mean_list_t = dblarr( nmc)
mean_list_g = dblarr( nmc)

n2       = long(nx_large)*long(ny_large)
noise    = dblarr( nx_large, ny_large)
;noise[*] = randomn( seed, n2, /double)
;dmn2     = fft( noise, /double)
;dmn2     = dmn2*norm
;index_plot = index_plot + 'a'


for imc=0L, nmc-1 do begin
   percent_status, imc, nmc, 10

   noise[*] = randomn( seed, n2, /double)
   name = "idl_arr/"+strtrim(imc,1)+"_A.fits"
   poker_writefits, name , noise, /silent
   ;noise    = noise - avg(noise)
   ;noise    = noise/stddev(noise)
   dmn1     = fft( noise, /double)
   dmn1     = dmn1*norm


   name = "idl_arr/"+strtrim(imc,1)+"_real_dmn1.fits"
   poker_writefits, name , real_part(dmn1), /silent
   name = "idl_arr/"+strtrim(imc,1)+"_imag_dmn1.fits"
   poker_writefits, name , imaginary(dmn1), /silent


   noise[*] = randomn( seed, n2, /double)
   name = "idl_arr/"+strtrim(imc,1)+"_B.fits"
   poker_writefits, name , noise, /silent
   ;noise    = noise - avg(noise)
   ;noise    = noise/stddev(noise)
   dmn2     = fft( noise, /double)
   dmn2     = dmn2*norm


   name = "idl_arr/"+strtrim(imc,1)+"_real_dmn2.fits"
   poker_writefits, name , real_part(dmn2), /silent
   name = "idl_arr/"+strtrim(imc,1)+"_imag_dmn2.fits"
   poker_writefits, name , imaginary(dmn2), /silent


   noise[*] = randomn( seed, n2, /double)
   ;noise    = noise - avg(noise)
   ;noise    = noise/stddev(noise)
   dmn3     = fft( noise, /double)
   dmn3     = dmn3 * norm



   noise[*] = randomn( seed, n2, /double)
   ;noise    = noise - avg(noise)
   ;noise    = noise/stddev(noise)
   dmn4     = fft( noise, /double)
   dmn4     = dmn4 * norm




   ;; amn amplitudes
   amn_t     = dmn1*0.0d0
   amn_g     = dmn2*0.0d0
   amn_n_t   = dmn3*sqrt( n_tt)
   amn_n_g   = dmn4*sqrt( n_gg)

   amn_t[wk] = dmn1[wk] * sqrt( cu_tt[wk])


   amn_g[wk] = a1*cu_tg[wk]/cu_tt[wk]*amn_t[wk] + a2*sqrt( cu_gg[wk] - cu_tg[wk]*cu_tg[wk]/cu_tt[wk])*dmn2[wk]
   
   name = "idl_arr/"+strtrim(imc,1)+"_real_amng.fits"
   poker_writefits, name , real_part(amn_g), /silent


   name = "idl_arr/"+strtrim(imc,1)+"_imag_amng.fits"
   poker_writefits, name , imaginary(amn_g), /silent


   
   ;poker_writefits, "imag_part_dmn1.fits",  imaginary(dmn1), /silent
   ;poker_writefits, "real_part_dmn1.fits",  real_part(dmn1), /silent

   ;poker_writefits, "imag_part_dmn2.fits",  imaginary(dmn2), /silent
   ;poker_writefits, "real_part_dmn2.fits",  real_part(dmn2), /silent

   ;poker_writefits, "imag_part_amnt.fits",  imaginary(amn_t), /silent
   ;poker_writefits, "real_part_amnt.fits",  real_part(amn_t), /silent

   ;poker_writefits, "imag_part_amng.fits",  imaginary(amn_g), /silent
   ;poker_writefits, "real_part_amng.fits",  real_part(amn_g), /silent

   ;; Include beam
   amn_t = amn_t * beam_t
   amn_g = amn_g * beam_g

   ;; Add noise
   amn_t = amn_t + amn_n_t
   amn_g = amn_g + amn_n_g

   ;; Turn into maps
   map_t = double( fft( amn_t, /inverse, /double))
   map_g = double( fft( amn_g, /inverse, /double))


   name = "idl_arr/"+strtrim(imc,1)+"_mapt.fits"
   poker_writefits, name , map_t, /silent

   name = "idl_arr/"+strtrim(imc,1)+"_mapg.fits"
   poker_writefits, name , map_g, /silent
   
   ;print, "after save maps"
   ;STOP


   map_t = map_t[0:nx-1, 0:ny-1]
   map_g = map_g[0:nx-1, 0:ny-1]

   mean_t =  mean(map_t)
   mean_list_t[imc] = mean_t
   meanless_t = map_t - mean_t
   mean_g = mean(map_g)
   mean_list_g[imc] = mean_g
   meanless_s = map_g - mean_g
   poker_writefits, "idl_arr/"+strtrim(imc,1)+"mapt_meanless.fits" ,meanless_t, /silent
   poker_writefits, "idl_arr/"+strtrim(imc,1)+"mapg_meanless.fits" ,meanless_s, /silent

   ;; Power spectra
   cmn_name = "idl_arr/"+strtrim(imc,1)+"_cmn_TG.fits"
   ipoker, map_t, res_arcmin, k, pk_tg, map1=map_g, in_params=params_tg, in_arrays=arrays_tg, out_arrays=out_arrays_tg
   cmn_name = "idl_arr/"+strtrim(imc,1)+"_cmn_T.fits"
   ipoker, map_t, res_arcmin, k, pk_tt,             in_params=params_tt, in_arrays=arrays_tt, out_arrays=out_arrays_tt
   cmn_name = "idl_arr/"+strtrim(imc,1)+"_cmn_G.fits"
   ipoker, map_g,  res_arcmin, k, pk_gg,             in_params=params_gg, in_arrays=arrays_gg, out_arrays=out_arrays_gg

   pk_tt_res[ imc,*] = pk_tt
   pk_gg_res[ imc,*] = pk_gg
   pk_tg_res[ imc,*] = pk_tg


   name = "idl_arr/"+strtrim(imc,1)+"_pkgg.fits"
   poker_writefits, name , pk_gg, /silent


   name = "idl_arr/"+strtrim(imc,1)+"_pktt.fits"
   poker_writefits, name , pk_tt, /silent

   name = "idl_arr/"+strtrim(imc,1)+"_pktg.fits"
   poker_writefits, name , pk_tg, /silent
   
endfor

poker_writefits, "idl_arr/mean_list_t.fits" , mean_list_t, /silent
poker_writefits, "idl_arr/mean_list_g.fits" , mean_list_g, /silent

mc_reduce, pk_tg_res, pk_tg_avg, sigma_pk_tg
mc_reduce, pk_tt_res, pk_tt_avg, sigma_pk_tt
mc_reduce, pk_gg_res, pk_gg_avg, sigma_pk_gg

print, noise_mean.mean(), noise_std.mean()
;print, mean_conj_1.mean(), mean_conj_2.mean(), mean_conj_3.mean(), mean_conj_4.mean()
;print, real_part_mean_1.mean(), real_part_mean_2.mean(), real_part_mean_3.mean(), real_part_mean_4.mean()
;print, imag_part_mean_1.mean(), imag_part_mean_2.mean(), imag_part_mean_3.mean(), imag_part_mean_4.mean()
;print, real_part_std_1.mean(), real_part_std_2.mean(), real_part_std_3.mean(), real_part_std_4.mean()
;print, imag_part_std_1.mean(), imag_part_std_2.mean(), imag_part_std_3.mean(), imag_part_std_4.mean()


;poker_writefits, "m_bb_m1_tt_scale12.fits", arrays_tt.x_mtt_bb_m1, /silent
print, 'Power spectra and error bars'
For ii=0, n_elements(k)-1 DO print, k(ii), pk_tt_avg(ii), pk_gg_avg(ii), pk_tg_avg(ii)
STOP

yra = [ min( abs( [pk_tt_avg[1:*], pk_gg_avg[1:*], pk_tg_avg[1:*]]))/10., 2]
yra = [min([ctg^2/ctt, ctt,cgg,ctg])/10., 2]
!p.multi=0
wind, 1, 1, /free, /large
outplot, file='Output_spectra_'+strtrim(index_plot,2), png=png
plot_oo, k_map[wk], cu_tt[wk], yra=yra, /ys, /xs, /nodata
oplot, kk, ctt + kk*0.0d0 + n_tt, col=col_t, line=2
oplot, kk, cgg + kk*0.0d0 + n_gg, col=col_g, line=2
oplot, k_map[wk], cu_tt[wk]*beam_t[wk]^2, line=2, col=col_t
oplot, k_map[wk], cu_tg[wk]^2/cu_tt[wk]
oplot, k_map[wk], cu_tt[wk], col=col_t
oplot, kk, ctg, col=col_tg
oploterror, k, pk_tt_avg, sigma_pk_tt, psym=4, thick=2, col=col_t, errcol=col_t
oplot, k_map[wk], cu_gg[wk], col=col_g
oplot, k_map[wk], cu_gg[wk]*beam_g[wk]^2, col=col_g, line=2
oploterror, k, pk_gg_avg, sigma_pk_gg, psym=4, thick=2, col=col_g, errcol=col_g
oplot, [1e-10, 1e10], [1,1], line=2
oplot, k_map[wk], cu_tg[wk]/sqrt(cu_tt[wk]*cu_gg[wk]), col=col_ratio
oplot, k, pk_tg_avg/sqrt( pk_tt_avg*pk_gg_avg), psym=4, thick=2, col=col_ratio
oplot, k_map[wk], cu_tg[wk], col=col_tg
oplot, k_map[wk], cu_tg[wk]*beam_t[wk]*beam_g[wk], line=2, col=col_tg
oploterror, k, pk_tg_avg, sigma_pk_tg, psym=4, thick=2, col=col_tg, errcol=col_tg
oplot, k_map[wk], cu_tg[wk]^2/cu_tt[wk], line=2
legendastro, ['TT', 'GG', 'TG', 'TG/sqrt(TT*GG)', 'TG!u2!n/TT'], col=[col_t, col_g, col_tg, col_ratio, !p.color], line=0
legendastro, ['C!uTT!n = '+string( ampl_tt,format="(F5.2)")+' k!u'+String(index_tt,format="(F5.2)")+'!n', $
         'C!uGG!n = '+string( ampl_gg,format="(F5.2)")+' k!u'+String(index_gg,format="(F5.2)")+'!n', $
         'C!uTG!n = '+string( ampl_tg,format="(F5.2)")+' k!u'+String(index_tg,format="(F5.2)")+'!n'], /bottom, chars=1.5
outplot, /close, /verb
;png, "plot_"+strtrim(index_plot,2)+".png"


cos_var_tt = k*0.0d0
cos_var_gg = k*0.0d0
cos_var_tg = k*0.0d0
bintab = arrays_tt.bintab

dk = (bintab - shift( bintab, 1))[1:*]
w = where( k gt 0)
fsky = nx*ny*(res_arcmin*!arcmin2rad)^2/(4*!dpi)
;;fsky = total(mask)*(res_arcmin*!arcmin2rad)^2/(4*!dpi)
fsky = nx_large*ny_large*(res_arcmin*!arcmin2rad)^2/(4*!dpi)
cos_var_tt[w] = (sqrt(2./(fsky*(2*k+1)*dk)) * pk_tt_avg)[w]
cos_var_gg[w] = (sqrt(2./(fsky*(2*k+1)*dk)) * pk_gg_avg)[w]
cos_var_tg[w] = (sqrt(1./(fsky*(2*k+1)*dk)) * sqrt( pk_tg_avg^2 + pk_tt_avg*pk_gg_avg))[w]

;; check intrinsic sigma_b ?
ipoker, map_t, res_arcmin, k, pk_tt,             sigma_poker_tt, in_params=params_tt, in_arrays=arrays_tt, out_arrays=out_arrays_tt
ipoker, map_g, res_arcmin, k, pk_gg,             sigma_poker_gg, in_params=params_gg, in_arrays=arrays_gg, out_arrays=out_arrays_gg
ipoker, map_t, res_arcmin, k, pk_tg, map1=map_g, sigma_poker_tg, in_params=params_tg, in_arrays=arrays_tg, out_arrays=out_arrays_tg

wind, 2, 2, /free, /large
!p.multi=[0,2,2]
yra = [min(cu_tt[wk])/10, max(cu_tt[wk])*10]
plot_oo, k_map[wk], cu_tt[wk], /xs, title='TT', yra=yra, /ys
oploterror, k, pk_tt_avg, sigma_pk_tt*sqrt(nmc), psym=4, col=col_t, errcol=col_t
oploterror, k*1.05, pk_tt_avg, cos_var_tt, psym=1, col=70, errcol=70
oploterror, k*0.95, pk_tt_avg, sigma_poker_tt, psym=1, col=150, errcol=150
legendastro, ['MC', 'Analytic', 'sigma_b'], col=[col_t, 70,150], line=0, chars=1.5

yra = [min(cu_gg[wk])/10, max(cu_gg[wk])*10]
plot_oo, k_map[wk], cu_gg[wk], /xs, title='GG', yra=yra, /ys
oploterror, k, pk_gg_avg, sigma_pk_gg*sqrt(nmc), psym=4, col=col_t, errcol=col_t
oploterror, k*1.05, pk_gg_avg, cos_var_gg, psym=1, col=70, errcol=70
oploterror, k*0.95, pk_gg_avg, sigma_poker_gg, psym=1, col=150, errcol=150
legendastro, ['MC', 'Analytic', 'sigma_b'], col=[col_t, 70,150], line=0, chars=1.5

yra = [min(cu_tg[wk])/10, max(cu_tg[wk])*10]
plot_oo, k_map[wk], cu_tg[wk], /xs, title='TG', yra=yra, /ys
oploterror, k, pk_tg_avg, sigma_pk_tg*sqrt(nmc), psym=4, col=col_t, errcol=col_t
oploterror, k*1.05, pk_tg_avg, cos_var_tg, psym=1, col=70, errcol=70
oploterror, k*0.95, pk_tg_avg, sigma_poker_tg, psym=1, col=150, errcol=150
legendastro, ['MC', 'Analytic', 'sigma_b'], col=[col_t, 70,150], line=0, chars=1.5
!p.multi=0

end
