
nx = 128
ny = nx
res_arcmin = 1.0d0
nx_large = nx
ny_large = ny
npix = long(nx_large)*long(ny_large)


nmc = 1000
ipoker, fltarr(nx_large,ny_large), res_arcmin, k, pk_tt, /bypass

nk = n_elements(k)
pk_tt_res = dblarr( nmc, nk)
pk_gg_res = dblarr( nmc, nk)
pk_tg_res  = dblarr( nmc, nk)
lmap_rad_x = double(nx_large*res_arcmin*!arcmin2rad)
lmap_rad_y = double(ny_large*res_arcmin*!arcmin2rad)
norm = sqrt( (nx_large/lmap_rad_x)*(ny_large/lmap_rad_y))

noise    = dblarr( nx_large, ny_large)
n2       = long(nx_large)*long(ny_large)
noise[*] = randomn( seed, n2, /double)
dmn1     = fft( noise, /double)
dmn1     = dmn1*norm

for imc=0L, nmc-1 do begin
   percent_status, imc, nmc, 10

   ;noise    = dblarr( nx_large, ny_large)
   ;n2       = long(nx_large)*long(ny_large)
   ;noise[*] = randomn( seed, n2, /double)
   ;dmn1     = fft( noise, /double)

   noise[*] = randomn( seed, n2, /double)
   dmn2     = fft( noise, /double)

   ;dmn1 = dmn1*norm
   dmn2 = dmn2*norm

   ;; Fill amn
   amn_t = dmn1
   amn_g = dmn2

   ;; Output map
   map_t = double( fft( amn_t, /inverse, /double))
   map_g = double( fft( amn_g, /inverse, /double))

   ipoker, map_t, res_arcmin, k, pk_tt, /bypass
   ipoker, map_g, res_arcmin, k, pk_gg, /bypass
   ipoker, map_t, res_arcmin, k, pk_tg, map1=map_g, /bypass

   pk_tt_res[ imc,*] = pk_tt
   pk_gg_res[ imc,*] = pk_gg
   pk_tg_res[ imc,*] = pk_tg
   
endfor
mc_reduce, pk_tg_res, pk_tg_avg, sigma_pk_tg
mc_reduce, pk_tt_res, pk_tt_avg, sigma_pk_tt
mc_reduce, pk_gg_res, pk_gg_avg, sigma_pk_gg

xra = [100, 2e4]
wind, 1, 1, /free, /large
!p.multi=[0,2,2]
plot_oo, k, pk_tt_avg, title='TT (1)', xra=xra, /xs
oploterror, k, pk_tt_avg, sigma_pk_tt, psym=8
legend, 'N!dMC!n = '+strtrim(nmc,2)
plot_oo, k, pk_gg_avg, title='GG (1)', xra=xra, /xs
oploterror, k, pk_gg_avg, sigma_pk_gg, psym=8
plot_oo, k, pk_tg_avg, title='TG (1)', xra=xra, /xs
oploterror, k, pk_tg_avg, sigma_pk_tg, psym=8
ploterror, k, pk_tg_avg, sigma_pk_tg, title='TG (1)', xra=xra, /xs
!p.multi=0
stop

pk_tt_res = dblarr( nmc, nk)
pk_gg_res = dblarr( nmc, nk)
pk_tg_res  = dblarr( nmc, nk)
for imc=0L, nmc-1 do begin
   percent_status, imc, nmc, 10

   dmn1 = dcomplexarr( nx_large, ny_large)
   dmn2 = dcomplexarr( nx_large, ny_large)
   rr = randomn( seed, 2*n_elements(map_t)) ;/sqrt(2.)
   ii = randomn( seed, 2*n_elements(map_t)) ;/sqrt(2.)
   dmn1[*] = dcomplex( rr[0:npix-1], ii[0:npix-1])
   dmn2[*] = dcomplex( rr[npix:*], ii[npix:*])
   norm = sqrt( 1.0d0/lmap_rad_x/lmap_rad_y)
   dmn1 = dmn1 * norm
   dmn2 = dmn2 * norm

   ;; Fill amn
   amn_t = dmn1
   amn_g = dmn2

   ;; Output map
   map_t = double( fft( amn_t, /inverse, /double))
   map_g = double( fft( amn_g, /inverse, /double))

   ipoker, map_t, res_arcmin, k, pk_tt, /bypass
   ipoker, map_g, res_arcmin, k, pk_gg, /bypass
   ipoker, map_t, res_arcmin, k, pk_tg, map1=map_g, /bypass

   pk_tt_res[ imc,*] = pk_tt
   pk_gg_res[ imc,*] = pk_gg
   pk_tg_res[ imc,*] = pk_tg
   
endfor
mc_reduce, pk_tg_res, pk_tg_avg, sigma_pk_tg
mc_reduce, pk_tt_res, pk_tt_avg, sigma_pk_tt
mc_reduce, pk_gg_res, pk_gg_avg, sigma_pk_gg

xra = [100, 2e4]
wind, 1, 1, /free, /large
!p.multi=[0,2,2]
plot_oo, k, pk_tt_avg, title='TT (2)', xra=xra, /xs
oploterror, k, pk_tt_avg, sigma_pk_tt, psym=8
legend, 'N!dMC!n = '+strtrim(nmc,2)
plot_oo, k, pk_gg_avg, title='GG (2)', xra=xra, /xs
oploterror, k, pk_gg_avg, sigma_pk_gg, psym=8
plot_oo, k, pk_tg_avg, title='TG (2)', xra=xra, /xs
oploterror, k, pk_tg_avg, sigma_pk_tg, psym=8
!p.multi=0

;; All real parts in one shot and imaginary parts in one shot
nmc = 100 ; force to a smaller number
rr = randomn( seed, 2*npix*nmc)
ii = randomn( seed, 2*npix*nmc)
pk_tt_res = dblarr( nmc, nk)
pk_gg_res = dblarr( nmc, nk)
pk_tg_res  = dblarr( nmc, nk)
for imc=0L, nmc-1 do begin
   percent_status, imc, nmc, 10

   dmn1 = dcomplexarr( nx_large, ny_large)
   dmn2 = dcomplexarr( nx_large, ny_large)
   dmn1[*] = dcomplex( rr[ 2*imc*npix     : 2*imc*npix      + npix-1], ii[2*imc*npix     : 2*imc*npix      + npix-1])
   dmn2[*] = dcomplex( rr[ 2*imc*npix+npix: 2*imc*npix+npix + npix-1], ii[2*imc*npix+npix: 2*imc*npix+npix + npix-1])
   norm = sqrt( 1.0d0/lmap_rad_x/lmap_rad_y)
   dmn1 = dmn1 * norm
   dmn2 = dmn2 * norm

   ;; Fill amn
   amn_t = dmn1
   amn_g = dmn2

   ;; Output map
   map_t = double( fft( amn_t, /inverse, /double))
   map_g = double( fft( amn_g, /inverse, /double))

   ipoker, map_t, res_arcmin, k, pk_tt, /bypass
   ipoker, map_g, res_arcmin, k, pk_gg, /bypass
   ipoker, map_t, res_arcmin, k, pk_tg, map1=map_g, /bypass

   pk_tt_res[ imc,*] = pk_tt
   pk_gg_res[ imc,*] = pk_gg
   pk_tg_res[ imc,*] = pk_tg
   
endfor
mc_reduce, pk_tg_res, pk_tg_avg, sigma_pk_tg
mc_reduce, pk_tt_res, pk_tt_avg, sigma_pk_tt
mc_reduce, pk_gg_res, pk_gg_avg, sigma_pk_gg

xra = [100, 2e4]
wind, 1, 1, /free, /large
!p.multi=[0,2,2]
plot_oo, k, pk_tt_avg, title='TT (3)', xra=xra, /xs
oploterror, k, pk_tt_avg, sigma_pk_tt, psym=8
legend, 'N!dMC!n = '+strtrim(nmc,2)
plot_oo, k, pk_gg_avg, title='GG (3)', xra=xra, /xs
oploterror, k, pk_gg_avg, sigma_pk_gg, psym=8
plot_oo, k, pk_tg_avg, title='TG (3)', xra=xra, /xs
oploterror, k, pk_tg_avg, sigma_pk_tg, psym=8
!p.multi=0

;; All in one shot
nmc = 100 ; force to a smaller number
pk_tt_res = dblarr( nmc, nk)
pk_gg_res = dblarr( nmc, nk)
pk_tg_res  = dblarr( nmc, nk)

rr = randomn( seed, 4*npix*nmc)
res_rad = res_arcmin*!arcmin2rad
ipoker, map_t, res_arcmin, k, pk_tt, /bypass, out_arrays = out_arrays
for imc=0L, nmc-1 do begin
   percent_status, imc, nmc, 10

   dmn1 = dcomplexarr( nx_large, ny_large)
   dmn2 = dcomplexarr( nx_large, ny_large)
   i1 = 2*imc*npix
   i2 = i1 + npix
   i3 = i2 + npix
   i4 = i3 + npix
   dmn1[*] = dcomplex( rr[ i1 : i1 + npix-1], rr[ i2 : i2 + npix-1])
   dmn2[*] = dcomplex( rr[ i3 : i3 + npix-1], rr[ i3 : i3 + npix-1])
   norm = sqrt( 1.0d0/lmap_rad_x/lmap_rad_y)
   dmn1 = dmn1 * norm
   dmn2 = dmn2 * norm

   ;; Fill amn
   amn_t = dmn1
   amn_g = dmn2

   ;; Output map
   map_t = double( fft( amn_t, /inverse, /double))
   map_g = double( fft( amn_g, /inverse, /double))

   ;ipoker, map_t, res_arcmin, k, pk_tt, /bypass
   ;ipoker, map_g, res_arcmin, k, pk_gg, /bypass
   ;ipoker, map_t, res_arcmin, k, pk_tg, map1=map_g, /bypass

   p2 = double( abs( dmn1)^2) * npix * res_rad^2
   cmn2cb, out_arrays.map_b, p2, pk_tt

   p2 = double( abs( dmn2)^2) * npix * res_rad^2
   cmn2cb, out_arrays.map_b, p2, pk_gg
   
   p2 = double( 0.5d0*( dmn1*conj(dmn2) + conj(dmn1)*dmn2)) * npix * res_rad^2
   cmn2cb, out_arrays.map_b, p2, pk_tg

   pk_tt_res[ imc,*] = pk_tt
   pk_gg_res[ imc,*] = pk_gg
   pk_tg_res[ imc,*] = pk_tg
   
endfor
mc_reduce, pk_tg_res, pk_tg_avg, sigma_pk_tg
mc_reduce, pk_tt_res, pk_tt_avg, sigma_pk_tt
mc_reduce, pk_gg_res, pk_gg_avg, sigma_pk_gg

xra = [100, 2e4]
wind, 1, 1, /free, /large
!p.multi=[0,2,2]
plot_oo, k, pk_tt_avg, title='TT (4)', xra=xra, /xs
oploterror, k, pk_tt_avg, sigma_pk_tt, psym=8
legend, 'N!dMC!n = '+strtrim(nmc,2)
plot_oo, k, pk_gg_avg, title='GG (4)', xra=xra, /xs
oploterror, k, pk_gg_avg, sigma_pk_gg, psym=8
plot_oo, k, pk_tg_avg, title='TG (4)', xra=xra, /xs
oploterror, k, pk_tg_avg, sigma_pk_tg, psym=8
!p.multi=0


end
