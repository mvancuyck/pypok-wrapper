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


pro ipoker_mpi, map, res_arcmin, k_out, pk_out, sigma_pk_out, $
            pb_out=pb_out, mask=mask, map1=map1, $
            a_beam=a_beam, b_beam=b_beam, $
            in_params=in_params, out_params=out_params, $
            in_arrays=in_arrays, out_arrays=out_arrays, $
            bypass_mtt_bb=bypass_mtt_bb, $
            remove_1st_bin=remove_1st_bin, log_binning=log_binning, $
            dir=dir, delta_l_over_l=delta_l_over_l, $
            nx_large=nx_large, ny_large=ny_large, apod_length=apod_length, $
            dk_min=dk_min, beta=beta, bintab_in=bintab_in, $
            noise_pseudo_spec=noise_pseudo_spec, $
            include_pix_window_function=include_pix_window_function, clean_poker=clean_poker, imc=imc
            

if n_params() lt 1 then begin
   dl_unix, 'ipoker'
   return
endif

if not keyword_set(bypass_mtt_bb) then bypass_mtt_bb = 0

;; Map parameters
nx = n_elements( map[*,0])
ny = n_elements( map[0,*])

if keyword_set(in_params) then begin
   out_params = in_params
endif else begin
   ;; Init params structure
   poker_types, out_params
endelse

;;If keywords are set, they override the default parameters or those in in_params
if keyword_set(beta)                        then out_params.beta           = beta
if keyword_set(delta_l_over_l)              then out_params.delta_l_over_l = delta_l_over_l
if keyword_set(include_pix_window_function) then out_params.include_pix_window_function = 1
if keyword_set(dir)                         then out_params.dir = dir
if keyword_set(remove_1st_bin)              then out_params.remove_1st_bin = remove_1st_bin
if keyword_set(apod_length)                 then out_params.apod_length = apod_length
if keyword_set(dk_min)                      then out_params.dk_min = dk_min
if keyword_set(nx_large)                    then out_params.nx_large = nx_large else out_params.nx_large=nx
if keyword_set(ny_large)                    then out_params.ny_large = ny_large else out_params.ny_large=ny

out_params.nx      = nx
out_params.ny      = ny
out_params.res_pix = res_arcmin*!arcmin2rad

out_params.k_nyquist = !dpi/out_params.res_pix

; Instrumental beam to include in mtt_bb
if keyword_set(a_beam) then begin
   s = size( a_beam)
   if (s[1] ne out_params.nx_large) or (s[2] ne out_params.ny_large) then begin
      message, /info, "Incompatible sizes : a_beam and mask_large"
      stop
   endif
endif else begin
   a_beam = dblarr( out_params.nx_large, out_params.ny_large) + 1.0d0
endelse

if keyword_set(b_beam) then begin
   s = size( b_beam)
   if (s[1] ne nx_large) or (s[2] ne ny_large) then begin
      message, /info, "Incompatible sizes : b_beam and mask_large"
      stop
   endif
endif else begin
   b_beam = a_beam
endelse



if keyword_set(in_arrays) then begin
   out_arrays = in_arrays
endif else begin
   if not keyword_set(mask) then mask = dblarr(nx,ny) + 1.0d0
   ipoker_arrays_mpi, map, mask, res_arcmin, a_beam, b_beam, out_params, out_arrays, out_params, $
                  bypass_mtt_bb=bypass_mtt_bb, bintab_in=bintab_in, clean_poker=clean_poker
endelse

; Embed data in the large patch
;; if out_params.keep_avg eq 0 then data_map = map - mean(map) else data_map = map
if out_params.keep_avg eq 0 then data_map = map - mean(map[where( out_arrays.mask gt 0.d0)]) else data_map = map
map_large     = out_arrays.w8*0.0d0
map_large[out_arrays.wp] = data_map



if keyword_set(map1) then begin
   map1_large = out_arrays.w8*0.0d0
   if not keyword_set(keep_avg) then data_map1 = map1 - mean(map1) else data_map1 = map1
   map1_large[out_arrays.wp] = data_map1
endif else begin
   map1_large = map_large
endelse


; Compute the data(+noise) pseudo-spectrum
npix = n_elements( out_arrays.w8)
my_p2, out_arrays.w8*map_large, npix, out_params.res_pix, c_mn_t, map1=map1_large*out_arrays.w8

;poker_writefits, "idl_arr/cmn.fits" , c_mn_t, /silent
cmn2cb, out_arrays.map_b, out_arrays.map_k_binning*c_mn_t, pb_mes, sigma_pb_mes
out_arrays.pseudo_pk = pb_mes


if not keyword_set(noise_pseudo_spec) then noise_pseudo_spec = out_arrays.noise_pseudo_spec



; Discard DC bin to improve Mtt_bb's conditioning
if out_params.remove_1st_bin then begin
   pb_mes            = pb_mes[1:*]
   noise_pseudo_spec = noise_pseudo_spec[1:*]
   sigma_pb_mes      = sigma_pb_mes[1:*]
endif

;; Binned k^beta*P(k) output power spectrum
nk           = n_elements( out_arrays.k_out)
pb_out       = reform( out_arrays.x_mtt_bb_m1##(pb_mes-noise_pseudo_spec), nk)


;dum_1        = reform( out_arrays.x_mtt_bb_m1##(pb_mes+sigma_pb_mes-noise_pseudo_spec), nk)
;dum_2        = reform( out_arrays.x_mtt_bb_m1##(pb_mes-sigma_pb_mes-noise_pseudo_spec), nk)
;sigma_pk_out = (dum_1-dum_2)/2.0d0
sigma_pk_out  = reform( out_arrays.x_mtt_bb_m1##sigma_pb_mes)

;; Effective power spectrum at bin average assuming the process is described by a
;; power law of index=beta


k_out        = out_arrays.k_out
pk_out       = pb_out * k_out^(-out_params.beta)
sigma_pk_out = sigma_pk_out * k_out^(-out_params.beta)

out_arrays.pk_out = pk_out



end
