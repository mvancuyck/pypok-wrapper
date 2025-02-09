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


pro ipoker_arrays, map, mask, res_arcmin, a_beam, b_beam, pars, out_arrays, out_params, $
                   bypass_mtt_bb=bypass_mtt_bb, bintab_in=bintab_in, clean_poker=clean_poker

if n_params() lt 1 then begin
   dl_unix, 'ipoker_arrays'
   return
endif

if not keyword_set(bypass_mtt_bb) then bypass_mtt_bb = 0

;; Map parameters
nx         = n_elements( map[*,0])
ny         = n_elements( map[0,*])

; Create the large mask
poker_make_mask, nx, ny, 1, 0, 0, w8, $
                 patch=patch, nx_large=pars.nx_large, ny_large=pars.ny_large, apod_length= 20
npix = long(pars.nx_large)*long(pars.ny_large)
wp = where( patch eq 1)
w8[wp] = mask * w8[wp] ; to apply apodization if requested


; Binning parameters
give_map_k, pars.res_pix, dblarr(pars.nx_large,pars.ny_large), k_map, k_mapx
k_map  = k_map  * 2.0d0*!dpi
wk     = where( k_map ne 0)

; Bins
kmin  = min( k_map[wk]) ;; on the "large" map
;;if not keyword_set(dk_min) then dk_min = 2*kmin
if pars.dk_min lt 0 then pars.dk_min = 2.0d0 * 2.0d0*!dpi/(nx*pars.res_pix) ;;2*min(k of the data)

if keyword_set(bintab_in) then begin
   bintab = bintab_in
endif else begin
   k_range   = [kmin, pars.k_nyquist]
   make_bintab, k_range, pars.dk_min, bintab, /float, delta_l_over_l=pars.delta_l_over_l, log=pars.log_binning
;; Add bin 0 and take margin with round off errors for the 1st bin
   bintab[0] = bintab[0]*0.99
   bintab    = [0, bintab]
endelse

; Init header
header    = strarr(9)
header[0] = "SIMPLE  =                    T / Written by IDL/ipoker.pro                      "
header[1] = "BITPIX  =                  -64 / IEEE double precision floating point           "
header[2] = "NAXIS   =                    2 /       /                                        "
header[3] = "NAXIS1  =                  1   /Number of positions along axis 1                "    ; place holder
header[4] = "NAXIS2  =                  1   /Number of positions along axis 2                "    ; place holder
header[5] = "BLOCKED =                    T         /                                        "
header[6] = "CDELT1  =        "+string( res_arcmin, "(F13.11)")+" /Resolution arcmin                               "
header[7] = "CDELT2  =        "+string( res_arcmin, "(F13.11)")+" /Resolution arcmin                               "
header[8] = "END                                                                             "

;; Write mask to disk
poker_writefits, pars.mask, w8, header, /silent

;; If no mask nor beam is provided, save computation time                                               
if max( abs(w8-1.0d0)) eq 0 and not keyword_set(a_beam) and not keyword_set(b_beam) then bypass_mtt_bb = 1

;; Write patch to disk
;poker_writefits, pars.patch, patch, header, /silent

;; Beam transfer function (no longer complex)
poker_writefits, pars.beam,  a_beam, /silent
poker_writefits, pars.beam1, b_beam, /silent

; Derive binning laws
map_k_binning = k_map*0.0d0 + 1.0d0
map_k_binning[wk] = k_map[wk]^pars.beta
poker_writefits, pars.map_k_binning, map_k_binning, /silent

;; Pass bintab to F90 (reduce units)
poker_writefits, pars.input_bintab, bintab/(2.0d0*!dpi/pars.res_pix), /silent
nbins = n_elements( bintab)-1

; Get mode mixing matrix
igive_map_b, pars.nx_large, pars.ny_large, bintab/(2.0d0*!dpi/pars.res_pix), map_b_junk, xbintab


if bypass_mtt_bb eq 1 then begin
   x_mtt_bb = identity( nbins)
   map_b    = map_b_junk
endif else begin
   
   parfile = "poker.par"
   params2ascii, pars, parfile
   
   if (!arch eq "franklin") or (!arch eq "parallel") then begin
      print, ""
      ;print, "=================="
      ;print, "ready to run poker_count_task and qsub run_poker.txt."
      ;print, "Then press .c to continue."
      ;print, "=================="
      print, " PARRALLELE "
      spawn, "mpirun poker_mbb_mpi"+parfile
   endif else begin
      print, "PAS PARRALELE"
      spawn, "poker_mbb "+parfile
      print, "" & print, "" & print, "MBB COMPUTED."
   endelse

   map_b    = readfits( pars.file_map_b,/silent)
   x_mtt_bb = readfits( pars.file_mtt_bb_x,/silent)

endelse

;; Compute bin addresses
cmn2cb, map_b, k_map, k_out
;cmn2cb_new, bintab, map_b, k_map, k_out_new
;stop
poker_writefits, pars.k_out, k_out, /silent


; Discard DC bin to improve Mtt_bb's conditioning
if pars.remove_1st_bin then begin
   x_mtt_bb_m1 = invert( x_mtt_bb[1:*,1:*])
   k_out       = k_out[1:*]
endif else begin
   x_mtt_bb_m1 = invert( x_mtt_bb)
endelse

out_arrays = {mask:mask, w8:w8, patch:patch, wp:wp, k_map:k_map, $
              bintab:bintab, a_beam:a_beam, b_beam:b_beam, $
              map_k_binning:map_k_binning, $
              x_mtt_bb:x_mtt_bb, $
              x_mtt_bb_m1:x_mtt_bb_m1, map_b:map_b, $
              k_out:k_out, pk_out:k_out*0.0d0, pseudo_pk:dblarr(nbins), noise_pseudo_spec:dblarr(nbins), $
              xbintab:xbintab}

if keyword_set(clean_poker) then poker_clean, out_params

end
