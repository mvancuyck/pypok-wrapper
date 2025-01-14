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

pro cu2maps, cu_tt, res_arcmin, map_t, $
             cu_gg=cu_gg, cu_tg=cu_tg, map_g=map_g, $
             beam_t=beam_t, beam_g=beam_g, threshold=threshold, $
             above_nyquist=above_nyquist, seed_tt=seed_tt, seed_gg=seed_gg

if n_params() lt 1 then begin
   message, /info, "Calling sequence: "
   print, "cu2maps, cu_tt, res_arcmin, map_t, $"
   print, "         cu_gg=cu_gg, cu_tg=cu_tg, map_g=map_g, $"
   print, "         beam_t=beam_t, beam_g=beam_g, seed_tt=seed_tt, seed_gg=seed_gg, threshold=threshold"
   print, "         above_nyquist=above_nyquist"
   return
endif

if not keyword_set(threshold) then threshold = 1e-5

nx = n_elements( cu_tt[*,0])
ny = n_elements( cu_tt[0,*])

lmap_rad_x = double(nx*res_arcmin*!arcmin2rad)
lmap_rad_y = double(ny*res_arcmin*!arcmin2rad)
norm       = sqrt( (nx/lmap_rad_x)*(ny/lmap_rad_y))

nyquist_mask = 1.0d0 ; init
if not keyword_set(above_nyquist) then begin
   k_nyquist = !dpi/(res_arcmin*!arcmin2rad)
   give_map_k, res_arcmin, cu_tt, map_k
   wa = where( 2*!dpi*map_k/!arcmin2rad gt k_nyquist, nwa)
   nyquist_mask     = map_k*0.0d0 + 1.0d0
   nyquist_mask[wa] = 0.0d0
endif

n2 = long(nx)*long(ny)
noise    = dblarr( nx, ny)
if keyword_set(seed_tt) then seed = seed_tt
noise[*] = randomn( seed, n2, /double)

dmn1 = fft( noise, /double)
dmn1 = dmn1*norm

if not keyword_set(beam_t) then beam_t = dblarr(nx,ny) + 1.0d0
if not keyword_set(beam_g) then beam_g = dblarr(nx,ny) + 1.0d0

amn_t = dmn1 * sqrt( cu_tt) * nyquist_mask
map_t = double( fft( amn_t*beam_t, /inverse, /double))

wk = lindgen( n2-1) + 1 ; only the first elements is k=0
if keyword_set(cu_gg) then begin

   if not keyword_set(cu_tg) then cu_tg = cu_tt*0.0d0

   sqrt_arg = cu_gg[wk] - cu_tg[wk]*cu_tg[wk]/cu_tt[wk]

   ww = where( sqrt_arg lt 0.d0, nww)
   if nww ne 0 then sqrt_arg[ww] = 0.0d0

   if keyword_set(seed_gg) then seed1 = seed_gg
   noise[*] = randomn( seed1, n2, /double)
   dmn2     = fft( noise, /double)
   dmn2     = dmn2*norm

   amn_g     = dmn2*0.0d0
   amn_g[wk] = cu_tg[wk]/cu_tt[wk]*amn_t[wk] + sqrt( sqrt_arg)*dmn2[wk]

   amn_g = amn_g * beam_g * nyquist_mask
   map_g = double( fft( amn_g, /inverse, /double))
endif

end
