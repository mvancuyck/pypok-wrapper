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


;+
; INPUTS:
;l   multipoles
;clt angular power spectrum (in units of microK^2 x sr to get a map in microK)
;nx  number of pixels in the x direction
;ny  number of pixels in the y direction
;res_arcmin physical resolution of the map in arcminutes
;
; OUTPUTs
;map_t: the map
;cu_t : interpolated values of the input power spectrum at each k value
;k_map: Fourier k modes values (can be input as well)
;k_mapx: Fourier k modes along x (can be input as well)
;k_mapy: Fourier k modes along y (can be input as well)
;
;KEYWORDS:
;zero_index: set to 1 t simulate an input power spectrum l^0
;index: set to any non zero value to simulate an input spectrum l^index
;l_cutoff: maximum multipole used for the simulation
;no_k_map_reset: set to 1 if k_map is provided in input to avoid to recompute it
;(saves time for Monte-Carlos on large maps)
;no_cu_t_reset: set to 1 with /no_k_map_reset for the same reasons
;ampl: extra amplitude factor in front of the power spectrum
;
;AUTHOR
; Nicolas Ponthieu, IAS
;-

pro cls2map, l, clt, nx, ny, res_arcmin, map_t, cu_t, k_map, k_mapx, k_mapy, $
             zero_index = zero_index, index=index, l_cutoff = l_cutoff, $
             no_k_map_reset=no_k_map_reset, no_cu_t_reset=no_cu_t_reset, ampl=ampl, $
             seed=seed, force=force, fwhm_arcmin=fwhm_arcmin


if n_params() lt 1 then begin
   message, /info, "Calling sequence: "
   print, "cls2map, l, clt, nx, ny, res_arcmin, map_t, cu_t, k_map, k_mapx, k_mapy, $"
   print, "         zero_index = zero_index, index=index, l_cutoff = l_cutoff, $"
   print, "         no_k_map_reset=no_k_map_reset, no_cu_t_reset=no_cu_t_reset, ampl=ampl, $"
   print, "         seed=seed, force=force, fwhm_arcmin=fwhm_arcmin"
   return
endif

nx = long(nx)
ny = long(ny)

if not keyword_set(index) then index = !undef

lmap_rad_x = double(nx*res_arcmin*!arcmin2rad)
lmap_rad_y = double(ny*res_arcmin*!arcmin2rad)


; Generate gaussian amplitudes
noise    = dblarr( nx, ny)
n2       = long(nx)*long(ny)
seed = !NULL
noise[*] = randomn( seed, n2, /double)



;noise    = noise - avg(noise)
;noise    = noise/stddev(noise)
dmn1     = fft( noise, /double)

noise[*] = randomn( seed, n2, /double)
;noise    = noise - avg(noise)
;noise    = noise/stddev(noise)
dmn2     = fft( noise, /double)

; Init amn fields
; nx and ny appear in the definition of norm because the amplitudes in Fourier space
; are generated by the fft of a white noise map in real space
norm = sqrt( (nx/lmap_rad_x)*(ny/lmap_rad_y))

;Define k modes
if not keyword_set( no_k_map_reset) then begin
   give_map_k, res_arcmin*!arcmin2rad, dblarr(nx,ny), k_map, k_mapx, k_mapy

  
   k_map  = k_map  * 2.0d0*!dpi
   k_mapx = k_mapx * 2.0d0*!dpi
   k_mapy = k_mapy * 2.0d0*!dpi
endif

if keyword_set(l_cutoff) then begin
   lmax = l_cutoff
endif else begin
   if (index ne !undef) or keyword_set(zero_index) then begin
      lmax = max( k_map)
   endif else begin
      lmax = max(l)
   endelse
endelse

if keyword_set(fwhm_arcmin) then begin
   sigma = fwhm_arcmin*!arcmin2rad*!fwhm2sigma
   dl    = double(l)
   bl    = exp(-dl*(dl+1)*sigma^2)
endif else begin
   bl = 1.0d0
endelse


;Interpolate input power spectrum
if not keyword_set( no_cu_t_reset) then begin
   cu_t = k_map*0.0d0

   w = where( k_map gt 0 and k_map le lmax, nw)
   if nw eq 0 then begin
      message, /info, "wrong k range"
      stop
   endif else begin

      ;;Power law spectrum
      if not keyword_set(ampl) then ampl = 1.0d0

      if (index ne !undef) or keyword_set(zero_index) then begin

         if keyword_set(fwhm_arcmin) then begin
            sigma = fwhm_arcmin*!arcmin2rad*!fwhm2sigma
            bk    = exp(-k_map*(k_map+1)*sigma^2)
         endif else begin
            bk = 1.0d0
         endelse

         if keyword_set(zero_index) then index = 0

         cu_t[w] = ampl * k_map[w]^index * bk[w]

      endif else begin
         cu_t[w] = ampl * interpol( clt*bl, l, k_map[w])
         w1 = where( cu_t lt 0, nw1)
         if nw1 ne 0 then begin
            if keyword_set(force) then begin
               cu_t[w1] = 0.0d0
            endif else begin
               print, "cu_t < 0 for "+strtrim( nw1,2)+" values"
               stop
            endelse
         endif

      endelse
   endelse
endif
   
;Fill amn_t
amn_t = dmn1 * norm * sqrt( cu_t)

;Output map
map_t = double( fft( amn_t, /inverse, /double))


end
