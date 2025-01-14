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



pro my_p2, map, npix, res_rad, p2, map1=map1

;+
; NAME:
;   my_p2
;
; PURPOSE:
;   Computes the 2D power spectrum of a 2D image or the cross-spectrum of 2
;   images of same size
;
; CALLING SEQUENCE:
;  my_p2, map, npix, res_rad, p2
;
; INPUTS:
;     map : image
;    npix : number of elements of image
;  res_rad: map resolution in radians
;
; OPTIONAL INPUTS:
;  none
;
; KEYWORD PARAMETERS:
;  none
;
; OUTPUTS:
;    p2 : 2D power spectrum of map
;
; OPTIONAL OUTPUTS:
;   none
;
; COMMON BLOCKS:
;   none
;
; SIDE EFFECTS:
;   none
;
; RESTRICTIONS:
;   none
;
; PROCEDURE:
;
; MODIFICATION HISTORY:
;    Dec. 2010, N. Ponthieu
;-

if n_params() lt 1 then begin
   message, /info, "Calling sequence: "
   print, "my_p2, map, npix, res_rad, p2, map1=map1"
   return
endif


if keyword_set(map1) then begin
   ft  = fft( map, /double)
   ft1 = fft( map1, /double)
   p2  = 0.5d0*( ft*conj(ft1) + conj(ft)*ft1) * npix * res_rad^2
endif else begin
   p2 = abs( fft( map, /double))^2 * npix * res_rad^2
endelse

p2 = double(p2)

end
