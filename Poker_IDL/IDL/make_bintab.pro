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



pro make_bintab, l, delta_l_min, bintab, delta_l_over_l=delta_l_over_l, $
                 log_bin=log_bin, delta_l_max = delta_l_max, $
                 bin_width = bin_width, float=float

;+
; NAME:
;   make_bintab
;
; PURPOSE:
;   Generates bin boundaries and central values from an input array of multipoles
;
; CALLING SEQUENCE:
;  make_bintab, l, delta_l_min, bintab, delta_l_over_l=delta_l_over_l, $
;                 log_bin=log_bin, delta_l_max = delta_l_max, $
;                 bin_width = bin_width, float=float
;
; INPUTS:
;     l : array (2 elements min) or range of multipoles to be binned
;  delta_l_min : minimum width of output bins
;
; OPTIONAL INPUTS:
;  none
;
; KEYWORD PARAMETERS:
;   delta_l_over_l : together with /log_bin, defines the relative width of each
;   bin
;   log_bin : set to 1 to build a logarithmic binning (default : 0, linear
;   binning)
;  delta_l_max : max width of output bins
;
; OUTPUTS:
;   map_cmn : 2D map of power spectrum values for each mode k(m,n)
;
; OPTIONAL OUTPUTS:
;  bin_width : width of bins
; float : to allow for non integer bin boundaries
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
;    Apr. 2011, N. ponthieu: removed lbin output for POKER.
;-

if n_params() lt 1 then begin
    message, /info, "call is:"
    print, "make_bintab, l, delta_l_min, bintab, delta_l_over_l=delta_l_over_l, $"
    print, "             log_bin=log_bin, delta_l_max = delta_l_max, $"
    print, "             bin_width = bin_width, float=float"
    return
endif

if keyword_set(delta_l_over_l) or keyword_set( log_bin) then begin
    if keyword_set(delta_l_over_l)*keyword_set(log_bin) eq 0 then begin
        message, /info, "delta_l_over_l and /log_bin must be set together"
        stop
    endif
endif

lmax      = double( max(l))
lmin      = double( min(l))
bintab    = [lmin]
bin_width = [0]
delta_l   = 0.0d0
l1        = lmin

dll = double( delta_l_over_l)

while (l1 + delta_l) lt lmax do begin
   
   if keyword_set(log_bin) then begin
      
      if keyword_set(float) then begin
         delta_l = ((l1*dll) > delta_l_min) < (lmax - l1)
      endif else begin
         delta_l = (long( (l1*dll) > delta_l_min)) < (lmax - l1)
      endelse
      
   endif else begin
      delta_l = delta_l_min
   endelse
   
   if keyword_set(delta_l_max) then delta_l = delta_l < delta_l_max
   
   l1        = l1 + delta_l
   bintab    = [bintab, l1]
   bin_width = [bin_width, delta_l]
endwhile

if max( bintab) lt max(l) then begin
;;    bin_width = [bin_width, max(l) - max(bintab)]
;;    bintab = [bintab, max(l)]

   ;; extend the last bin instead of adding a new small one
   n = n_elements( bintab)
   bintab[n-1]    = max(l)
   bin_width[n-1] = bintab[n-1]-bintab[n-2]

endif
bin_width = bin_width[1:*]

;bintab2lbin, bintab, lbin

end
