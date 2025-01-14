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





; Generate power spectrum flat map (2D) from input binned power spectrum (1D)

pro cb2cmn, map_b, cb, map_k_binning, map_cmn

if n_params() lt 1 then begin
   message, /info, "Calling sequence: "
   print, "cb2cmn, map_b, cb, map_k_binning, map_cmn"
   return
endif


map_cmn = map_b * 0.0d0

for ib=0, n_elements(cb)-1 do begin
   w = where( map_b eq ib, nw)
   if nw ne 0 then begin
      map_cmn[w] = 1.0d0/map_k_binning[w] * cb[ib]
   endif
endfor

end
