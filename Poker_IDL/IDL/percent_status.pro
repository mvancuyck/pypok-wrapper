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



pro percent_status, imc, nmc, step

if n_params() lt 1 then begin
   message, /info, "Calling sequence: "
   print, "percent_status, imc, nmc, step"
   return
endif


if long(imc) mod (long(nmc)/long(step)) eq 0 then print, strtrim( long(100.*imc/float(nmc)),2)+"%"

end
