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




pro params2ascii, params, parfile, html=html

if n_params() lt 1 then begin
   print, "Calling sequence: "
   print, "params2ascii, params, parfile"
   return
endif

tags = tag_names( params)

if keyword_set(html) then begin
   suff = "<br>"
endif else begin
   suff = ""
endelse

openw, 1, parfile
for i=0, n_elements(tags)-1 do printf, 1, strlowcase(tags[i])+' = '+strtrim( params.(i),2)+" "+suff
close, 1

end
