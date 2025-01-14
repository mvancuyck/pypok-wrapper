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




; Bin a power spectrum square map into one dimensional power spectrum

pro cmn2cb, map_b, map_cmn, cb, sigma_b, sigma_b_1

if n_params() lt 1 then begin
   message, /info, "Calling sequence:"
   print, "cmn2cb, map_b, map_cmn, cb"
   return
endif

w  = where( map_b ne !undef)
m1 = map_b[w]
m2 = map_cmn[w]

h         = histogram( m1, bin=1.0d0, reverse_ind=R)
cb        = double(h*0.0d0)
sigma_b   = double(h*0.0d0)
sigma_b_1 = double(h*0.0d0)
for i=0, n_elements(h)-1 do begin
   ;; IF r[i] NE r[i+1] then begin
   ;;    cb[i] = avg( m2[R[R[i] : R[i+1]-1]])
   ;;    sigma_b_1[i] = stddev( m2[R[R[i] : R[i+1]-1]])
   ;;    sigma_b[i]   = stddev( m2[R[R[i] : R[i+1]-1]])/sqrt((r[i+1]-r[i]))
   ;; endif

   if (R[i+1]-R[i]) ge 2 then begin
      cb[i]        = avg(    m2[ R[ R[i] : R[i+1]-1]])
      sigma_b_1[i] = stddev( m2[ R[ R[i] : R[i+1]-1]])
      sigma_b[i]   = stddev( m2[ R[ R[i] : R[i+1]-1]])/sqrt((R[i+1]-R[i]))
   endif

   if (R[i+1]-R[i]) eq 1 then begin
      cb[i]        = m2[ R[ R[i]]]
      sigma_b_1[i] = !values.f_nan
      sigma_b[i]   = !values.f_nan
   endif

endfor

end
