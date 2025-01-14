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



pro mc_reduce, tab_res, tab_avg, sigma_avg, cov_mat, xcorr, quick=quick

;+
; NAME:
;   mc_reduce
;
; PURPOSE:
;   Computes average, stddev, covariance and correlations on one dimension of a
;   2D array.
;
; CALLING SEQUENCE:
;  mc_reduce, tab_res, tab_avg, sigma_avg, cov_mat, xcorr, quick=quick
;
; INPUTS:
;     tab_res : array of values (2D)
;
; OPTIONAL INPUTS:
;  none
;
; KEYWORD PARAMETERS:
;  quick : to skip covariance and correlation estimations
;
; OUTPUTS:
;    tab_avg : average of input values on the 1st dimension
;  sigma_avg : standard deviation associated to tab_avg
;    cov_mat : covariance matrix between elements of tab_avg
;      xcorr : correlation matrix between elements of tab_Avg
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
   message, /info, "Call is:"
   print, "mc_reduce, tab_res, tab_avg, sigma_avg, cov_mat, xcorr"
   return
endif

nmc   = n_elements( tab_res[*,0])
nbins = n_elements( tab_res[0,*])

; Compute average 
tab_avg = avg( tab_res, 0)

; Stddev
sigma_avg = tab_avg*0.
for i=0, n_elements(tab_res[0,*])-1 do sigma_avg[i] = stddev( tab_res[*,i])/sqrt( nmc)



if not keyword_set(quick) then begin
; Compute covariance matrix
   cov_mat = dblarr( nbins, nbins)
   for b=0, nbins-1 do begin
      for b1=0, nbins-1 do begin
         cov_mat[b,b1] = avg( (tab_res[*,b] - tab_avg[b])*(tab_res[*,b1] - tab_avg[b1]))
      endfor
   endfor

;Cross-correlations
   xcorr = cov_mat*0.0d0
   for b=0, nbins-1 do begin
      for b1=0, nbins-1 do begin
         xcorr[b,b1] = cov_mat[b, b1]/sqrt( cov_mat[b,b]*cov_mat[b1,b1])
      endfor
   endfor
endif


end
