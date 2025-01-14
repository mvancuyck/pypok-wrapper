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

pro poker_clean, params

if n_params() lt 1 then begin
   poker_types, params
endif

;; test file existence to avoid anoying error messages from the system
if file_test(params.map            ) then spawn, "\rm "+params.map
if file_test(params.map_k_binning  ) then spawn, "\rm "+params.map_k_binning
if file_test(params.report_card    ) then spawn, "\rm "+params.report_card
if file_test(params.mask           ) then spawn, "\rm "+params.mask
if file_test(params.beam           ) then spawn, "\rm "+params.beam
if file_test(params.beam1          ) then spawn, "\rm "+params.beam1
if file_test(params.file_map_b     ) then spawn, "\rm "+params.file_map_b
if file_test(params.file_mtt_bb    ) then spawn, "\rm "+params.file_mtt_bb
if file_test(params.file_mtt_bb_1  ) then spawn, "\rm "+params.file_mtt_bb_1
if file_test(params.file_mtt_bb_x  ) then spawn, "\rm "+params.file_mtt_bb_x
if file_test(params.input_bintab   ) then spawn, "\rm "+params.input_bintab
if file_test(params.outparfile     ) then spawn, "\rm "+params.outparfile
if file_test(params.patch          ) then spawn, "\rm "+params.patch
if file_test(params.k_out          ) then spawn, "\rm "+params.k_out          

end
