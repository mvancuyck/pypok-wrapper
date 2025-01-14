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



pro poker_types, str


;Init str
str = create_struct("date", '111022_23:33:52')
str = create_struct( str, "map", "map.fits")
str = create_struct( str, "map_k_binning", "map_k_binning.fits")
str = create_struct( str, "res_pix",    0.00029088821d0 )
str = create_struct( str, "report_card", "poker_report_card.txt")
str = create_struct( str, "mask", "mask.fits")
str = create_struct( str, "nbeams",    1 )
str = create_struct( str, "beam", "beam.fits")
str = create_struct( str, "beam1", "beam1.fits")
str = create_struct( str, "file_map_b", "map_b.fits")
str = create_struct( str, "file_mtt_bb", "poker_mtt_bb_out.fits")
str = create_struct( str, "file_mtt_bb_1", "poker_mtt_bb_out_1.fits")
str = create_struct( str, "file_mtt_bb_x", "poker_mtt_bb_out_x.fits")
str = create_struct( str, "simul_type",    1 )
str = create_struct( str, "input_bintab", "bintab.fits")
str = create_struct( str, "include_pix_window_function",    0 )
str = create_struct( str, "verb",    0 )
str = create_struct( str, "n_cpu_max",    100 )
str = create_struct( str, "scale",    1.0d0 )
str = create_struct( str, "delta_l_over_l",    0.1d0 )
str = create_struct( str, "beta",    0.d0 )
str = create_struct( str, "log_binning",    1 )
str = create_struct( str, "dir", ".")
str = create_struct( str, "outparfile", "poker_out.par")
str = create_struct( str, "nx",    100 )
str = create_struct( str, "ny",    100 )
str = create_struct( str, "nx_large",    100 )
str = create_struct( str, "ny_large",    100 )
str = create_struct( str, "patch", "patch.fits")
str = create_struct( str, "remove_1st_bin",    0 )
str = create_struct( str, "k_out", "k_out.fits")
str = create_struct( str, "keep_avg",    0 )
str = create_struct( str, "k_nyquist",    0.0d0 )
str = create_struct( str, "apod_length",    0 )
str = create_struct( str, "dk_min",    -1.0d0 )


end
