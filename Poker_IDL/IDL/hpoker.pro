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


pro hpoker, map, res_arcmin, x1_list, x2_list, y1_list, y2_list, $
            k_out, pk_out, sigma_pk_out, $
            pb_out=pb_out, mask=mask, map1=map1, $
            a_beam=a_beam, b_beam=b_beam, $
            in_params=in_params, out_params=out_params, $
            in_arrays=in_arrays, out_arrays=out_arrays, $
            bypass_mtt_bb=bypass_mtt_bb, $
            remove_1st_bin=remove_1st_bin, log_binning=log_binning, $
            dir=dir, delta_l_over_l=delta_l_over_l, $
            nx_large=nx_large, ny_large=ny_large, apod_length=apod_length, $
            dk_min=dk_min, beta=beta, bintab_in=bintab_in, $
            noise_pseudo_spec=noise_pseudo_spec, $
            include_pix_window_function=include_pix_window_function, clean_poker=clean_poker


message, /info, "See check_hpoker for now"            
end
