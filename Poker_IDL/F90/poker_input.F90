! Copyright 2010, Nicolas Ponthieu, Julien Grain, Guilaine Lagache
! 
! This file is part of Poker.
! 
! Poker is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! Poker is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with Poker.  If not, see <http://www.gnu.org/licenses/>.
! 
! For more information about Poker, see http://www.ias.u-psud.fr/poker
! ======================================================================




module poker_input


contains

subroutine poker_read_params(i_parfile, myid, pars, parafile)
use healpix_types
use paramfile_io
use extension
use poker_types
use misc_utils
implicit none

type(paramfile_handle) :: handle
type(poker_str), intent(inout) :: pars
character(len=200), intent(inout) :: parafile
integer, intent(in) :: i_parfile, myid
call getArgument(i_parfile, parafile)
!if (myid == 0) then
!   handle = parse_init(parafile)
!else
   handle = parse_init( parafile, silent=.true.)
!endif

pars%map = parse_string( handle, "map", "map.fits") ! data maps for check-kernel
pars%map_k_binning = parse_string( handle, "map_k_binning", "map_k_binning") ! typically map_k^bin_exp or anyother function
pars%res_pix = parse_double( handle, "res_pix",  0.00029088821d0 ) ! map resolution (radians) (1arcmin default)
pars%report_card = parse_string( handle, "report_card", "poker_report_card.txt") ! short execution summary
pars%mask = parse_string( handle, "mask", "mask.fits") ! Mask real weights
pars%nbeams = parse_int( handle, "nbeams",  1 ) ! set to 2 if the two cross-correlated maps have different beams (beam and beam1)
pars%beam = parse_string( handle, "beam", "beam.fits") ! beam transfer function (real numbers)
pars%beam1 = parse_string( handle, "beam1", "beam1.fits") ! beam transfer function (real numbers) of the second map for cross-correlations
pars%file_map_b = parse_string( handle, "file_map_b", "map_b.fits") ! bin addresses
pars%file_mtt_bb = parse_string( handle, "file_mtt_bb", "poker_mtt_bb_out.fits") ! mode mixing matrix for temperature only
pars%file_mtt_bb_1 = parse_string( handle, "file_mtt_bb_1", "poker_mtt_bb_out_1.fits") ! mode mixing matrix for temperature only (beam1)
pars%file_mtt_bb_x = parse_string( handle, "file_mtt_bb_x", "poker_mtt_bb_out_x.fits") ! mode mixing matrix for temperature only (cross-beam)
pars%simul_type = parse_int( handle, "simul_type",  1 ) ! 1 for temperature only, 2 for temperature and polarization
pars%input_bintab = parse_string( handle, "input_bintab", "bintab.fits") ! bin table for poker_mbb.f90
pars%include_pix_window_function = parse_int( handle, "include_pix_window_function",  0 ) ! set to 1 to include pixel window function in Mbb (warning : be consistent in cls2map)
pars%verb = parse_int( handle, "verb",  0 ) ! verbose
pars%n_cpu_max = parse_int( handle, "n_cpu_max",  100 ) ! maximum number of CPUs to be used
! #IDL parameters and keywords
pars%scale = parse_double( handle, "scale",  1.0d0 ) ! scale for zero padding
pars%delta_l_over_l = parse_double( handle, "delta_l_over_l",  0.1d0 ) ! binning width
pars%beta = parse_double( handle, "beta",  1.d0 ) ! binning index
pars%log_binning = parse_int( handle, "log_binning",  1 ) ! set to 1 to have a linear binning
pars%dir = parse_string( handle, "dir", ".") ! output directory
pars%outparfile = parse_string( handle, "outparfile", "poker_out.par") ! copy of the parameter file passed from IDL to F90
pars%nx = parse_int( handle, "nx",  100 ) ! size of obs patch in pixels
pars%ny = parse_int( handle, "ny",  100 ) ! size of obs patch in pixels
pars%nx_large = parse_int( handle, "nx_large",  100 ) ! size of obs patch in pixels
pars%ny_large = parse_int( handle, "ny_large",  100 ) ! size of obs patch in pixels
pars%patch = parse_string( handle, "patch", "patch.fits") ! obs patch
pars%remove_1st_bin = parse_int( handle, "remove_1st_bin",  0 ) ! set to 1 to remove the 1st bin of Mbb to improve its conditionning (DC level)
pars%k_out = parse_string( handle, "k_out", "k_out.fits") ! output frequency bins
pars%keep_avg = parse_int( handle, "keep_avg",  0 ) ! set to 1 to keep the map average
pars%k_nyquist = parse_double( handle, "k_nyquist",  0.0d0 ) ! Nyquist mode.
pars%apod_length = parse_int( handle, "apod_length",  0 ) ! apodization length in number of pixels
pars%dk_min = parse_double( handle, "dk_min",  -1.0d0 ) ! minimum width of frequency bins


end subroutine poker_read_params

end module poker_input
