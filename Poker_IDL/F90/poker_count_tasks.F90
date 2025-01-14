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



program poker_count_tasks
  use healpix_types
  use extension
  use misc_utils

  use poker_types
  use poker_input
  use poker_tools

  implicit none

  character(len=200) :: parafile
  integer :: i_parfile
  type(poker_str) :: pars

  integer(I4B) :: i, j, npix, nproc, npix_proc, ipix_start, ipix_end, naxis1, naxis2
  integer(I4B) :: n_ops_avg, n_ops_min, ios, r, n_ops, nproc_opt
  character(len=80) :: code
  character(len=200) :: file

  code = "poker_count_tasks"
  i_parfile = 1
  call poker_read_params(i_parfile, 0, pars, parafile)

  call get_naxis12( trim(pars%mask), naxis1, naxis2)
  npix = naxis1 * naxis2

  n_ops_min = npix !init
  do nproc = 1, pars%n_cpu_max
     npix_proc = npix/nproc
     r = npix - npix_proc*nproc

     n_ops = npix_proc + r
     print*, "nproc, npix_proc, r, n_ops: ", nproc, npix_proc, r, n_ops

     if (n_ops < n_ops_min) then
        nproc_opt = nproc
        n_ops_min = n_ops
     endif
  enddo


  print*, "nproc_opt, n_ops_min: ", nproc_opt, n_ops_min

end program poker_count_tasks
