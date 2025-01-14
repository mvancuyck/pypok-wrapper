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



program poker_give_map_b
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

  integer(I4B) :: i, j, ios, n_iter, npix, b, b1, m, n, m1, n1, np, ii, jj, i_alpha, j_alpha, itap, nbins
  integer(I4B) :: myid, nproc, ierr, npix_proc, i_start, i_end, naxis1, naxis2
  real(DP),     dimension(:,:), allocatable :: mtt_bb_loc, mtt_bb, mpp_bb_loc, mpp_bb, mtp_bb_loc, mtp_bb, cu_t
  real(DP),     dimension(:), allocatable :: nb, nb_loc, bintab, pseudo_cb_loc, pseudo_cb

  real(DP), dimension(0:1) :: kv, k1v, x, y
  real(DP) :: cosphi, sinphi, cos2phi, sin2phi, cos2phi2, sin2phi2, sin4phi, sigma, bk

  real(DP),     dimension(:), allocatable :: temp
  real(DP),     dimension(:,:),   allocatable :: pix_mn, map_k_binning, map_b, map_b_loc
  complex(DPC), dimension(:,:), allocatable :: w_mn, beam_mn!, kernel
  real(DP) :: err, ww, kx, ky, k, k1, p, q, p1, kx1, ky1, kk, r, rm1
  complex(DPC) :: ww_c

  character(len=2) :: si, sj
  character(len=30) :: code
  character(len=200) :: file

  real(SP)     :: time0, time1

  !fftw parameter
  integer(I8B) :: plan
  !-------------------------------------------------------------------------------------------------

  INCLUDE 'fftw3.f'

!!#ifdef PARA
!!  INCLUDE 'mpif.h'
!!  CALL mpi_init( ierr)
!!  CALL mpi_comm_size( MPI_COMM_WORLD, nproc, ierr)
!!  CALL mpi_comm_rank( MPI_COMM_WORLD, myid,  ierr)
!!#else
  myid = 0
  nproc = 1
!!#endif

  code = "poker_mbb"

  if (myid == 0) call cpu_time(time0)

  i_parfile = 1
  call poker_read_params(i_parfile, myid, pars, parafile)
  if ((myid == 0) .and. (pars%verb /= 0)) print*, "input parameters done"

  !get map size parameters
  call get_naxis12( trim(pars%mask), naxis1, naxis2)
  npix = naxis1 * naxis2

  !get binning parameters
  call get_naxis1( trim(pars%input_bintab), nbins)
  nbins = nbins - 1 !bintab has nbins+1 elements

  allocate( bintab(0:nbins+1-1), stat=ios)
  call assert_alloc( ios, code, "bintab")

  allocate( map_b( 0:naxis1-1, 0:naxis2-1), stat=ios)
  call assert_alloc( ios, code, "map_b")
  map_b = 0.0d0

  !Read bin boundaries
  call read_fits_1d( trim(pars%input_bintab), nbins+1, bintab)

  !Compute map_b
  do i=0, npix-1

     call pix2mn( i, naxis1, naxis2, m, n)
     call mn2k( m, n, naxis1, naxis2, kx, ky, k)
     call k2bin( k, nbins, bintab, b)
     
     map_b(m,n) = b
  enddo

  !call writeimage( trim("/tmp/map_b.fits"), naxis1, naxis1, map_b)
  call writeimage( trim(pars%file_map_b), naxis1, naxis2, map_b)


!!#ifdef PARA
!!  CALL mpi_finalize(ierr)
!!#endif

  deallocate( bintab)
  deallocate( map_b)

end program poker_give_map_b
