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

program poker_mbb
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
  real(DP),     dimension(:,:), allocatable :: mtt_bb_loc, mtt_bb
  !real(DP),     dimension(:,:), allocatable :: mtt_bb_loc_1, mtt_bb_loc_x, mtt_bb_1, mtt_bb_x
  real(DP),     dimension(:), allocatable :: nb, nb_loc, bintab

  real(DP), dimension(0:1) :: kv, k1v, x, y
  real(DP) :: cosphi, sinphi, cos2phi, sin2phi, cos2phi2, sin2phi2, sin4phi, sigma, bk

  real(DP),     dimension(:,:),   allocatable :: pix_mn, map_k_binning, map_b
  complex(DPC), dimension(:,:), allocatable :: w_mn!, beam_mn, beam1_mn!, kernel
  real(DP),     dimension(:,:), allocatable :: beam_mn, beam1_mn
  real(DP) :: err, ww, kx, ky, k, k1, p, q, p1, kx1, ky1, kk, r, rm1, dummy_min, dummy_max
  complex(DPC) :: ww_c

  character(len=2) :: si, sj
  character(len=30) :: code
  character(len=200) :: file

  real(SP)     :: time0, time1

  !fftw parameter
  integer(I8B) :: plan
  !-------------------------------------------------------------------------------------------------

!  INCLUDE 'fftw3.f'
  INCLUDE '/home/nponthieu/Soft/fftw-3.3.8/api/fftw3.f'

#ifdef PARA
  INCLUDE 'mpif.h'
  CALL mpi_init( ierr)
  CALL mpi_comm_size( MPI_COMM_WORLD, nproc, ierr)
  CALL mpi_comm_rank( MPI_COMM_WORLD, myid,  ierr)
  if (myid == 0) then
     print*, "parallel version, myid=", myid
     print*, "#proc nproc=", nproc
  endif
#else
  myid = 0
  nproc = 1
#endif

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

  allocate( nb_loc(0:nbins-1), stat=ios)
  call assert_alloc( ios, code, "nb_loc")

  allocate( bintab(0:nbins+1-1), stat=ios)
  call assert_alloc( ios, code, "bintab")

  allocate( mtt_bb_loc(0:nbins-1,0:nbins-1), stat=ios)
  call assert_alloc( ios, code, "mtt_bb_loc")
  mtt_bb_loc = 0.0d0 !init

  allocate( w_mn( 0:naxis1-1, 0:naxis2-1), stat=ios)
  call assert_alloc( ios, code, "w_mn")
  w_mn = 0.0d0

  allocate( beam_mn(0:naxis1-1, 0:naxis2-1), stat=ios)
  call assert_alloc( ios, code, "beam_mn")
  beam_mn = 0.0d0

  allocate( beam1_mn(0:naxis1-1, 0:naxis2-1), stat=ios)
  call assert_alloc( ios, code, "beam1_mn")
  beam1_mn = 0.0d0

  allocate( pix_mn(0:naxis1-1, 0:naxis2-1), stat=ios)
  call assert_alloc( ios, code, "pix_mn")
  pix_mn = 0.0d0

  allocate( map_k_binning( 0:naxis1-1, 0:naxis2-1), stat=ios)
  call assert_alloc( ios, code, "map_k_binning")
  map_k_binning = 0.0d0

  if( myid == 0) then
     
     ! Useful quantity to bin 2D power spectra
     allocate( map_b( 0:naxis1-1, 0:naxis2-1), stat=ios)
     call assert_alloc( ios, code, "map_b")
     map_b = 0.0d0

     !Read bin boundaries
     call read_fits_1d( trim(pars%input_bintab), nbins+1, bintab)

     !Read binning weights
     call readimage( trim(pars%map_k_binning), naxis1, naxis2, map_k_binning)
     if (pars%verb /= 0) print*, "read "//trim(pars%map_k_binning)//" done."

     !Read mask (use pix_mn array to save memory)
     call readimage( trim(pars%mask), naxis1, naxis2, pix_mn)
     if (pars%verb /= 0) print*, "read "//trim(pars%mask)//" done."
     w_mn = pix_mn !passing R8 values to complex array

     ! Fourier transform the mask
     call dfftw_plan_dft_2d( plan, naxis1, naxis2, w_mn, w_mn, -1, fftw_estimate)
     call dfftw_execute_dft( plan, w_mn, w_mn)
     call dfftw_destroy_plan(plan) 

     !Correct for FFTW normalization
     w_mn = w_mn/npix

     !read beam transfer function (use pix_mn array to save memory)
     call readimage( trim(pars%beam), naxis1, naxis2, pix_mn)
     if (pars%verb /= 0) print*, "read "//trim(pars%beam)//" done."
     beam_mn = pix_mn !passing R8 values to complex array

     call readimage( trim(pars%beam1), naxis1, naxis2, pix_mn)
     if (pars%verb /= 0) print*, "read "//trim(pars%beam1)//" done."
     beam1_mn = pix_mn !passing R8 values to complex array

  endif

#ifdef PARA
  call MPI_Barrier( MPI_COMM_WORLD, ierr)
  call MPI_Bcast( bintab,        nbins+1, MPI_Double_Precision, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST( w_mn,          npix,      MPI_double_complex,   0, MPI_COMM_WORLD, ierr)
  !call MPI_BCAST( beam_mn,       npix,      MPI_double_complex,   0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST( beam_mn,       npix,      MPI_double_Precision,   0, MPI_COMM_WORLD, ierr)
  !call MPI_BCAST( beam1_mn,   npix,      MPI_double_complex,   0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST( beam1_mn,   npix,      MPI_double_Precision,   0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST( map_k_binning, npix,      MPI_double_precision, 0, MPI_COMM_WORLD, ierr)
  call MPI_Barrier( MPI_COMM_WORLD, ierr)
#endif

  if ((myid == 0) .and. (pars%verb /= 0)) print*, "MPI_Bcast done"

  ! Square pixel window function
  pix_mn = 1.0d0 !init to default
  if (pars%include_pix_window_function /= 0) then
     do i=0, npix-1
        call pix2mn( i, naxis1, naxis2, m, n)
        call mn2k(  m, n, naxis1, naxis2, kx, ky, k)
        if ( (kx /= 0.) .AND. (ky /= 0.)) then
           pix_mn(m,n) = dsin( PI*kx)*dsin( PI*ky) / (PI*kx * PI*ky)
        endif
     enddo
  endif

  ! Distribute pixels amond CPUs
  npix_proc  = npix/nproc !integer
  i_start = myid*npix_proc
  if (myid /= nproc-1) then
     i_end = i_start + npix_proc -1
  else
     i_end = npix-1
  endif

  !----------------------------------------- Main loop -----------------------------------------
  nb_loc       = 0.0d0 ! Init normalization factor in binning matrix P

  if ((myid == 0) .and. (pars%verb /= 0)) print*, "Starting main loop"
  if (pars%simul_type == 1) then

     do i=i_start, i_end
        call pix2mn( i, naxis1, naxis2, m, n)
        call mn2k( m, n, naxis1, naxis2, kx, ky, k)
        call k2bin( k, nbins, bintab, b)

        if (b >= 0) then
           p = map_k_binning( m, n)
           
           nb_loc(b) = nb_loc(b) + 1.0d0

           do j=0, npix-1
              call pix2mn( j, naxis1, naxis2, m1, n1)
              call mn2k( m1, n1, naxis1, naxis2, kx, ky, k)
              call k2bin( k, nbins, bintab, b1)
              
              if (b1 >= 0) then
                 
                 q = 1.0d0/map_k_binning(m1,n1)

                 if ( (m1 <= m) .and. (n1 <= n)) then
                    ww = abs( pix_mn(m1,n1) * sqrt( beam_mn( m1, n1) * beam1_mn( m1, n1)) * w_mn( m-m1,              n-n1))**2
                 endif
                 if ( (m1 <= m) .and. (n1 >  n)) then
                    ww = abs( pix_mn(m1,n1) * sqrt( beam_mn( m1, n1) * beam1_mn( m1, n1)) * w_mn( m-m1,       naxis2+n-n1))**2
                 endif
                 if ( (m1 >  m) .and. (n1 <= n)) then
                    ww = abs( pix_mn(m1,n1) * sqrt( beam_mn( m1, n1) * beam1_mn( m1, n1)) * w_mn( naxis1+m-m1,       n-n1))**2
                 endif
                 if ( (m1 >  m) .and. (n1 >  n)) then
                    ww = abs( pix_mn(m1,n1) * sqrt( beam_mn( m1, n1) * beam1_mn( m1, n1)) * w_mn( naxis1+m-m1,naxis2+n-n1))**2
                 endif
                 mtt_bb_loc(b1,b) = mtt_bb_loc(b1,b) + p*ww*q

              endif !b1>=0
           enddo!j=0, npix-1
        endif !b>=0
     enddo !i=i_start, i_end

  else
     print*, "poker_mbb : simul_type /= 1"
     print*, "polarization not implemented yet"
     stop
  endif!simul_type==1

  if ((myid == 0) .and. (pars%verb /= 0)) print*, "Main loop done."

  if (myid == 0) then
     allocate( mtt_bb(0:nbins-1,0:nbins-1), stat=ios)
     call assert_alloc( ios, code, "mtt_bb")
     mtt_bb = 0.0d0

     allocate( nb(0:nbins-1), stat=ios)
     call assert_alloc( ios, code, "nb")
     nb = 0.0d0
  endif

#ifdef PARA
  call MPI_Barrier( MPI_COMM_WORLD, ierr)
  call MPI_Reduce( mtt_bb_loc, mtt_bb, nbins**2,     MPI_Double_Precision, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce( nb_loc,     nb,     nbins,        MPI_Double_Precision, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Barrier( MPI_COMM_WORLD, ierr)
#else
  mtt_bb    = mtt_bb_loc
  nb        = nb_loc
#endif

  if ((myid == 0) .and. (pars%verb /= 0)) print*, "Reduce done."

  if (myid == 0) then

     !Compute map_b
     do i=0, npix-1
        call pix2mn( i, naxis1, naxis2, m, n)
        call mn2k( m, n, naxis1, naxis2, kx, ky, k)
        call k2bin( k, nbins, bintab, b)
        map_b(m,n) = b
     enddo

     !Output
     call writeimage( trim(pars%file_map_b), naxis1, naxis2, map_b)

     ! Account for normalization factor in P (number of modes per bin, see Hivon et al, Eq. (20))
     do b=0, nbins-1
        if (nb(b) /= 0.) then
           mtt_bb(:,b) = mtt_bb(:,b)/nb(b)
        endif
     enddo

     !Output
     call writeimage( trim(pars%file_mtt_bb_x), nbins, nbins, mtt_bb)
     if ((myid == 0) .and. (pars%verb /= 0)) print*, "wrote mtt_bb (x)"

     call cpu_time(time1)
     !!open(  2, file=trim(pars%report_card), status='unknown')
     !!write( 2, *) "npix  = ", npix
     !!write( 2, *) "nproc = ", nproc
     !!write( 2, *) "CPU time (sec): ", time1 - time0
     !!close( 2)

     if (myid ==0) print*, "CPU time (sec): ", time1 - time0

  endif

#ifdef PARA
  CALL mpi_finalize(ierr)
#else

  deallocate( nb_loc)
  deallocate( bintab)
  deallocate( mtt_bb_loc)
  deallocate( w_mn)
  deallocate( beam_mn)
  deallocate( pix_mn)
  deallocate( map_k_binning)
  deallocate( mtt_bb)
  deallocate( nb)

#endif

end program poker_mbb
