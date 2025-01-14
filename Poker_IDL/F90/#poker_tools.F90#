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



module poker_tools

  USE healpix_types
  use misc_utils
  IMPLICIT none

contains

subroutine dot_product( n, v1, v2, res)
  implicit none
  integer(I4B), intent(in) :: n
  real(DP), dimension(0:n-1), intent(in) :: v1, v2
  real(DP), intent(out) :: res
  integer(I4B) :: i
  res = 0.0d0
  do i=0, n-1 
     res = res + v1(i)*v2(i)
  enddo
end subroutine dot_product

!========================================================================================
subroutine bcg( n, a1, b, itmax, tol, x, n_iter, err, verb)
  implicit none
!INPUT:
!n     : matrix a is n x n
!a1    : matrix a in vector form
!b     : the routine solves Ax = b
!itmax : maximum number of iterations
!tol   : tolerance on convergence
!
!OUTPUT:
!x      : output solution of the system
!n_iter : number of iterations actually performed
!err    : actual error on the solution

  !input/output parameters
  integer(I4B), intent(in) :: n, itmax, verb
  real(DP), dimension(0:n*n-1), intent(in) :: a1
  real(DP), dimension(0:n-1), intent(in) :: b
  real(DP), intent(in) :: tol

  real(DP), dimension(0:n-1), intent(out) :: x
  integer(I4B), intent(out) :: n_iter
  real(DP), intent(out) :: err


  !Local variables
  character(len=3) :: code
  integer(I4B) :: n2, i, j
  real(DP), dimension(0:n-1) :: r, r_bar, p, p_bar, res
  real(DP), dimension(0:n-1,0:n-1) :: a

  real(DP) :: alpha, num, denom, num1, denom1, beta

  code = "bcg" 
  n2 = n*n

  x = 0.0d0 !initial guess
  
  do i=0, n-1
     do j=0, n-1
        a(i,j) = a1(n*i + j)
     enddo
  enddo

  !Initial values
  x = 0.0d0 !initial guess
  
  r = b !- matmul(a,x)
  p = r
  
  r_bar = r
  p_bar = p

  n_iter = 0
  err = 2*tol !init to larger value
  do while( (n_iter <= itmax) .and. (err > tol))

     call dot_product( n, r, r_bar, num1)
     call dot_product( n, p_bar, matmul(a,p), denom1)
     alpha = num1/denom1

     r     = r - alpha*matmul(a,p)
     r_bar = r_bar - alpha*matmul( transpose(a), p_bar)

     call dot_product( n, r_bar, r, num)
     beta = num/num1

     x     = x + alpha*p
     p     = r + beta*p
     p_bar = r_bar + beta*p_bar

     !update err and iter number
     res = b - matmul(a,x)
     call dot_product( n, res, res, err)

     if (verb /= 0) print*, "iter, err: ", n_iter, err
     n_iter = n_iter + 1
  enddo

end subroutine bcg

!======================================================================================================
subroutine pix2mnpq( ipix, nx, ny, nbins, bintab, map_k_binning, m, n, kx, ky, k, p, q)
  use healpix_types
  implicit none

  integer(I4B), intent(in)                        :: ipix, nx, ny, nbins
  
  !real(DP),     intent(in)                        :: bin_exp
  real(DP),     dimension(0:nx*ny-1) , intent(in) :: map_k_binning

  real(DP),     intent(in), dimension(0:nbins)    :: bintab !size of bintab = nbins+1

  integer(I4B), intent(out)                       :: m, n!, b
  real(DP),     intent(out)                       :: kx, ky, k, p, q

  integer(I4B) :: i

  n = ipix/ny
  m = ipix - n*ny
  if (m <= nx/2) then
     kx = dble(m)/nx
  else
     kx = dble(nx-m)/nx
  endif
  if (n <= ny/2) then
     ky = dble(n)/ny
  else
     ky = dble(ny-n)/ny
  endif
  k = sqrt( kx**2 + ky**2) !reduced units
  
  ! Determine binning weight and bin address
  p = map_k_binning(ipix)
  if (p /= 0.) then
     q = 1.0d0/p
  else
     q = 0.0d0
  endif

end subroutine pix2mnpq

!=========================================================
subroutine pix2mn( ipix, nx, ny, m, n)
  use healpix_types
  implicit none

  integer(I4B), intent(in)  :: ipix, nx, ny
  integer(I4B), intent(out) :: m, n

  n = ipix/nx
  m = ipix - n*nx

end subroutine pix2mn

!=========================================================
subroutine mn2k( m, n, nx, ny, kx, ky, k)
  use healpix_types
  implicit none

  integer(I4B), intent(in) :: m, n, nx, ny
  real(DP), intent(out)    :: kx, ky, k

  if (m <= nx/2) then
     kx = dble(m)/nx
  else
     kx = dble(nx-m)/nx
  endif
  if (n <= ny/2) then
     ky = dble(n)/ny
  else
     ky = dble(ny-n)/ny
  endif
  k = sqrt( kx**2 + ky**2) !reduced units
  
end subroutine mn2k

!===========================================================================
! Determines to which bin belongs a given k wavevector
subroutine k2bin( k, nbins, bintab, b)
  use healpix_types
  implicit none

  real(DP), intent(in) :: k
  integer(I4B), intent(in) :: nbins
  real(DP), intent(in), dimension(0:nbins) :: bintab !bintab is nbins+1 long
  integer(I4B), intent(out) :: b

  integer(I4B) :: i, test

  test = 0
  b    = 0
  do while( (test==0) .and. (b <nbins))
     if ( (bintab(b) <= k) .and. ( k < bintab(b+1))) test = 1
     b = b+1
  enddo

  !correct for extra addition and the end of "do while" loop
  b = b-1

  !if k is outside bintab bounds, set b to undef value
  if (test == 0) then
     b = -32768
  endif

end subroutine k2bin

!=========================================================================
subroutine get_naxis12( filename, naxis1, naxis2)
!! from cookbook.f in cfitstio library
  use healpix_types
  implicit none

  character(len=*), intent(in) :: filename
  integer(I4B), intent(out) :: naxis1, naxis2
  integer(I4B) ::  status,unit,readwrite,blocksize,nkeys,nspace,hdutype,i,j, naxes(2), naxis, nfound

  character(len=80) :: record, comment

  !  The STATUS parameter must always be initialized.
  status=0

  !  Get an unused Logical Unit Number to use to open the FITS file.
  call ftgiou(unit,status)

  !     open the FITS file, with read-only access.  The returned BLOCKSIZE
  !     parameter is obsolete and should be ignored. 
  readwrite=0
  call ftopen(unit,filename,readwrite,blocksize,status)
  if (status > 0) then
     call printerror(status)
     call ftclos(unit, status)
     return
  endif

  !     determines the presence of image
  call ftgkyj(unit,'NAXIS', naxis, comment, status)
  if (naxis <= 0) then
     print*, "cound not find NAXIS"
     call ftclos(unit, status)
     call fatal_error
  else
     !---------------------------------
     ! there is an image
     !---------------------------------
     !        determine the size of the image (look naxis1 and naxis2)
     call ftgknj(unit,'NAXIS',1_i4b,2_i4b,naxes,nfound,status)
     if (nfound /= 2) then
        print*, trim(filename)//" does not seem to contain a 2D image"
        call ftclos(unit, status)
        call fatal_error
     endif
     if (nfound < 1) then
        call printerror(status)
        print *,'can not find NAXIS1.'
        call ftclos(unit, status)
        call fatal_error
     endif
     
     naxis1 = naxes(1)
     naxis2 = naxes(2)
  endif
  
  !  The FITS file must always be closed before exiting the program. 
  !  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
  call ftclos(unit, status)
  call ftfiou(unit, status)
  
  !  Check for any error, and if so print out error messages.
  !  The PRINTERROR subroutine is listed near the end of this file.
  if (status .gt. 0)call printerror(status)
  
end subroutine get_naxis12

!=========================================================================
subroutine get_naxis1( filename, naxis1)
!to read a 1D array
!! from cookbook.f in cfitstio library
  use healpix_types
  implicit none

  character(len=*), intent(in) :: filename
  integer(I4B), intent(out) :: naxis1
!  integer(I4B) :: naxis2
  integer(I4B) ::  status,unit,readwrite,blocksize,nkeys,nspace,hdutype,i,j, naxes(2), naxis, nfound

  character(len=80) :: record, comment

  !  The STATUS parameter must always be initialized.
  status=0

  !  Get an unused Logical Unit Number to use to open the FITS file.
  call ftgiou(unit,status)

  !     open the FITS file, with read-only access.  The returned BLOCKSIZE
  !     parameter is obsolete and should be ignored. 
  readwrite=0
  call ftopen(unit,filename,readwrite,blocksize,status)
  if (status > 0) then
     call printerror(status)
     call ftclos(unit, status)
     return
  endif

  !     determines the presence of image
  call ftgkyj(unit,'NAXIS', naxis, comment, status)
  if (naxis <= 0) then
     print*, "cound not find NAXIS"
     call ftclos(unit, status)
     call fatal_error
  else
     call ftgknj(unit,'NAXIS',1_i4b,2_i4b,naxes,nfound,status)
     if (nfound < 1) then
        call printerror(status)
        print *,'can not find NAXIS1.'
        call ftclos(unit, status)
        call fatal_error
     endif
     
     naxis1 = naxes(1)
  endif
  
  !  The FITS file must always be closed before exiting the program. 
  !  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
  call ftclos(unit, status)
  call ftfiou(unit, status)
  
  !  Check for any error, and if so print out error messages.
  !  The PRINTERROR subroutine is listed near the end of this file.
  if (status .gt. 0)call printerror(status)
  
end subroutine get_naxis1

!============================================================
subroutine writeimage( filename, naxis1, naxis2, map)
!! cookbook.f in the Cfitsio distribution
!! see also http://heasarc.gsfc.nasa.gov/docs/software/fitsio/c/f_user/node8.html
  use healpix_types
  implicit none

!     Create a FITS primary array containing a 2-D image

  character(len=*), intent(in) :: filename
  integer(I4B),     intent(in) :: naxis1, naxis2
  real(DP), dimension(0:naxis1-1,0:naxis2-1), intent(in) :: map

  character(len=200) :: file_out
  integer(I4B) :: status,unit,blocksize,bitpix,naxis,naxes(2)
  integer(I4B) :: i,j,group,fpixel,nelements
  logical :: simple, extend

  status=0
  !     Name of the FITS file to be created:
  !filename='ATESTFILE.FITS'

  file_out = "!"//trim(filename)

  !     Get an unused Logical Unit Number to use to create the FITS file
  call ftgiou(unit,status)

  !     create the new empty FITS file
  blocksize=1
  call ftinit(unit,file_out,blocksize,status)

  !     initialize parameters about the FITS image (300 x 200 16-bit integers)
  !bitpix=16
  
  !32-bit double
  bitpix = -64
  simple=.true.
  naxis=2
  naxes(1) = naxis1 !7 !300
  naxes(2) = naxis2 !3 !200
  extend=.true.

  !     write the required header keywords
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  
  !     write the map to the FITS file
  group=1
  fpixel=1
  nelements=naxes(1)*naxes(2)
  !call ftpprj(unit,group,fpixel,nelements,map,status)
  call ftpprd(unit,group,fpixel,nelements,map,status)
  
  !     close the file and free the unit number
  call ftclos(unit, status)
  call ftfiou(unit, status)

end subroutine writeimage

!=================================================================
subroutine readimage( filename, naxis1, naxis2, map)
  !! cookbook.f in the Cfitsio distribution
  !! see also http://heasarc.gsfc.nasa.gov/docs/software/fitsio/c/f_user/node8.html
  use healpix_types
  implicit none
 
  !  Read a FITS image and determine the minimum and maximum pixel value.
  !  Rather than reading the entire image in
  !  at once (which could require a very large array), the image is read
  !  in pieces, 100 pixels at a time.  
  
  character(len=*), intent(in)  ::  filename
  integer(I4B),     intent(in)  :: naxis1, naxis2
  !real(DP),         intent(out) :: map(1:naxis1,1:naxis2)
  real(DP),         intent(out) :: map(0:naxis1-1,0:naxis2-1)

  integer(I4B) :: status,unit,readwrite,blocksize,naxes(2),nfound
  integer(I4B) :: group,firstpix,nbuffer,npixels,i,j,k, ipix
  real(DP)     :: datamin,datamax,nullval,buffer(100)
  !real(SP)     :: datamin,datamax,nullval,buffer(100)
  logical      :: anynull
 
  !  The STATUS parameter must always be initialized.
  status=0
 
  !  Get an unused Logical Unit Number to use to open the FITS file.
  call ftgiou(unit,status)
  
  !  Open the FITS file previously created by WRITEIMAGE
  !filename='ATESTFILEZ.FITS'
  readwrite=0
  call ftopen(unit,filename,readwrite,blocksize,status)
  
  !  Determine the size of the image.
  call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)
 
  !  Check that it found both NAXIS1 and NAXIS2 keywords.
  if (nfound .ne. 2)then
     print *,'READIMAGE failed to read the NAXISn keywords.'
     return
  end if
  
  !  Initialize variables
  npixels=naxes(1)*naxes(2)
  group=1
  firstpix=1
  nullval=-999
  datamin=1.0E30
  datamax=-1.0E30
 
  ! Fill output map
  map = -32768 !init
  do while (npixels .gt. 0)
     !         read up to 100 pixels at a time 
     nbuffer = min(100,npixels)
     
     !call ftgpve(unit,group,firstpix,nbuffer,nullval,buffer,anynull,status)
     call ftgpvd(unit,group,firstpix,nbuffer,nullval,buffer,anynull,status)
     
     !         find the min and max values
     do k=1,nbuffer

        !pixel indices
        ipix = ( firstpix- 1 + k-1)
        j = ipix/naxis1
        i = ipix - j*naxis1

        map(i,j) = buffer(k)
     end do
     
     !         increment pointers and loop back to read the next group of pixels
     npixels=npixels-nbuffer
     firstpix=firstpix+nbuffer
  end do

  !  The FITS file must always be closed before exiting the program. 
  !  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
  call ftclos(unit, status)
  call ftfiou(unit, status)
  
  !  Check for any error, and if so print out error messages.
  !  The PRINTERROR subroutine is listed near the end of this file.
  if (status .gt. 0)call printerror(status)
  
end subroutine readimage


!=================================================================
subroutine read_fits_1d( filename, naxis1, array)
  !! cookbook.f in the Cfitsio distribution
  !! see also http://heasarc.gsfc.nasa.gov/docs/software/fitsio/c/f_user/node8.html
  use healpix_types
  implicit none
 
  character(len=*), intent(in)  ::  filename
  integer(I4B),     intent(in)  :: naxis1
  real(DP),         intent(out) :: array(0:naxis1-1)

  integer(I4B) :: status,unit,readwrite,blocksize,naxes(2),nfound
  integer(I4B) :: group,firstpix,nbuffer,npixels,i,j,k, ipix
  real(DP)     :: datamin,datamax,nullval,buffer(100)
  logical      :: anynull
 
  !  The STATUS parameter must always be initialized.
  status=0
 
  !  Get an unused Logical Unit Number to use to open the FITS file.
  call ftgiou(unit,status)
  
  !  Open the FITS file previously created by WRITEIMAGE
  !filename='ATESTFILEZ.FITS'
  readwrite=0
  call ftopen(unit,filename,readwrite,blocksize,status)
  
  !  Determine the size of the image.
  call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)
 
  !  Initialize variables
  npixels=naxes(1)
  group=1
  firstpix=1
  nullval=-999
  datamin=1.0E30
  datamax=-1.0E30
 
  ! Fill output map
  array = -32768 !init
  do while (npixels .gt. 0)
     !         read up to 100 pixels at a time 
     nbuffer = min(100,npixels)
     
     call ftgpvd(unit,group,firstpix,nbuffer,nullval,buffer,anynull,status)
     
     !         find the min and max values
     do k=1,nbuffer
        ipix = ( firstpix-1 + k-1)
        array(ipix) = buffer(k)
     end do
     
     !         increment pointers and loop back to read the next group of pixels
     npixels=npixels-nbuffer
     firstpix=firstpix+nbuffer
  end do

  !  The FITS file must always be closed before exiting the program. 
  !  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
  call ftclos(unit, status)
  call ftfiou(unit, status)
  
  !  Check for any error, and if so print out error messages.
  !  The PRINTERROR subroutine is listed near the end of this file.
  if (status .gt. 0)call printerror(status)
  
end subroutine read_fits_1d

!============================================================
subroutine write_fits_1d( filename, n, array)
!! cookbook.f in the Cfitsio distribution
!! see also http://heasarc.gsfc.nasa.gov/docs/software/fitsio/c/f_user/node8.html
  use healpix_types
  implicit none

!     Create a FITS primary array containing a 2-D image

  character(len=*), intent(in) :: filename
  integer(I4B),     intent(in) :: n
  real(DP), dimension(0:n-1), intent(in) :: array

  character(len=200) :: file_out
  integer(I4B) :: status,unit,blocksize,bitpix,naxis,naxes(2)
  integer(I4B) :: i,j,group,fpixel,nelements
  logical :: simple, extend

  status=0
  !     Name of the FITS file to be created:
  !filename='ATESTFILE.FITS'

  file_out = "!"//trim(filename)

  !     Get an unused Logical Unit Number to use to create the FITS file
  call ftgiou(unit,status)

  !     create the new empty FITS file
  blocksize=1
  call ftinit(unit,file_out,blocksize,status)

  !     initialize parameters about the FITS image (300 x 200 16-bit integers)
  !bitpix=16
  
  !32-bit double
  bitpix = -64
  simple=.true.
  naxis=1
  naxes(1) = n
  extend=.true.

  !     write the required header keywords
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  
  !     write the array to the FITS file
  group=1
  fpixel=1
  nelements=naxes(1)
  !call ftpprj(unit,group,fpixel,nelements,map,status)
  call ftpprd(unit,group,fpixel,nelements,array,status)
  
  !     close the file and free the unit number
  call ftclos(unit, status)
  call ftfiou(unit, status)

end subroutine write_fits_1d

!======================================================
subroutine printerror(status)

!  This subroutine prints out the descriptive text corresponding to the
!  error status value and prints out the contents of the internal
!  error message stack generated by FITSIO whenever an error occurs.
  
  integer status
  character errtext*30,errmessage*80
  
  !  Check if status is OK (no error); if so, simply return
  if (status .le. 0)return
  
  !  The FTGERR subroutine returns a descriptive 30-character text string that
  !  corresponds to the integer error status number.  A complete list of all
  !  the error numbers can be found in the back of the FITSIO User's Guide.
  call ftgerr(status,errtext)
  print *,'FITSIO Error Status =',status,': ',errtext
  
  !  FITSIO usually generates an internal stack of error messages whenever
  !  an error occurs.  These messages provide much more information on the
  !  cause of the problem than can be provided by the single integer error
  !  status value.  The FTGMSG subroutine retrieves the oldest message from
  !  the stack and shifts any remaining messages on the stack down one
  !  position.  FTGMSG is called repeatedly until a blank message is
  !  returned, which indicates that the stack is empty.  Each error message
  !  may be up to 80 characters in length.  Another subroutine, called
  !  FTCMSG, is available to simply clear the whole error message stack in
  !  cases where one is not interested in the contents.
  call ftgmsg(errmessage)
  do while (errmessage .ne. ' ')
     print *,errmessage
     call ftgmsg(errmessage)
  end do
end subroutine printerror

end module poker_tools
