! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!   solves tridiagonal matrix -  convection scheme

module tridiag_all_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'TRIDIAG_ALL_MOD'
contains

subroutine tridiag_all(n,nvec,nvec_len,a,b,c,r,u)

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none
!------------------------------------------------------------------------
! Description:
!  Solves the equations A.X = Y,  where A is a tridiagonal matrix
!          for several matrices A at once.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.3 programming standards.
!-----------------------------------------------------------------------
! Subroutine arguments

integer, intent(in) ::                                                         &
  n                    & ! maximum size of vectors X and Y
 ,nvec                 & ! Number of vectors/matrices to solve
 ,nvec_len(nvec)         ! length of each vector.

real(kind=real_umphys), intent(in) ::                                          &
  a(nvec,n)         & ! Components of tridiagonal matrix
 ,b(nvec,n)         & ! Components of tridiagonal matrix
 ,c(nvec,n)         & ! Components of tridiagonal matrix
 ,r(nvec,n)           ! vector Y, R.H.S of linear equation

real(kind=real_umphys), intent(out) ::                                         &
  u(nvec,n)            ! solution vectors


!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

integer ::                                                                     &
  i,j            ! loop counter

real(kind=real_umphys) ::                                                      &
  gam(nvec,n)  & ! work array
 ,bet(nvec)

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='TRIDIAG_ALL'
!-----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------
do i=1,nvec
  bet(i) = b(i,1)
  u(i,1) = r(i,1)/bet(i)
end do

do j=2,n
  do i=1, nvec
    if (j <= nvec_len(i)) then
      gam(i,j) = c(i,j-1)/bet(i)
      bet(i)   = b(i,j) - a(i,j)*gam(i,j)
      !for info !    if (bet(i) == 0.0) stop   ! was in original code
      u(i,j)   = (r(i,j) - a(i,j)*u(i,j-1))/bet(i)
    end if
  end do
end do

do j=n-1,1,-1
  do i=1, nvec
    if (j <= (nvec_len(i)-1)) then
      u(i,j) = u(i,j) - gam(i,j+1)*u(i,j+1)
    end if
  end do
end do

!-----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

!-----------------------------------------------------------------------
return
end subroutine tridiag_all
end module tridiag_all_mod
