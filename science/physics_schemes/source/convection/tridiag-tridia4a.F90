! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Solves the equation A.X=Y where A is a tridiagonal matrix
!

module tridiag_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'TRIDIAG_MOD'
contains

subroutine tridiag(a,b,c,r,u,n)

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

!------------------------------------------------------------------------
! Description:
!   Solves the equation A.X=Y where A is a tridiagonal matrix
!
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.3 programming standards.
!-----------------------------------------------------------------------
! Subroutine arguments

integer, intent(in) ::                                                         &
  n                      ! size of vectors


real(kind=real_umphys), intent(in)    ::                                       &
  a(n)                 & ! Components of tridiagonal matrix
 ,b(n)                 & !
 ,c(n)                 & !
 ,r(n)                   ! RHS of linear equation

real(kind=real_umphys), intent(out) ::                                         &
  u(n)                   ! Solution vector

! Local variables

integer ::                                                                     &
  j                ! Loop counter

real(kind=real_umphys) ::                                                      &
  gam(n)         & ! Dimension of gam, should be the same as n
 ,bet

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='TRIDIAG'
!------------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!------------------------------------------------------------------------

bet=b(1)
u(1)=r(1)/bet
do j=2,n
  gam(j)=c(j-1)/bet
  bet=b(j)-a(j)*gam(j)
  u(j)=(r(j)-a(j)*u(j-1))/bet
end do
do j=n-1,1,-1
  u(j)=u(j)-gam(j+1)*u(j+1)
end do
!------------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

!------------------------------------------------------------------------

return
end subroutine tridiag
end module tridiag_mod
