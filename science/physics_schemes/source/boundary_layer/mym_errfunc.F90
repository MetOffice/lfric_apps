! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Purpose: To fastly calculate values of error function in the MY
!           model.
!           The calculation is based on the expansion up to
!           the 13th order.

!  Programming standard : UMDP 3

!  Documentation: UMDP 025

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
module mym_errfunc_mod

use um_types, only: r_bl

implicit none

character(len=*), parameter, private :: ModuleName = 'MYM_ERRFUNC_MOD'
contains

subroutine mym_errfunc(nn, x, y)

use conversions_mod, only: pi
use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
implicit none

integer, intent(in) :: nn       ! size of array

real(kind=r_bl), intent(in)    :: x(nn)    ! input array

real(kind=r_bl), intent(out)   :: y(nn)    ! output array

! Local Variables
integer             :: i        ! Loop index

real(kind=r_bl) ::                                                             &
   x01,                                                                        &
       ! x with upper limit
   x02,                                                                        &
       ! x powered by 2
   x04,                                                                        &
       ! x powered by 4
   x06,                                                                        &
       ! x powered by 6
   x08,                                                                        &
       ! x powered by 8
   x10,                                                                        &
       ! x powered by 10
   x12
       ! x powered by 12

real(kind=r_bl), parameter ::                                           &
   erfmax = 1.0
       ! upper limit of the value to avoid it outside domain

real(kind=r_bl), parameter ::                                           &
   argmax = 100.0
       ! upper limit of the arguments to avoid floating overflow
       ! Given the precision of this Taylor expansion, calculations for
       ! |x|>1.65 are sticked to erfmax and yield no meaningful results,
       ! so this poses no problem.

real(kind=r_bl), parameter :: factor =  2.0 / sqrt(pi)
real(kind=r_bl), parameter :: c01 = factor * 1.0
real(kind=r_bl), parameter :: c03 = factor * 1.0 /    3.0
real(kind=r_bl), parameter :: c05 = factor * 1.0 /   10.0
real(kind=r_bl), parameter :: c07 = factor * 1.0 /   42.0
real(kind=r_bl), parameter :: c09 = factor * 1.0 /  216.0
real(kind=r_bl), parameter :: c11 = factor * 1.0 / 1320.0
real(kind=r_bl), parameter :: c13 = factor * 1.0 / 9360.0

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='MYM_ERRFUNC'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do i = 1, nn
  x01 = max(min(x(i), argmax), -argmax)
  x02 = x01 * x01
  x04 = x02 * x02
  x06 = x04 * x02
  x08 = x06 * x02
  x10 = x08 * x02
  x12 = x10 * x02
  y(i) = x01 * (                                                               &
        + c01                                                                  &
        - c03 * x02                                                            &
        + c05 * x04                                                            &
        - c07 * x06                                                            &
        + c09 * x08                                                            &
        - c11 * x10                                                            &
        + c13 * x12)
  if (x01 > 0) then
    y(i) = min(y(i), erfmax)
  else
    y(i) = max(y(i), -erfmax)
  end if
end do
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine mym_errfunc
end module mym_errfunc_mod
