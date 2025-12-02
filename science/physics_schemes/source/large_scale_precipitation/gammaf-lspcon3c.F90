! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Calculates constants used in large-scale precipitation scheme.
!  subroutine GAMMAF--------------------------------------------------
!   PURPOSE: CALCULATES COMPLETE GAMMAF function BY
!   A POLYNOMIAL APPROXIMATION
!
!   Notes: The original routine was only designed for arguments >1.
!          From VN7.6 this has modified this to work with arguments <1
!          following guidance from Ben Shipway.
!
!          NB For negative integers the solution is undefined and
!          a divide by zero will occur, but model settings should prevent
!          this happening.
! --------------------------------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation
module gammaf_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName='GAMMAF_MOD'

contains

subroutine gammaf(y,gam)

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
implicit none
real(kind=real_umphys) ::                                                      &
                            !, intent(in)
  y
real(kind=real_umphys) ::                                                      &
                            !, intent(out)
  gam
! Gamma function of Y

! LOCAL VARIABLE
integer :: i,m
real(kind=real_umphys) :: gg,g,pare,x

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='GAMMAF'

! --------------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
gg=1.0
m=floor(y)
x=y-m
if (m > 1) then
  do i = 1, m-1
    g=y-i
    gg=gg*g
  end do
else if (m < 1) then
  do i = m, 0
    g=y-i
    gg=gg/g
  end do
end if
pare=-0.5748646*x+0.9512363*x*x-0.6998588*x*x*x                                &
+0.4245549*x*x*x*x-0.1010678*x*x*x*x*x+1.0
gam=pare*gg
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine gammaf

! --------------------------------------------------------------------
! --------------------------------------------------------------------

subroutine gammaf_lower_incomp(s_in,x_in,gamma_out,nn)

use yomhook,    only: lhook, dr_hook
use parkind1,   only: jprb, jpim

implicit none

! Input parameters
real(kind=real_umphys), intent(in)    :: s_in,x_in
! Output value of lower incomplete gamma function
real(kind=real_umphys), intent(out)   :: gamma_out
! Number of elements in summation.
integer, intent(in) :: nn

! Local variables
integer :: k
real(kind=real_umphys) :: factor, BigGamma, tmpsum

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='GAMMAF_LOWER_INCOMP'

! --------------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

call gammaf(s_in, BigGamma)

factor = (x_in**s_in) * BigGamma * exp(-x_in)

tmpsum = 0.0

do k = 0, nn
  call gammaf(s_in+real(k, kind=real_umphys)+1.0_real_umphys, BigGamma)

  tmpsum = tmpsum + ( (x_in**k) / BigGamma )

end do

gamma_out = factor * tmpsum

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

return

end subroutine gammaf_lower_incomp

end module gammaf_mod
