! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Purpose: To calculate fundamental values such as non-dimensional
!           diffusion coefficients in the MY model with level 2
!           scheme.
!  Programming standard : UMDP 3

!  Documentation: UMDP 025

!  This code is based on the code provided by the authors who wrote
!  the following papers.
!   * Nakanishi, M. and H. Niino, 2009: Development of an improved
!      turbulence closure model for the atmospheric boundary layer.
!      J. Meteor. Soc. Japan, 87, 895-912.
!   * Nakanishi, M. and H. Niino, 2006: An improved Mellor-Yamada
!      Level-3 model: Its numerical stability and application to
!      a regional prediction of advection fog.
!      Boundary-Layer Meteor., 119, 397-407.
!   * Nakanishi, M. and H. Niino, 2004: An improved Mellor-Yamada
!      Level-3 model with condensation physics: Its design and
!      verification.
!      Boundary-Layer Meteor., 112, 1-31.
!   * Nakanishi, M., 2001: Improvement of the Mellor-Yamada
!      turbulence closure model based on large-eddy simulation data.
!      Boundary-Layer Meteor., 99, 349-378.
!   The web site publicising their code:
!    http://www.nda.ac.jp/~naka/MYNN/index.html

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
module mym_level2_mod

use um_types, only: r_bl

implicit none

character(len=*), parameter, private :: ModuleName = 'MYM_LEVEL2_MOD'
contains

subroutine mym_level2(                                                         &
      bl_levels, dbdz, dvdzm, gm, gh, sm, sh)

use atm_fields_bounds_mod, only: tdims
use mym_const_mod, only: ri1, ri2, ri3, ri4, rfc, rf1, rf2, shc, smc
use mym_option_mod, only: tke_levels
use parkind1, only: jprb, jpim
use yomhook, only: lhook, dr_hook
implicit none

! Intent IN Variables
integer, intent(in) ::                                                         &
   bl_levels
                  ! Max. no. of "boundary" levels

real(kind=r_bl), intent(in) ::                                          &
   dbdz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        2:tke_levels),                                                         &
                  ! Buoyancy gradient across layer
                  ! interface interpolated to theta levels.
                  ! (:,:,K) repserents the value on theta level K-1
   dvdzm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
         2:bl_levels)
                  ! Modulus of wind shear at theta levels.
                  ! (:,:,K) repserents the value on theta level K-1

! Intent OUT Variables
real(kind=r_bl), intent(out) ::                                         &
   gm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),         &
                 ! square of wind shear on theta level K-1
                 ! (a denominator of gradient Richardson number)
   gh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),         &
                 ! - buoyancy gradient on theta level K-1
                 ! (a numerator of gradient Richardson number)
   sm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),         &
                 ! Non-dimensional diffusion coefficients for
                 ! momentum from level 2 scheme
                 ! defined on theta level K-1
   sh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels)
                 ! Non-dimensional diffusion coefficients for
                 ! scalars from level 2 scheme
                 ! define on theta level K-1
! Local variables
integer ::                                                                     &
   i, j, k
                 ! Loop indexes

real(kind=r_bl) ::                                                      &
   ri,                                                                         &
                 ! gradient Richardson Number
   rf
                 ! flux Richardson Number

real(kind=r_bl), parameter ::                                           &
   ri_max = 1.0e5
                 ! upper limit of gradient Richardson number to avoid
                 ! loss of significance in single precision
                 ! When Ri exceeds 10^6, stability function gets non-zero value
                 ! even though Ri exceeds the critical Richardson number
                 ! in single precision.

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='MYM_LEVEL2'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do k = 2, tke_levels
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      gm(i, j, k) =   dvdzm(i, j, k) * dvdzm(i, j, k)
      gh(i, j, k) = - dbdz(i, j, k)
      !   Gradient Richardson number
      ri = - gh(i, j, k) / max( gm(i, j, k), 1.0e-10 )
      ri = min( ri, ri_max )
      !   Flux Richardson number
      rf = min(ri1 * (ri + ri2 - sqrt(ri ** 2 - ri3 * ri + ri4)),              &
               rfc )
      sh(i, j, k) = shc * (rfc - rf) / (1.0 - rf)
      sm(i, j, k) = smc * (rf1 - rf) / (rf2 - rf) * sh(i, j, k)
    end do
  end do
end do

do j = tdims%j_start, tdims%j_end
  do i = tdims%i_start, tdims%i_end
    gm(i, j, 1) = 0.0
    gh(i, j, 1) = 0.0
    sh(i, j, 1) = 0.0
    sm(i, j, 1) = 0.0
  end do
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine mym_level2
end module mym_level2_mod
