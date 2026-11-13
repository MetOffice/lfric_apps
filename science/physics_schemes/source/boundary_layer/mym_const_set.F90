! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Module mym_const_set----------------------------------------------

!  Purpose: To set constants defined in mym_const_mod

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
module mym_const_set_mod

use um_types, only: r_bl

implicit none

character(len=*), parameter, private :: ModuleName = 'MYM_CONST_SET_MOD'
contains

subroutine mym_const_set

use mym_const_mod, only: g1,b1,b2,c2,c3,c4,c5,pr,a1,c1,a2,g2,a1_2,             &
    rfc,f1,f2,rf1,rf2,smc,shc,ri1,ri2,ri3,ri4,cc2,cc3,e1c,e2c,e3c,             &
    e4c,e5c,my_alpha1,my_alpha2,my_alpha3,my_alpha4,elt_min,                   &
    one_third,two_thirds,coef_trbvar_diff_tke,coef_trbvar_diff,                &
    qke_max,e_trb_max
use mym_option_mod, only: l_my3_improved_closure
use parkind1, only: jprb, jpim
use yomhook, only: lhook, dr_hook
implicit none

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='MYM_CONST_SET'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
if (l_my3_improved_closure) then
  ! closure constants in the improved MY model (MYNN)
  g1  =  0.235
  b1  = 24.0
  b2  = 15.0
  c2  =  0.7
  c3  =  0.323
  c4  =  0.0
  c5  =  0.2
  pr  =  0.74
else
  ! closure constants in the original MY model
  g1 = 0.222
  b1 = 16.6
  b2 = 10.1
  c2 = 0.0
  c3 = 0.0
  c4 = 0.0
  c5 = 0.0
  pr  =  0.80
end if
! Combined constants
a1  = b1 * ( 1.0 - 3.0 * g1 ) / 6.0
c1  = g1 - 1.0 / ( 3.0 * a1 * b1 ** (1.0 / 3.0))
a2  = a1 *( g1 - c1 ) / ( g1 * pr )
g2  = b2 / b1 * ( 1.0 - c3 ) + 2.0 * a1 / b1 * ( 3.0 - 2.0 * c2 )
a1_2 = a1 / a2

rfc = g1 / (g1 + g2)
f1 = b1 * (g1 - c1)  + 3.0 * a2 *(1.0 - c2 ) * ( 1.0 - c5 )                    &
                       +2.0 * a1 *(3.0 -2.0 * c2)
f2  = b1 * (g1 + g2) - 3.0 * a1 *(1.0 - c2)
rf1 = b1 * (g1 - c1) / f1
rf2 = b1 * g1 / f2
smc = a1 / a2 * f1 / f2
shc = 3.0 * a2 * (g1 + g2)

ri1 = 0.5 /smc
ri2 = rf1 * smc
ri3 = 4.0 * rf2 * smc - 2.0 * ri2
ri4 = ri2 ** 2

cc2 =  1.0 - c2
cc3 =  1.0 - c3
e1c =  3.0 * a2 * b2 * cc3
e2c =  9.0 * a1 * a2 * cc2
e3c =  9.0 * a2 * a2 * cc2 * (1.0 - c5)
e4c = 12.0 * a1 * a2 *cc2
e5c =  6.0 * a1 *a1

! Other parameters
my_alpha1 = 0.23
my_alpha2 = 1.0
my_alpha3 = 5.0
my_alpha4 = 100.0
elt_min = 20.0

one_third = 1.0 / 3.0
two_thirds = 2.0 / 3.0

coef_trbvar_diff_tke = 3.0
coef_trbvar_diff = 1.0

qke_max = 500.0
e_trb_max = 0.5 * qke_max

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine mym_const_set
end module mym_const_set_mod
