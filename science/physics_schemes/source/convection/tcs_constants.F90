! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Module holding collection of physical constants and derived
! combinations thereof.
!
module tcs_constants

use planet_constants_mod, only:                                                &
    r, cp, repsilon, c_virtual, kappa, pref, rv,                               &
    recip_epsilon, recip_kappa, g, cv,                                         &
    lc_o_cp => lcrcp, gamma_dry => grcp, ra2 => recip_a2

use water_constants_mod, only: lc, lf

use cv_diag_param_mod, only:                                                   &
   a_bolton, b_bolton, c_bolton, d_bolton

implicit none
!
! Description:
!   This module holds a collection of physical constants and derived
! combinations thereof.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.1 programming standards.
!

end module tcs_constants
