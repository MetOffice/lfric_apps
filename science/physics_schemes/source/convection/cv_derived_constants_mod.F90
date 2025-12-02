! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module holding commonly derived constants used in convection
!
module cv_derived_constants_mod

use planet_constants_mod, only:                                                &
  ls,                           & ! Latent heat of sublimation
  lsrcp,                                                                       &
  lcrcp,                                                                       &
  lfrcp,                                                                       &
  gamma_dry => grcp,            & ! dry adiabatic lapse rate
  cv,                           & ! specific heat of dry air at constant volume
  ra2 => recip_a2

implicit none

!-------------------------------------------------------------------
! Description:
! This module calculate commonly used constants derived from
! UM constants
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.3 programming standards.
!-------------------------------------------------------------------

end module cv_derived_constants_mod
