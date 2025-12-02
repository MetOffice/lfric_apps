! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Constants necessary for model level heights in advection and other schemes.

module level_heights_Mod

! Description:
! This module is used to hold levels values set in atmos_init.
!
! Method:
! The height levels arrays are calculated in routine control/top_level/
! set_levels called from atmos_init,
! and used in dynamics_advection and other routines.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: dynamics_advection
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

implicit none

! Arguments

! Heights Arrays
! - eta values of theta levels
real, allocatable, target  ::  eta_theta_levels(:)
! - height of theta levels
real, allocatable, target  ::  r_theta_levels (:,:,:)
! - height of rho levels
real, allocatable, target  ::  r_rho_levels   (:,:,:)

end module level_heights_Mod
