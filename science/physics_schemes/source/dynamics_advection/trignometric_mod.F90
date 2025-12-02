! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Constants necessary for Coriolis terms in advection scheme.

module trignometric_Mod

! Description:
! This module is used to hold trignometric values set in atmos_init.
!
! Method:
! The trignometric arrays are calculated in routine control/top_level/
! set_trig called from atmos_init,
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

! Trignometric Arrays
! Most names are self-defining

! --- Variables relating to grid latitude or grid longitude ---
! Values may differ from true latitude and true longitude (e.g. in
! the case of a limited area model (LAM) where a rotated grid is used).
! Variables that are targets, are so as they are targets in
! multivariate swap_bounds.
real, allocatable, target ::  cos_theta_latitude (:,:)


end module trignometric_Mod
