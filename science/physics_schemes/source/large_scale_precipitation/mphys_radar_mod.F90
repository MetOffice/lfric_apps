! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module mphys_radar_mod

! Description:
! Holds reflectivity constants required by the microphysics
! scheme

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

use mphys_constants_mod, only: mprog_min

use um_types, only: real_umphys

implicit none

real(kind=real_umphys), parameter :: kliq  = 0.93
real(kind=real_umphys), parameter :: kice  = 0.174
real(kind=real_umphys), parameter :: mm6m3 = 1.0e18

! Define reflectivity limit
real(kind=real_umphys), parameter :: ref_lim = -35.0 ! dBZ
! Convert this to linear units (mm6 m-3) using 10.0**p
! Where p = -35 dBZ / 10.0 .
real(kind=real_umphys), parameter :: ref_lim_lin = 3.1623e-4

! Mixing ratio limit (below which we ignore the species)
! Set this to be the same as the absolute value used in the
! rest of the microphysics.
real(kind=real_umphys), parameter :: mr_lim = mprog_min

! Cloud fraction limit (below which we ignore to avoid massive
! reflectivity values)- currently set as 1% of the grid box
real(kind=real_umphys), parameter :: cf_lim = 0.01

! Cloud drop number concentration limit (below which we ignore liquid cloud)
! Equivalent to 5 per cm3, used by most of the UM routines
real(kind=real_umphys), parameter :: nd_lim  = 5.0e6

real(kind=real_umphys), parameter :: rho_g   = 500.0
! (This was incorrectly set to 900.0 in the original HWT code.)
real(kind=real_umphys), parameter :: rho_i   = 900.0
real(kind=real_umphys), parameter :: rho_i2  = 900.0

real(kind=real_umphys), parameter :: ref_mom = 4.0 ! Radar reflectivity moment


end module mphys_radar_mod
