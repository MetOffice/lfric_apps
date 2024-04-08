! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lfricinp_physics_constants_mod

USE constants_mod, ONLY: r_def

IMPLICIT NONE

PRIVATE

PUBLIC :: density_h2o

! Density of water in kg/m3
REAL(KIND=r_def), PARAMETER :: density_h2o = 1000.0_r_def

END MODULE lfricinp_physics_constants_mod
