! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module lfricinp_physics_constants_mod

use constants_mod, only: r_def

implicit none

private

public :: density_h2o

! Density of water in kg/m3
real(kind=r_def), parameter :: density_h2o = 1000.0_r_def

end module lfricinp_physics_constants_mod
