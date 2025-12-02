! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module lfricinp_surface_parameters_mod

use constants_mod, only: r_def, i_def

implicit none

private

public :: sm_levels, dzsoil

! Number of soil moisture levels - this corresponds to jules_physics_init
integer(kind=i_def), parameter :: sm_levels = 4

! Soil moisture layer thickness - this corresponds to jules_physics_init
real(kind=r_def), parameter    :: dzsoil(1:sm_levels) =                        &
                            (/ 0.1_r_def, 0.25_r_def, 0.65_r_def, 2.0_r_def /)

end module lfricinp_surface_parameters_mod
