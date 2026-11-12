! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
! Description:
!   Module containing general physics runtime options
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: top_level
!
! Code Description:
!   Language: FORTRAN 95

module gen_phys_inputs_mod

implicit none

! Use mixing ratio in atmos_physics1, atmos_physics2 and end of atm_step_4A
logical :: l_mr_physics = .false.

! Switch to use volume-based interpolation of rho to theta-levels is true
! Note that this is only valid when the rho levels lie halfway between the
! theta levels
logical :: l_vol_interp_rho

end module gen_phys_inputs_mod
