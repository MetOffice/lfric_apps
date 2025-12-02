! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module set_qsat_mod

! Routines for calculating saturation water-vapour mixing-ratio qsat...

! If running CoMorph inside the UM, the host model already has
! its own qsat functions, so this module just contains
! wrapper routines around those

! Note: important to point at the mixing-ratio versions of the
! UM qsat routines!

implicit none

contains


! Wrapper around UM routine for qsat w.r.t. liquid water
subroutine set_qsat_liq( n_points, temperature, pressure, qsat )

use comorph_constants_mod, only: real_cvprec
use qsat_mod, only: qsat_wat_mix

implicit none

! Number of points
integer, intent(in) :: n_points

! Input temperature and pressure
real(kind=real_cvprec), intent(in) :: temperature(n_points)
real(kind=real_cvprec), intent(in) :: pressure(n_points)

! Output saturation mixing-ratio
real(kind=real_cvprec), intent(out) :: qsat(n_points)

! New version of qsat actually points to module procedure interface
! statements in the UM's qsat_mod, which should automatically
! select the version of the code with the correct precision
! (real32 or real64).  So no need to convert to host-model precision.
call qsat_wat_mix( qsat, temperature, pressure, n_points )

return
end subroutine set_qsat_liq


! Wrapper around UM routine for qsat w.r.t. ice
! (the UM routines actually return qsat w.r.t. liquid water when
! above the melting point, and w.r.t. ice when below; but comorph
! only uses this when below the melting point).
subroutine set_qsat_ice( n_points, temperature, pressure, qsat )

use comorph_constants_mod, only: real_cvprec
use qsat_mod, only: qsat_mix

implicit none

! Number of points
integer, intent(in) :: n_points

! Input temperature and pressure
real(kind=real_cvprec), intent(in) :: temperature(n_points)
real(kind=real_cvprec), intent(in) :: pressure(n_points)

! Output saturation mixing-ratio
real(kind=real_cvprec), intent(out) :: qsat(n_points)

! New version of qsat actually points to module procedure interface
! statements in the UM's qsat_mod, which should automatically
! select the version of the code with the correct precision
! (real32 or real64).  So no need to convert to host-model precision.
call qsat_mix( qsat, temperature, pressure, n_points )

return
end subroutine set_qsat_ice


end module set_qsat_mod
