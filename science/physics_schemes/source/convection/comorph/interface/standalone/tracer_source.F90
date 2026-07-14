! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

#if !defined(LFRIC)
module tracer_source_mod

implicit none

contains

! Subroutine adds on in-parcel source or sink terms for tracer species
! (increments to in-parcel values over the current level-step).
!
! The rest of comorph treats tracers as generic passively-transported
! scalars; any calculations that are specific to the actual variables
! represented by the tracers in the stand-alone test should be done in here
!
subroutine tracer_source( n_points, n_points_super,                            &
                          massflux_d, dq_prec, par_fields, par_tracers )

use comorph_constants_mod, only: real_cvprec, n_tracers
use fields_type_mod, only: n_fields

implicit none

! Number of points
integer, intent(in) :: n_points
! Array size for super-arrays (maybe larger than n_points where array
! is reused to save on having to deallocate / reallocate to smaller size)
integer, intent(in) :: n_points_super

! Parcel dry-mass flux / kg m-2
real(kind=real_cvprec), intent(in) :: massflux_d(n_points)

! In-parcel precip mixing-ratio production increment over the current
! level-step (summed over all precip species) / kg kg-1
real(kind=real_cvprec), intent(in) :: dq_prec(n_points)

! In-parcel primary fields (latest values)
real(kind=real_cvprec), intent(in) :: par_fields                               &
                                      ( n_points_super, n_fields )

! In-parcel tracer concentrations to be updated
real(kind=real_cvprec), intent(in out) :: par_tracers                          &
                                          ( n_points_super, n_tracers )


! So far this is just a dummy routine that does nothing, since no
! source or sink terms are yet applied to tracers in the standalone test.


return
end subroutine tracer_source

end module tracer_source_mod
#endif
