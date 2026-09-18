! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module microphysics_2_mod

implicit none

contains

! Subroutine to do final calculations for BOTEMS
! (Back-Of-The-Envelope Microphysical Scheme)
!
! This includes any microphysical processes that need to be done
! after the implicit solution of phase-changes.
! Currently only does autoconversion of liquid-cloud to rain
subroutine microphysics_2( n_points, n_points_super, nc, index_ic,             &
                           delta_t, vert_len, wf_cond, q_cond,                 &
                           l_diags, moist_proc_diags,                          &
                           n_points_diag, n_diags, diags_super )

use comorph_constants_mod, only: real_cvprec, n_cond_species, i_cond_cl,       &
                                 i_cond_rain, l_cv_rain
use moist_proc_diags_type_mod, only: moist_proc_diags_type
use autoconversion_mod, only: autoconversion

implicit none

! Number of points
integer, intent(in) :: n_points
! Number of points in the condensed water species super-array
! (this might be reused on subsequent calls with bigger size
!  than needed, to save having to re-allocate it)
integer, intent(in) :: n_points_super

! Number of points where each hydrometeor species is non-zero
integer, intent(in) :: nc(n_cond_species)
! Indices of those points
integer, intent(in) :: index_ic(n_points,n_cond_species)

! Time interval for converting process rates to increments.
real(kind=real_cvprec), intent(in) :: delta_t(n_points)

! Vertical length-scale of the parcel.
real(kind=real_cvprec), intent(in) :: vert_len(n_points)

! Fall-speed of each hydrometeor species
real(kind=real_cvprec), intent(in) :: wf_cond                                  &
                                   ( n_points, n_cond_species )

! Mixing ratios of condensed water species
real(kind=real_cvprec), intent(in out) :: q_cond                               &
                             ( n_points_super, n_cond_species )

! Master switch for whether or not to calculate any diagnostics
logical, intent(in) :: l_diags
! Structure containing diagnostic flags etc.
type(moist_proc_diags_type), intent(in) :: moist_proc_diags
! Number of points in the diagnostics super-array
integer, intent(in) :: n_points_diag
! Total number of diagnostics in the super-array
integer, intent(in) :: n_diags
! Super-array to store all the diagnostics
real(kind=real_cvprec), intent(in out) :: diags_super                          &
                                         ( n_points_diag, n_diags )


if ( nc(i_cond_cl) > 0 .and. l_cv_rain) then
  ! If any liquid cloud present
  ! Autoconversion of liquid cloud to rain
  call autoconversion( n_points, nc(i_cond_cl), index_ic(:,i_cond_cl),         &
                       delta_t, vert_len, wf_cond(:,i_cond_cl),                &
                       q_cond(:,i_cond_cl), q_cond(:,i_cond_rain),             &
                       l_diags, moist_proc_diags % diags_cl,                   &
                                moist_proc_diags % diags_rain,                 &
                       n_points_diag, n_diags, diags_super )
end if


return
end subroutine microphysics_2

end module microphysics_2_mod
