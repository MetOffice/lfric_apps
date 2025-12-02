! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module calc_turb_perts_mod

implicit none

contains


!----------------------------------------------------------------
! Subroutine to calculate the parcel initial perturbations
! for winds and scalars, based on the turbulent fluxes on
! their native half-level grid.
!----------------------------------------------------------------
subroutine calc_turb_perts( n_points, n_points_super,                          &
                            turb_super, turb_pert_fields )

use comorph_constants_mod, only: real_cvprec, zero, one, sqrt_min_float,       &
                                 par_gen_pert_fac, par_gen_w_fac
use turb_type_mod, only: n_turb, i_w_var,                                      &
                         i_f_templ, i_f_q_tot,                                 &
                         i_f_wind_u, i_f_wind_v
use fields_type_mod, only: i_wind_u, i_wind_v, i_wind_w,                       &
                           i_q_vap, i_temperature

implicit none

! Number of points
integer, intent(in) :: n_points

! Number of points in super-arrays (might be > n_points due to
! re-use of same arrays for differing-length compression lists).
integer, intent(in) :: n_points_super

! Super-array containing turbulence fields
real(kind=real_cvprec), intent(in) :: turb_super                               &
                                      ( n_points_super, n_turb )

! Super-array containing turbulence-based perturbations to each
! primary field
real(kind=real_cvprec), intent(out) :: turb_pert_fields                        &
                                       ( n_points, i_wind_u:i_q_vap )

! w-scale from TKE
real(kind=real_cvprec) :: work_w(n_points)
! Inverse w-scale for converting fluxes into perturbations
real(kind=real_cvprec) :: recip_w(n_points)

! Loop counter
integer :: ic


! Compute w scaling for converting fluxes <w'phi'> into
! parcel perturbations
do ic = 1, n_points
  work_w(ic) = sqrt( turb_super(ic,i_w_var) )
end do

! Check for w too close to zero for safe division
do ic = 1, n_points
  if ( work_w(ic) < sqrt_min_float ) then
    ! Set values to make all outputs zero if w is too near zero
    work_w(ic)  = zero
    recip_w(ic) = zero
  else
    ! Safe to divide
    recip_w(ic) = one / work_w(ic)
  end if
end do

! Set w perturbation
do ic = 1, n_points
  turb_pert_fields(ic,i_wind_w) = par_gen_w_fac * work_w(ic)
end do

! Scale by factor for calculating phi' from <w'phi'>
do ic = 1, n_points
  recip_w(ic) = par_gen_pert_fac * recip_w(ic)
end do

! Add perturbations to u, v, T, q
do ic = 1, n_points
  turb_pert_fields(ic,i_wind_u) = recip_w(ic) * turb_super(ic,i_f_wind_u)
end do
do ic = 1, n_points
  turb_pert_fields(ic,i_wind_v) = recip_w(ic) * turb_super(ic,i_f_wind_v)
end do
do ic = 1, n_points
  turb_pert_fields(ic,i_temperature) = recip_w(ic) * turb_super(ic,i_f_templ)
end do
do ic = 1, n_points
  turb_pert_fields(ic,i_q_vap) = recip_w(ic) * turb_super(ic,i_f_q_tot)
end do


return
end subroutine calc_turb_perts


end module calc_turb_perts_mod
