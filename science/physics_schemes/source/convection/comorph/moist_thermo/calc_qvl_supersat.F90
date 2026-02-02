! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module calc_qvl_supersat_mod

implicit none

contains

! Subroutine to calculate total-water supersaturation
! supersat = q_vap + q_cl - qsat_liq(Tl)
subroutine calc_qvl_supersat( n_points, n_points_super,                        &
                              pressure, temperature, q_vap, q_cond,            &
                              qvl_supersat, virt_temp_noliq, linear_qs )

use comorph_constants_mod, only: real_cvprec, zero, one,                       &
                                 n_cond_species, i_cond_cl
use linear_qs_mod, only: n_linear_qs_fields, i_ref_temp, i_qsat_liq_ref,       &
                         i_dqsatdT_liq
use set_cp_tot_mod, only: set_cp_tot
use lat_heat_mod, only: lat_heat_incr, i_phase_change_evp
use set_qsat_mod, only: set_qsat_liq
use calc_virt_temp_mod, only: calc_virt_temp

implicit none

! Number of points
integer, intent(in) :: n_points
! Dimensions of condensate super-array
integer, intent(in) :: n_points_super

! Air pressure
real(kind=real_cvprec), intent(in) :: pressure(n_points)

! Temperature and moisture mixing-ratios
real(kind=real_cvprec), intent(in) :: temperature(n_points)
real(kind=real_cvprec), intent(in) :: q_vap(n_points)
real(kind=real_cvprec), intent(in) :: q_cond(n_points_super,n_cond_species)

! Output super-saturation that the air would have if all liquid-cloud
! were evaporated
real(kind=real_cvprec), intent(out) :: qvl_supersat(n_points)

! Optionally also output the virtual temperature the air would have
real(kind=real_cvprec), intent(out) :: virt_temp_noliq(n_points)

! Super-array for computing qsat using linearisation about an existing value
real(kind=real_cvprec), optional, intent(in) :: linear_qs                      &
                                              ( n_points, n_linear_qs_fields )

! Liquid-water temperature
real(kind=real_cvprec) :: temperature_l(n_points)
! q_vap + q_cl
real(kind=real_cvprec) :: q_vap_l(n_points)
! Condensed water super-array
real(kind=real_cvprec) :: q_cond_l ( n_points, n_cond_species )

! Total heat capacity
real(kind=real_cvprec) :: cp_tot(n_points)
! Saturation vapour mixing-ratio
real(kind=real_cvprec) :: qsat(n_points)

! Loop counter
integer :: ic, i_cond


! Set total heat capacity of the air
call set_cp_tot( n_points, n_points_super, q_vap, q_cond, cp_tot )

! Increment the temperature with latent heating from evaporating q_cl
! (i.e. calculate liquid-water temperature Tl)
do ic = 1, n_points
  temperature_l(ic) = temperature(ic)
end do
call lat_heat_incr( n_points, n_points, i_phase_change_evp,                    &
                    cp_tot, temperature_l, dq=q_cond(:,i_cond_cl) )

! Increment q_vap and q_cl consistently
do i_cond = 1, n_cond_species
  do ic = 1, n_points
    q_cond_l(ic,i_cond) = q_cond(ic,i_cond)
  end do
end do
do ic = 1, n_points
  q_vap_l(ic) = q_vap(ic) + q_cond(ic,i_cond_cl)
  q_cond_l(ic,i_cond_cl) = zero
end do

! Calculate saturation vapour mixing ratio w.r.t. liquid water
! at the liquid-water temperature
if ( present(linear_qs) ) then
  ! Do cheap calculation using linearisation about a reference temperature
  ! if the fields for this are present
  do ic = 1, n_points
    qsat(ic) = linear_qs(ic,i_qsat_liq_ref)                                    &
             + ( temperature_l(ic) - linear_qs(ic,i_ref_temp) )                &
               * linear_qs(ic,i_dqsatdT_liq)
  end do
else
  ! Otherwise calculate qsat from scratch
  call set_qsat_liq( n_points, temperature_l, pressure, qsat )
end if

! Set supersaturation = q_vap + q_cl - qsat(Tl)
do ic = 1, n_points
  qvl_supersat(ic) = q_vap_l(ic) - qsat(ic)
end do

! Calculate virtual temperature from fields with q_cl evaporated
call calc_virt_temp( n_points, n_points,                                       &
                     temperature_l, q_vap_l, q_cond_l, virt_temp_noliq )


return
end subroutine calc_qvl_supersat


end module calc_qvl_supersat_mod
