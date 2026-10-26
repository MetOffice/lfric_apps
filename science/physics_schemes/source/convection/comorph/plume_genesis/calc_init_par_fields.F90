! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module calc_init_par_fields_mod

implicit none

contains

! Subroutine to set the parcel initial perturbations for temperature
! and water-vapour, and apply safety-checks to the perturbations.
subroutine calc_init_par_fields( n_points, nc2, index_ic2, l_down, i_region,   &
                                 pressure_k, cf_liq,                           &
                                 turb_tl, turb_qt,                             &
                                 delta_temp_neut, delta_qvap_neut,             &
                                 rhpert_t, fields_par, pert_tl, pert_qt )

use comorph_constants_mod, only: real_cvprec, zero, one, sqrt_min_delta,       &
                                 l_turb_par_gen, par_gen_core_fac,             &
                                 par_gen_rhpert
use subregion_mod, only: i_liq, i_mph
use fields_type_mod, only: n_fields, i_temperature, i_q_vap
use set_qsat_mod, only: set_qsat_liq

implicit none

! Total number of points being processed in init_mass_moist_frac
integer, intent(in) :: n_points

! Number of points where initiating mass is occuring in the
! current sub-grid region
integer, intent(in) :: nc2

! Indices of the above points, for compression / decompression
integer, intent(in) :: index_ic2(nc2)

! Flag for downdrafts vs updrafts
logical, intent(in) :: l_down

! Indicator for which sub-grid partition we're in
integer, intent(in) :: i_region

! Level k pressure compressed onto initiating points
real(kind=real_cvprec), intent(in) :: pressure_k(nc2)

! Total liquid-cloud fraction
real(kind=real_cvprec), intent(in) :: cf_liq(n_points)

! Turbulence-based perturbations to Tl and qt
real(kind=real_cvprec), intent(in) :: turb_tl(n_points)
real(kind=real_cvprec), intent(in) :: turb_qt(n_points)

! Neutrally buoyant perturbations to T,q so-as to yield
! a specified RH perturbation
real(kind=real_cvprec), intent(in) :: delta_temp_neut(n_points)
real(kind=real_cvprec), intent(in) :: delta_qvap_neut(n_points)

! RH perturbation for current convection type
real(kind=real_cvprec), intent(in) :: rhpert_t(n_points)

! Unperturbed initiating parcel properties
real(kind=real_cvprec), intent(in) :: fields_par                               &
                            ( n_points, i_temperature:n_fields )
! Perturbations applied to initiating parcel Tl and qt
real(kind=real_cvprec), intent(out) :: pert_tl(nc2)
real(kind=real_cvprec), intent(out) :: pert_qt(nc2)

! Relative Humidity of the parcel
real(kind=real_cvprec) :: rh_par(nc2)
! Relative Humidity and temperature of the neutrally-buoyant perturbed parcel
real(kind=real_cvprec) :: rh_par_pert(nc2)
real(kind=real_cvprec) :: temperature(nc2)
real(kind=real_cvprec) :: q_vap(nc2)
! Factor for scaling the neutrally-buoyant moisture perturbations
real(kind=real_cvprec) :: factor(nc2)
! Factor for the neutrally buoyant core moisture perturbations
real(kind=real_cvprec) :: fac_core

! Loop counters
integer :: ic, ic2


if ( l_turb_par_gen ) then
  ! Using turbulence-based perturbations...

  ! Set Tl and qt perturbations equal to turbulence-based ones.
  ! (sign of perturbations is reversed for downdrafts)
  if ( l_down ) then
    do ic2 = 1, nc2
      ic = index_ic2(ic2)
      pert_tl(ic2) = -turb_tl(ic)
      pert_qt(ic2) = -turb_qt(ic)
    end do
  else
    do ic2 = 1, nc2
      ic = index_ic2(ic2)
      pert_tl(ic2) = turb_tl(ic)
      pert_qt(ic2) = turb_qt(ic)
    end do
  end if

  ! If triggering from liquid-cloud with low grid-mean RH / small
  ! cloud-fraction, the moisture excess of the liquid-cloud relative
  ! to the grid-mean will likely already reflect the turbulent moisture
  ! fluctuations.  So adding the turbulent q perturbation onto the
  ! in-cloud q is double-counting.  In the opposite limit
  ! (grid-mean saturation with liquid cloud fraction of 1), the in-cloud
  ! q doesn't include any turbulent fluctuation, so we aren't
  ! double-counting.  Therefore, blend the turbulent perturbation as a
  ! function of the liquid cloud fraction.
  if ( i_region==i_liq .or. i_region==i_mph ) then
    do ic2 = 1, nc2
      ic = index_ic2(ic2)
      pert_tl(ic2) = pert_tl(ic2) * cf_liq(ic)*cf_liq(ic)
      pert_qt(ic2) = pert_qt(ic2) * cf_liq(ic)*cf_liq(ic)
    end do
  else
    do ic2 = 1, nc2
      ic = index_ic2(ic2)
      pert_tl(ic2) = pert_tl(ic2) * (one-cf_liq(ic))*(one-cf_liq(ic))
      pert_qt(ic2) = pert_qt(ic2) * (one-cf_liq(ic))*(one-cf_liq(ic))
    end do
  end if

else  ! ( l_turb_par_gen )
  ! Not using turbulence-based perturbations...

  ! Set Tl and qt perturbations equal to zero
  do ic2 = 1, nc2
    pert_tl(ic2) = zero
    pert_qt(ic2) = zero
  end do

end if  ! ( l_turb_par_gen )


if ( par_gen_rhpert > zero ) then
  ! Add on neutrally-buoyant T,q perturbations consistent
  ! with a specified Relative Humidity perturbation.
  ! T,q perturbations which increase RH by a fixed amount whilst
  ! preserving Tv were calculated in calc_env_regions
  ! and passed in here as delta_temp_neut, delta_qvap_neut...

  ! Set factor to scale them by to get the input variable RH perturbation
  do ic2 = 1, nc2
    ic = index_ic2(ic2)
    factor(ic2) = rhpert_t(ic) / par_gen_rhpert
  end do

  ! Calculate relative humidity of the parcel
  call set_qsat_liq( nc2, fields_par(:,i_temperature), pressure_k, rh_par )
  do ic2 = 1, nc2
    rh_par(ic2) = fields_par(ic2,i_q_vap) / rh_par(ic2)
  end do

  ! Calculate relative humidity of the parcel with the T,q perturbations
  ! that will be applied to the parcel core (scaled up by core_fac)
  do ic2 = 1, nc2
    ic = index_ic2(ic2)
    fac_core = factor(ic2) * par_gen_core_fac
    ! TEMPORARY CODE TO FORCE BIT-REPRODUCIBILITY WITH PREVIOUS VERSION
    if ( abs(fac_core-one) < sqrt_min_delta )  fac_core = one
    temperature(ic2) = fields_par(ic2,i_temperature) + delta_temp_neut(ic)     &
                                                       * fac_core
    q_vap(ic2)       = fields_par(ic2,i_q_vap)       + delta_qvap_neut(ic)     &
                                                       * fac_core
  end do
  call set_qsat_liq( nc2, temperature, pressure_k, rh_par_pert )
  do ic2 = 1, nc2
    rh_par_pert(ic2) = q_vap(ic2) / rh_par_pert(ic2)
  end do

  do ic2 = 1, nc2
    ! TEMPORARY CODE TO FORCE BIT-REPRODUCIBILITY WITH PREVIOUS VERSION
    factor(ic2) = factor(ic2) * par_gen_core_fac
    if ( abs(factor(ic2)-one) < sqrt_min_delta )  factor(ic2) = one
    ! Apply only a fraction of the neutrally-buoyant relative humidity
    ! perturbation if it would make the parcel core supersaturated
    if ( rh_par(ic2) >= one ) then
      ! If parcel is already saturated, set RH perturbation to zero
      factor(ic2) = zero
    else if ( rh_par_pert(ic2) > one ) then
      ! If parcel was subsaturated, but perturbation makes it supersaturated
      ! calculate fraction of perturbation needed to reach saturation
      factor(ic2) = factor(ic2) * ( one - rh_par(ic2) )                        &
                                / ( rh_par_pert(ic2) - rh_par(ic2) )
      ! Note: due to loss of precision, it is possible to end up with
      ! rh_par_pert == rh_par.
      ! The above logic ensures we cannot enter this branch and get a
      ! div-by-zero when this happens.
    end if
    ! TEMPORARY CODE TO FORCE BIT-REPRODUCIBILITY WITH PREVIOUS VERSION
    factor(ic2) = factor(ic2) / par_gen_core_fac
  end do

  ! Add on the scaled neutrally-buoyant RH perturbations
  do ic2 = 1, nc2
    ic = index_ic2(ic2)
    pert_tl(ic2) = pert_tl(ic2) + delta_temp_neut(ic) * factor(ic2)
    pert_qt(ic2) = pert_qt(ic2) + delta_qvap_neut(ic) * factor(ic2)
  end do

end if  ! ( par_gen_rhpert > zero )

! Limit the q perturbation to avoid creating negative or
! excessively large moisture contents.
! We want:
! par_gen_core_fac ABS(pert_qt) <= par_q (1 - sqrt(epsilon))
do ic2 = 1, nc2
  factor(ic2) = fields_par(ic2,i_q_vap) * (one-sqrt_min_delta)                 &
                                        / par_gen_core_fac
end do
do ic2 = 1, nc2
  if ( abs(pert_qt(ic2)) > factor(ic2) ) then
    ! Fortran SIGN intrinsic returns factor but with same sign as pert_qt.
    pert_qt(ic2) = sign( factor(ic2), pert_qt(ic2) )
  end if
end do

return
end subroutine calc_init_par_fields


end module calc_init_par_fields_mod
