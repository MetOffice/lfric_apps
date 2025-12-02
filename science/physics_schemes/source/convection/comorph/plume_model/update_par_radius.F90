! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module update_par_radius_mod

implicit none

contains


! Subroutine to update the parcel's thermal radius following
! changes in both mass-flux and density over a level-step.
subroutine update_par_radius( n_points, det_mass_d,                            &
                              prev_temperature, next_temperature,              &
                              prev_q_vap,       next_q_vap,                    &
                              prev_pressure,    next_pressure,                 &
                              prev_massflux_d,  next_massflux_d,               &
                              prev_radius,      next_radius )

use comorph_constants_mod, only: real_cvprec, one, third, min_float,           &
                                 min_radius, par_radius_evol_method,           &
                                 par_radius_evol_const,                        &
                                 par_radius_evol_volume,                       &
                                 par_radius_evol_no_decrease,                  &
                                 par_radius_evol_no_detrain
use calc_rho_dry_mod, only: calc_rho_dry

implicit none

! Number of points
integer, intent(in) :: n_points

! Detrained mass
real(kind=real_cvprec), intent(in) :: det_mass_d(n_points)

! Parcel temperature, water-vapour mixing-ratio and pressure at prev and next
real(kind=real_cvprec), intent(in) :: prev_temperature(n_points)
real(kind=real_cvprec), intent(in) :: next_temperature(n_points)
real(kind=real_cvprec), intent(in) :: prev_q_vap(n_points)
real(kind=real_cvprec), intent(in) :: next_q_vap(n_points)
real(kind=real_cvprec), intent(in) :: prev_pressure(n_points)
real(kind=real_cvprec), intent(in) :: next_pressure(n_points)

! prev and next convective fluxes of dry-mass
real(kind=real_cvprec), intent(in) :: prev_massflux_d(n_points)
real(kind=real_cvprec), intent(in) :: next_massflux_d(n_points)

! Previous parcel radius
real(kind=real_cvprec), intent(in) :: prev_radius(n_points)
! Next parcel radius, to be calculated here
real(kind=real_cvprec), intent(out) :: next_radius(n_points)

! Parcel dry densities before and after the level-step
real(kind=real_cvprec) :: prev_rho_dry(n_points)
real(kind=real_cvprec) :: next_rho_dry(n_points)

! Loop counters
integer :: ic


if ( par_radius_evol_method==par_radius_evol_const ) then
  ! First and simplest option: just copy the previous parcel radius
  ! so that it stays constant with height

  do ic = 1, n_points
    next_radius(ic) = prev_radius(ic)
  end do

else
  ! Other options...

  ! Compute parcel dry-density at previous and next half-levels
  call calc_rho_dry( n_points, prev_temperature, prev_q_vap, prev_pressure,    &
                     prev_rho_dry )
  call calc_rho_dry( n_points, next_temperature, next_q_vap, next_pressure,    &
                     next_rho_dry )

  select case( par_radius_evol_method )
  case ( par_radius_evol_volume )
    ! Parcel radius just changes in proportion with the expected volume^(1/3)

    do ic = 1, n_points
      ! Thermal volume r^3 scales with
      !    flux of dry-mass * volume per unit dry-mass (1/rho_dry)
      next_radius(ic) = prev_radius(ic)                                        &
           * ( ( next_massflux_d(ic) / max( prev_massflux_d(ic), min_float ) ) &
             * ( prev_rho_dry(ic) / next_rho_dry(ic) ) )**third
    end do

  case ( par_radius_evol_no_decrease )
    ! As above, but not allowed to decrease with height
    ! (assumed that where this factor starts to net decrease with height,
    !  detrainment will start to selectively terminate some thermals but not
    !  others, so that the mean thermal size in the bulk plume stays the same).

    do ic = 1, n_points
      next_radius(ic) = prev_radius(ic)                                        &
        * max( ( next_massflux_d(ic) / max( prev_massflux_d(ic), min_float ) ) &
             * ( prev_rho_dry(ic) / next_rho_dry(ic) ), one )**third
    end do

  case ( par_radius_evol_no_detrain )
    ! Similar, but now we always allow mass increase by entrainment to
    ! lead to increased radius, and ignore any detrainment occuring
    ! simultaneously (assumed that detrainment acts to just selectively
    ! terminate some thermals in the ensemble, rather than homogeneously
    ! making them all smaller).

    do ic = 1, n_points
      ! Thermal volume r^3 scales with
      ! flux of dry-mass * volume per unit dry-mass (1/rho_dry)
      ! Using next mass-flux we would have had without any detrainment.
      next_radius(ic) = prev_radius(ic)                                        &
           * ( ( ( next_massflux_d(ic) + det_mass_d(ic) )                      &
                 / max( prev_massflux_d(ic), min_float ) )                     &
             * ( prev_rho_dry(ic) / next_rho_dry(ic) ) )**third
    end do

  end select  ! ( par_radius_vol_method )

end if  ! ( par_radius_evol_method==par_radius_evol_const )

! Limit the parcel radius to not fall below a hardwired
! safety-limit (some calculations that use it run into
! trouble if it gets too near zero)
do ic = 1, n_points
  next_radius(ic) = max( next_radius(ic), min_radius )
end do


return
end subroutine update_par_radius

end module update_par_radius_mod
