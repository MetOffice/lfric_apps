!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief module that contains routines taken to initialise fields based upon
!> DCMIP test 31 - non-hydrostatic gravity waves
!> @details The non-hydrostatic gravity wave test examines the response of models to short time-scale wavemotion
!> triggered by a localized perturbation. The formulation presented in this document is new,
!> but is based on previous approaches by Skamarock et al. (JAS 1994), Tomita and Satoh (FDR 2004), and
!> Jablonowski et al. (NCAR Tech Report 2008)
!>
module generate_global_gw_fields_mod

use constants_mod,                  only : r_def, pi
use initial_pressure_config_mod,    only : surface_pressure
use initial_temperature_config_mod, only : bvf_square, &
                                           pert_width_scaling, &
                                           pert_centre
use planet_config_mod,              only : gravity, &
                                           scaled_radius, scaled_omega, &
                                           Rd, Cp, p_zero, kappa, scaling_factor
use initial_wind_config_mod,        only : u0
use formulation_config_mod,         only : rotating
use log_mod,                        only : log_scratch_space, log_event, &
                                           LOG_LEVEL_ERROR

implicit none

private

public :: generate_global_gw_fields
public :: generate_global_gw_pert

contains

!> @brief Function to generate the background fields for
!> the global gravity wave test
!> @param[in] lat Latitude (radians) of the point to compute the fields at
!> @param[in] z Height (m) above the mean surface of the point to compute the fields at
!> @param[out] exner Exner pressure at the desired point
!> @param[out] u Velocity vector at the desired point
!> @param[out] theta Potential temperature at the desired point
!> @param[out] rho Density at the desired point
!>
subroutine generate_global_gw_fields (lat, z, exner, u, theta, rho)

implicit none

  real(kind=r_def), intent(in)  :: lat, z ! Latitude (radians) and Height (m)

  real(kind=r_def), intent(out) :: u(3), &               ! (Zonal,Meridional,Vertical) wind (m s^-1)
                                   theta, &              ! potential Temperature (K)
                                   exner, &              ! exner pressure
                                   rho                   ! density (kg m^-3)

  real(kind=r_def), parameter :: T_EQUATOR = 300.0_r_def,   &     ! Temperature at Equator
                                 ZTOP      = 10000.0_r_def        ! Model Top

  real(kind=r_def) :: bigG                                    ! G constant from DCMIP formulation
  real(kind=r_def) :: tsurf, psurf                            ! Surface temperature (k) and pressure (Pa)
  real(kind=r_def) :: temperature, pressure                   ! temperature(k) and pressure (Pa)
  real(kind=r_def) :: exp_fac
  real(kind=r_def) :: p_equator
  real(kind=r_def) :: u00

  p_equator = surface_pressure

! Calculate bigG
  bigG = gravity**2/(bvf_square*Cp)

! Initialise wind field
  u(1) = u0 * cos(lat)
  u(2) = 0.0_r_def
  u(3) = 0.0_r_def

!
  u00 = u0
  if ( rotating ) u00 = u00 + 2.0_r_def*scaled_omega*scaled_radius
  exp_fac = u00*(cos(2.0_r_def*lat)-1.0_r_def)

! Compute surface temperture
  tsurf = bigG + (T_EQUATOR - bigG)*exp( -(u0*bvf_square/(4.0_r_def*gravity**2))*exp_fac )


! Raise a helpful error here to prevent crashes if Tsurf is negative
  if (tsurf <= 0.0_r_def) then
    write(log_scratch_space,'(A)') 'The choice of u0 is too high for this '// &
                                   'test, it is resulting in negative ' // &
                                   'surface temperature (in Kelvins).'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
  end if

! Compute surface pressure
  psurf = p_equator*exp( (u0/(4.0_r_def*bigG*Rd))*exp_fac  ) * (tsurf/T_EQUATOR)**(Cp/Rd)

! Compute pressure and temperature
  pressure = psurf * ( (bigG/tsurf)*exp(-bvf_square*z/gravity)+1.0_r_def - (bigG/tsurf)  )**(Cp/Rd)

  temperature = bigG*(1.0_r_def - exp(bvf_square*z/gravity)) &
                + tsurf*exp(bvf_square*z/gravity)

! Compute density from equation of state
  rho = pressure/(Rd*temperature)

! Convert pressure to exner pressure and temperature to potential temperature
  exner = (pressure/p_zero)**kappa
  theta = temperature/exner

end subroutine generate_global_gw_fields

!=================================================================================

!> @brief Function to generate the potential temperature pertubation for
!> the global gravity wave test
!> @param[in] lon Longtitude (radians) of the point to compute the field at
!> @param[in] lat Latitude (radians) of the point to compute the field at
!> @param[in] z Height (m) above the mean surface of the point to compute the fields at
!> @return theta Potential temperature perturbation at the desired point
!>
pure function generate_global_gw_pert(lon, lat, z) result(theta)

implicit none

  real(kind=r_def)              :: theta
  real(kind=r_def), intent(in)  :: lon, lat, z

  real(kind=r_def) :: sin_tmp, cos_tmp, r, shape_function

  real(kind=r_def), parameter :: D       = 625000.0_r_def,             &     ! Width for Pert
                                 PHIC    = 0.0_r_def,                  &     ! Lat of Pert Center
                                 DELTA_THETA = 1.0_r_def,              &     ! Max Amplitude of Pert
                                 LZ      = 10000.0_r_def                     ! Vertical half-Wavelength of Pert

 real(kind=r_def) :: D_scaled, lambda_c

  sin_tmp = sin(lat) * sin(PHIC)
  cos_tmp = cos(lat) * cos(PHIC)

! great circle distance
  lambda_c = pert_centre/180.0_r_def * pi
  r  = scaled_radius * acos (sin_tmp + cos_tmp*cos(lon-lambda_c))
  D_scaled = pert_width_scaling*D/scaling_factor
  shape_function = (D_scaled**2)/(D_scaled**2 + r**2)

  theta = DELTA_THETA*shape_function*sin(pi*z/LZ)
end function generate_global_gw_pert

end module generate_global_gw_fields_mod

