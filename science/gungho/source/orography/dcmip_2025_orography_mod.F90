!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Calculates DCMIP 2025 orography.
!>
!> @details This module contains type definition and routines to calculate
!>          analytic orography profile of DCMIP200 case mountain function from
!>          spherical polar coordinates: longitude (lambda) and latitude (phi).
!>          Reference: Ulrich et al. (2012), Section 2.0.
!>          DCMIP200 mountain parameters in (lon,lat) coordinates are:
!>          mountain_height - Height of DCMIP200 mountain function (m),
!>          radius - Radius of DCMIP200 mountain function (radian),
!>          osc_half_width - Oscillation half-width of DCMIP200 mountain (radian),
!>          lambda_centre - Longitudinal centre of DCMIP200 mountain function (radian),
!>          phi_centre - Latitudinal centre of DCMIP200 mountain function (radian).
!-------------------------------------------------------------------------------
module dcmip_2025_orography_mod

  use constants_mod,          only : r_def, i_def
  use analytic_orography_mod, only : analytic_orography_type
  use log_mod,                only : log_event, LOG_LEVEL_ERROR

  implicit none

  private

  type, public, extends(analytic_orography_type) :: dcmip_breaking_gw_type

    private

  contains

    procedure, public, pass(self) :: analytic_orography => dcmip_breaking_gw_orography

  end type dcmip_breaking_gw_type

  type, public, extends(analytic_orography_type) :: dcmip_gap_type

    private

  contains

    procedure, public, pass(self) :: analytic_orography => dcmip_gap_orography

  end type dcmip_gap_type

  type, public, extends(analytic_orography_type) :: dcmip_vortex_type

    private

  contains

    procedure, public, pass(self) :: analytic_orography => dcmip_vortex_orography

  end type dcmip_vortex_type

  interface dcmip_breaking_gw_type
    module procedure dcmip_breaking_gw_constructor
  end interface

  interface dcmip_gap_type
    module procedure dcmip_gap_constructor
  end interface

  interface dcmip_vortex_type
    module procedure dcmip_vortex_constructor
  end interface

contains

  !=============================================================================
  type(dcmip_breaking_gw_type) function dcmip_breaking_gw_constructor( )  result(self)

    implicit none

    return
  end function dcmip_breaking_gw_constructor

  !=============================================================================
  function dcmip_breaking_gw_orography(self, chi_1, chi_2) result(chi_surf)

    use constants_mod,     only : PI

    implicit none

    ! Arguments
    class(dcmip_breaking_gw_type), intent(in) :: self
    real(kind=r_def),              intent(in) :: chi_1, chi_2

    real(kind=r_def) :: chi_surf


    call log_event('Breaking GW orography not yet implemented', LOG_LEVEL_ERROR)
    chi_surf = 0.0_r_def

    return
  end function dcmip_breaking_gw_orography


  !=============================================================================
  type(dcmip_gap_type) function dcmip_gap_constructor( )  result(self)

    implicit none

    return
  end function dcmip_gap_constructor

  !=============================================================================
  function dcmip_gap_orography(self, chi_1, chi_2) result(chi_surf)

    use constants_mod,     only : PI
    use planet_config_mod, only : scaling_factor, scaled_radius

    implicit none

    ! Arguments
    class(dcmip_gap_type), intent(in) :: self
    real(kind=r_def),      intent(in) :: chi_1, chi_2

    real(kind=r_def), parameter :: e1 = 10.0_r_def
    real(kind=r_def), parameter :: e2 = 10.0_r_def
    real(kind=r_def), parameter :: e3 = 10.0_r_def
    real(kind=r_def), parameter :: h0 = 1500.0_r_def
    real(kind=r_def), parameter :: H = 0.1_r_def*h0
    real(kind=r_def), parameter :: lon_c = PI/2.0_r_def
    real(kind=r_def), parameter :: lat_c = 0.0_r_def
    real(kind=r_def), parameter :: x_lon = 800000.0_r_def
    real(kind=r_def), parameter :: x_lat = 6000000.0_r_def
    real(kind=r_def), parameter :: x_gap = 1000000.0_r_def

    real(kind=r_def) :: d1, d2, d3
    real(kind=r_def) :: scaled_x_lon, scaled_x_lat, scaled_x_gap
    real(kind=r_def) :: chi_surf

    scaled_x_lon = x_lon / scaling_factor
    scaled_x_lat = x_lat / scaling_factor
    scaled_x_gap = x_gap / scaling_factor

    d1 = 0.5_r_def * scaled_x_lon / scaled_radius * LOG(h0/H)**(-1.0_r_def/e1)
    d2 = 0.5_r_def * scaled_x_lat / scaled_radius * LOG(h0/H)**(-1.0_r_def/e2)
    d3 = 0.5_r_def * scaled_x_gap / scaled_radius * LOG(h0/H)**(-1.0_r_def/e3)

    chi_surf = h0 * EXP(-((chi_1 - lon_c)/d1)**e1 - ((chi_2 - lat_c)/d2)**e2)  &
               * (1.0_r_def - EXP(-((chi_2 - lat_c)/d3)**e3))

    return
  end function dcmip_gap_orography

  !=============================================================================
  type(dcmip_vortex_type) function dcmip_vortex_constructor( )  result(self)

    implicit none

    return
  end function dcmip_vortex_constructor

  !=============================================================================
  function dcmip_vortex_orography(self, chi_1, chi_2) result(chi_surf)

    use constants_mod,     only : PI
    use planet_config_mod, only : scaling_factor, scaled_radius

    implicit none

    ! Arguments
    class(dcmip_vortex_type), intent(in) :: self
    real(kind=r_def),         intent(in) :: chi_1, chi_2

    real(kind=r_def), parameter :: h0 = 1500.0_r_def
    real(kind=r_def), parameter :: d = 250000.0_r_def
    real(kind=r_def), parameter :: lon_c = PI / 2.0_r_def
    real(kind=r_def), parameter :: lat_c = PI / 9.0_r_def

    real(kind=r_def) :: chi_surf, scaled_d, r

    ! Great arc distance
    r = scaled_radius * ACOS(                                                  &
        SIN(lat_c) * SIN(chi_2) + COS(lat_c) * COS(chi_2) * COS(chi_1 - lon_c) &
    )
    scaled_d = d / scaling_factor

    chi_surf = h0 * EXP(-(r/scaled_d)**2)

    return
  end function dcmip_vortex_orography

end module dcmip_2025_orography_mod

