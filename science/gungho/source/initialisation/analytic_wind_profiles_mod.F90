!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

!> @brief Functions to compute a variety of analytic profiles
!> @details Collection of functions to return the value of a field at a given
!!          point based upon a specified analytic formula
module analytic_wind_profiles_mod

use constants_mod,            only: r_def, pi, EPS, i_def
use initial_wind_config_mod,  only: profile_none,                         &
                                    profile_solid_body_rotation,          &
                                    profile_solid_body_rotation_alt,      &
                                    profile_constant_uv,                  &
                                    profile_constant_shear_uv,            &
                                    profile_dcmip301,                     &
                                    profile_sin_uv,                       &
                                    profile_deep_baroclinic_steady,       &
                                    profile_deep_baroclinic_perturbed,    &
                                    profile_vortex,                       &
                                    profile_xy_NL_case_1,                 &
                                    profile_yz_NL_case_1,                 &
                                    profile_NL_case_1,                    &
                                    profile_NL_case_2,                    &
                                    profile_NL_case_3,                    &
                                    profile_NL_case_4,                    &
                                    profile_hadley_like_dcmip,            &
                                    profile_curl_free_reversible,         &
                                    profile_sbr_with_vertical,            &
                                    profile_dcmip_101,                    &
                                    profile_vertical_deformation,         &
                                    profile_four_part_sbr,                &
                                    wind_time_period
use planet_config_mod,        only: scaled_radius
use log_mod,                  only: log_event,                         &
                                    log_scratch_space,                 &
                                    LOG_LEVEL_ERROR
use deep_baroclinic_wave_mod, only: deep_baroclinic_wave
use formulation_config_mod,   only: shallow

implicit none

private :: vortex_wind
private :: NL_wind_case_1
private :: NL_wind_case_2
private :: NL_wind_case_3
private :: NL_wind_case_4
private :: xy_NL_wind_case_1
private :: yz_NL_wind_case_1
private :: hadley_like_dcmip
private :: xy2longlat
private :: yz2longlat
private :: curl_free_reversible

public :: analytic_wind

contains

!> @brief Compute wind field for Hadley-like DCMIP test
!> @param[in] lat Latitudinal position
!> @param[in] height Radial distance
!> @param[in] time Time
!> @result u Result wind field vector (u,v,w)
function hadley_like_dcmip(lat,height,time) result(u)
  implicit none
  real(kind=r_def), intent(in)    :: lat
  real(kind=r_def), intent(in)    :: height
  real(kind=r_def), intent(in)    :: time
  real(kind=r_def), dimension(3)  :: u

  real(kind=r_def)             :: u0, w0, k, top_of_atmosphere, l, z
  real(kind=r_def)             :: tau, time_period, rho, scale_height

  ! Equations below have been taken from Allen and Zerroukat, "A deep
  !> non-hydrostatic compressible atmospheric model on a Yin-Yang grid", JCP, 2016,
  !> or equivalently Kent, Ullrich, Jablonowski, "Dynamical core intercomparison project:
  !> Tracer transport test cases", QJRMS 2014.
  !> Note that there is a missing negative sign in Allen and Zerroukat,
  !> equation (5.2).

  u0 = 40.0_r_def
  w0 = 0.03_r_def
  k = 5.0_r_def
  top_of_atmosphere = 12000.0_r_def
  l = pi/top_of_atmosphere
  time_period = 24.0_r_def*60.0_r_def*60.0_r_def
  tau = pi/time_period
  z = height-scaled_radius
  scale_height = 300.0_r_def*287.0_r_def/9.80616_r_def
  rho = exp(-z/scale_height)

  u(1) = u0*cos(lat)
  u(2) = (1.0_r_def/rho)*(-scaled_radius*w0*l)*cos(lat)*sin(k*lat)*cos(l*z)*cos(tau*time)
  u(3) = (1.0_r_def/rho)*w0*(-2.0_r_def*sin(k*lat)*sin(lat)+k*cos(k*lat)*cos(lat))*sin(l*z)*cos(tau*time)

end function hadley_like_dcmip

!> @brief Compute a vortex wind field
!> @param[in] lat Latitudinal position
!> @param[in] long Longitudinal position
!> @param[in] radius Radial distance
!> @result u Result wind field vector (u,v,w)
function vortex_wind(lat,long,radius) result(u)
  implicit none
  real(kind=r_def), intent(in)    :: lat
  real(kind=r_def), intent(in)    :: long
  real(kind=r_def), intent(in)    :: radius
  real(kind=r_def), dimension(3)  :: u

  real(kind=r_def) :: lat_pole, lon_pole, lat_dash
  real(kind=r_def) :: r0, v0, V
  real(kind=r_def) :: radial_distance, omega

  ! Equations below have been taken from Nair and Jablonowski, "Moving vortices
  ! on the sphere: A test case for horizontal advection problems", AMS 2008,
  ! equations (20) and (21)
  ! lat_pole and lon_pole denote the position of one of the vortices, the other
  ! vortex lies on the opposite side of the sphere
  ! Parameter values have been taken from the paper and are currently
  ! hard-wired

  lat_pole = 0.0_r_def
  lon_pole = 0.0_r_def

  lat_dash = asin(sin(lat)*sin(lat_pole) + cos(lat)*cos(lat_pole)*cos(long-lon_pole))

  r0 = 3.0_r_def
  v0 = 38.61073731_r_def
  radial_distance = r0*cos(lat_dash)
  V = v0*3.0_r_def*(sqrt(3.0_r_def)/2.0_r_def)*(sinh(radial_distance)/cosh(radial_distance)**3)

  if (abs(radial_distance)<1E-10_r_def) then
    omega = 0.0_r_def
  else
    omega = V/(radius*radial_distance)
  end if

  u(1) = radius*omega*(sin(lat_pole)*cos(lat) - cos(lat_pole)*cos(long-lon_pole)*sin(lat))
  u(2) = radius*omega*cos(lat_pole)*sin(long-lon_pole)
  u(3) = 0.0_r_def
end function vortex_wind

!> @brief Compute Case 1 of Nair and Lauritzen JCP 229 (2010)
!> @param[in] long Longitudinal position in spherical coordinates
!> @param[in] lat  Latitudinal position in spherical coordinates
!> @param[in] time Current time (timestep multiplied by dt)
!> @result    u    Wind field vector (u,v,w) in spherical coordinates
function NL_wind_case_1(long,lat,time) result(u)
  use initial_wind_config_mod, only : wind_time_period, nl_constant

  implicit none
  real(kind=r_def), intent(in)    :: long
  real(kind=r_def), intent(in)    :: lat
  real(kind=r_def), intent(in)    :: time
  real(kind=r_def), dimension(3)  :: u

  real(kind=r_def) :: u0

  ! Equations below have been taken from Case 1 of Nair and Lauritzen, 2010
  ! "A class of deformational flow test cases for linear transport problems on the sphere"

  ! Coefficents for any sized planet
  u0 = nl_constant*scaled_radius/wind_time_period

  u(1) = u0*sin(long/2.0_r_def)**2*sin(2.0_r_def*lat)*cos(pi*time/wind_time_period)
  u(2) = (u0/2.0_r_def)*sin(long)*cos(lat)*cos(pi*time/wind_time_period)
  u(3) = 0.0_r_def

end function NL_wind_case_1

!> @brief Compute Case 2 of Nair and Lauritzen JCP 229 (2010)
!> @param[in] long Longitudinal position in spherical coordinates
!> @param[in] lat  Latitudinal position in spherical coordinates
!> @param[in] time Current time (timestep multiplied by dt)
!> @result    u    Wind field vector (u,v,w) in spherical coordinates
function NL_wind_case_2(long,lat,time) result(u)
  use initial_wind_config_mod, only : wind_time_period, nl_constant

  implicit none
  real(kind=r_def), intent(in)    :: long
  real(kind=r_def), intent(in)    :: lat
  real(kind=r_def), intent(in)    :: time
  real(kind=r_def), dimension(3)  :: u

  real(kind=r_def) :: u0

  ! Equations below have been taken from Case 2 of Nair and Lauritzen, 2010
  ! "A class of deformational flow test cases for linear transport problems on the sphere"

  ! Coefficents for any sized planet
  u0 = nl_constant*scaled_radius/wind_time_period

  u(1) = u0*sin(long)**2*sin(2.0_r_def*lat)*cos(pi*time/wind_time_period)
  u(2) = u0*sin(2.0_r_def*long)*cos(lat)*cos(pi*time/wind_time_period)
  u(3) = 0.0_r_def

end function NL_wind_case_2

!> @brief Compute Case 3 of Nair and Lauritzen JCP 229 (2010)
!> @param[in] long Longitudinal position in spherical coordinates
!> @param[in] lat  Latitudinal position in spherical coordinates
!> @param[in] time Current time (timestep multiplied by dt)
!> @result    u    Wind field vector (u,v,w) in spherical coordinates
function NL_wind_case_3(long,lat,time) result(u)
  use initial_wind_config_mod, only : wind_time_period, nl_constant

  implicit none
  real(kind=r_def), intent(in)    :: long
  real(kind=r_def), intent(in)    :: lat
  real(kind=r_def), intent(in)    :: time
  real(kind=r_def), dimension(3)  :: u

  real(kind=r_def) :: long_dash
  real(kind=r_def) :: u0, v0

  ! Equations below have been taken from Case 3 of Nair and Lauritzen, 2010
  ! "A class of deformational flow test cases for linear transport problems on the sphere"
  ! with an additional background flow

  long_dash = long-2.0_r_def*pi*time/wind_time_period

  ! Coefficents for any sized planet
  u0 = nl_constant*scaled_radius/wind_time_period
  v0 = 2.0_r_def*pi*scaled_radius/wind_time_period

  u(1) = -u0*sin(long_dash/2.0_r_def)**2*sin(2.0_r_def*lat)*cos(lat)**2 &
            *cos(pi*time/wind_time_period) + v0*cos(lat)
  u(2) = (u0/2.0_r_def)*sin(long_dash)*cos(lat)**3*cos(pi*time/wind_time_period)
  u(3) = 0.0_r_def

end function NL_wind_case_3

!> @brief Compute Case 4 of Nair and Lauritzen JCP 229 (2010)
!> @param[in] long Longitudinal position in spherical coordinates
!> @param[in] lat  Latitudinal position in spherical coordinates
!> @param[in] time Current time (timestep multiplied by dt)
!> @result    u    Wind field vector (u,v,w) in spherical coordinates
function NL_wind_case_4(long,lat,time) result(u)
  use initial_wind_config_mod, only : wind_time_period, nl_constant

  implicit none
  real(kind=r_def), intent(in)    :: long
  real(kind=r_def), intent(in)    :: lat
  real(kind=r_def), intent(in)    :: time
  real(kind=r_def), dimension(3)  :: u

  real(kind=r_def) :: long_dash
  real(kind=r_def) :: u0, v0

  ! Equations below have been taken from Case 4 of Nair and Lauritzen, 2010
  ! "A class of deformational flow test cases for linear transport problems on the sphere"

  long_dash = long-2.0_r_def*pi*time/wind_time_period

  ! Coefficents for any sized planet
  u0 = nl_constant*scaled_radius/wind_time_period
  v0 = 2.0_r_def*pi*scaled_radius/wind_time_period

  u(1) = u0*sin(long_dash)**2*sin(2.0_r_def*lat)*cos(pi*time/wind_time_period) + v0*cos(lat)
  u(2) = u0*sin(2.0_r_def*long_dash)*cos(lat)*cos(pi*time/wind_time_period)
  u(3) = 0.0_r_def

end function NL_wind_case_4

!> @brief Remap biperiodic (x,y) values to longitude/latitude coordinates
!> @param[in] x x position in biperiodic mesh
!> @param[in] y y position in biperiodic mesh
!> @param[in] domain_max_x Domain maximum x-coordinate.
!> @param[in] domain_max_y Domain maximum y-coordinate.
!> @result xy2longlat Vector of (long,lat) values which have been remapped
function xy2longlat(x, y, domain_max_x, domain_max_y)
  ! Rescales biperiodic x,y coordinates into lat,lon coordinates to enable
  ! reuse of functions for initialising winds defined in lat, lon coordinates.
  ! The x coordinate is modified such that it is in the interval [0,2*pi]
  ! The y coordinate is modified such that it is in the interval [-pi/2,pi/2]

  implicit none
  real(kind=r_def), intent(in) :: x
  real(kind=r_def), intent(in) :: y
  real(kind=r_def), intent(in) :: domain_max_x
  real(kind=r_def), intent(in) :: domain_max_y

  real(kind=r_def) :: xy2longlat(2)

  xy2longlat(1) = pi*(x/domain_max_x+1.0_r_def)     ! x --> longitude
  xy2longlat(2) = (pi/2.0_r_def)*y/domain_max_y     ! y --> latitude

end function xy2longlat

!> @brief Case 1 of Nair and Lauritzen test cases but on a x-y slice of the
!>        biperiodic domain
!> @param[in] x x position in biperiodic mesh
!> @param[in] y y position in biperiodic mesh
!> @param[in] time Time (timestep multiplied by dt)
!> @param[in] domain_max_x Domain maximum x-coordinate.
!> @param[in] domain_max_y Domain maximum y-coordinate.
!> @result u Wind field vector (u,v,w)
function xy_NL_wind_case_1(x, y, time, domain_max_x, domain_max_y) result(u)
  ! This function is designed to be used on the biperiodic mesh and defines a
  ! wind in the x-y plane. Inputs are (x,y) which are the coordinates on the
  ! biperiodic mesh.

  implicit none
  real(kind=r_def), intent(in) :: x
  real(kind=r_def), intent(in) :: y
  real(kind=r_def), intent(in) :: time
  real(kind=r_def), intent(in) :: domain_max_x
  real(kind=r_def), intent(in) :: domain_max_y

  real(kind=r_def), dimension(3)  :: u

  real(kind=r_def)                :: longlat(2)

  longlat = xy2longlat(x, y, domain_max_x, domain_max_y)
  u = NL_wind_case_1(longlat(1),longlat(2),time)

end function xy_NL_wind_case_1

!> @brief Remap biperiodic (y,z) values to longitude/latitude coordinates
!> @param[in] y y position in biperiodic mesh
!> @param[in] z z position in biperiodic mesh
!> @param[in] domain_max_y Domain maximum y-coordinate.
!> @result yz2longlat Vector of (long,lat) values which have been remapped
function yz2longlat(y, z, domain_max_y)
  ! Rescales biperiodic y,z coordinates into lat,lon coordinates to enable
  ! reuse of functions for initialising winds defined in lat, lon coordinates.
  ! The z coordinate is modified such that it is in the interval [-pi/2,pi/2]
  ! The y coordinate is modified such that it is in the interval [0,2*pi]

  use extrusion_config_mod, only : domain_height

  implicit none
  real(kind=r_def), intent(in) :: y
  real(kind=r_def), intent(in) :: z
  real(kind=r_def), intent(in) :: domain_max_y

  real(kind=r_def) :: yz2longlat(2)

  yz2longlat(2) = pi*z/domain_height - pi/2.0_r_def ! z --> latitude
  yz2longlat(1) = pi*(y/domain_max_y+1.0_r_def)  ! y --> longitude

end function yz2longlat

!> @brief Case 1 of Nair and Lauritzen test cases but on a y-z slice of the
!>        biperiodic domain
!> @param[in] x x position in biperiodic mesh
!> @param[in] y y position in biperiodic mesh
!> @param[in] time Time (timestep multiplied by dt)
!> @param[in] domain_max_y Domain maximum y-coordinate.
!> @result u Wind field vector (u,v,w)
function yz_NL_wind_case_1(y, z, time, domain_max_y) result(u)
  ! This function is designed to be used on the biperiodic mesh and defines a
  ! wind in the y-z plane. Inputs are (y,z) which are the coordinates on the
  ! biperiodic mesh.

  implicit none
  real(kind=r_def), intent(in)    :: y
  real(kind=r_def), intent(in)    :: z
  real(kind=r_def), intent(in)    :: time
  real(kind=r_def), intent(in)    :: domain_max_y

  real(kind=r_def), dimension(3)  :: u
  real(kind=r_def), dimension(3)  :: u_temp

  real(kind=r_def)    :: longlat(2)

  longlat = yz2longlat(y, z, domain_max_y)

  u_temp = NL_wind_case_1(longlat(1),longlat(2),time)
  u(1) = 0.0_r_def
  u(2) = u_temp(1)
  u(3) = u_temp(2)

end function yz_NL_wind_case_1

!> @brief Reversible time-varying flow in a vertical slice,
!>        whose winds are the gradient of a scalar potential
!> @param[in] domain_max_x Domain maximum x-coordinate.
!> @result u Wind field vector (u,v,w)
function curl_free_reversible(x, z, t, domain_max_x) result(u)

  use extrusion_config_mod,   only : domain_height

  implicit none
  real(kind=r_def), intent(in)    :: x
  real(kind=r_def), intent(in)    :: z
  real(kind=r_def), intent(in)    :: t
  real(kind=r_def), intent(in)    :: domain_max_x
  real(kind=r_def), dimension(3)  :: u
  real(kind=r_def)                :: domain_length, time_period

  domain_length = 2.0_r_def * domain_max_x
  time_period = domain_length  ! Profile advected round once

  ! This profile will only work if the domain height is big enough
  if (domain_length > sqrt(2.0_r_def) * domain_height + EPS) then
    call log_event('Curl-free reversible wind profile only ' // &
                   'works when height >= length', LOG_LEVEL_ERROR)
  end if

  u(1) = (domain_length / time_period) * (1.0_r_def - 0.5_r_def               &
          * (t / time_period - 0.5_r_def)                                     &
          * sin(4.0_r_def * pi * (x / domain_length - t / time_period))       &
          * cos(2.0_r_def * pi * z / domain_height))
  u(2) = 0.0_r_def
  u(3) = - 0.25_r_def * domain_length**2.0_r_def / (domain_height * time_period) &
          * (t / time_period - 0.5_r_def)                                     &
          * cos(4.0_r_def * pi * (x / domain_length - t / time_period))       &
          * sin(2.0_r_def * pi * z / domain_height)

end function curl_free_reversible


!> @brief Compute an analytic wind field
!> @param[in] chi Position in physical coordinates
!> @param[in] time Time (timestep multiplied by dt)
!> @param[in] choice Integer defining which specified formula to use
!> @param[in] num_options Number of scalar options to supply
!> @param[in] domain_max_x Domain maximum x-coordinate.
!> @param[in] domain_max_y Domain maximum y-coordinate.
!> @param[in] option Array of real values used to generate the initial profile
!> @result u Result wind field vector (u,v,w)
function analytic_wind( chi, time, choice, num_options, &
                        domain_max_x, domain_max_y,     &
                        option_arg ) result(u)

  use extrusion_config_mod,   only : domain_height

  implicit none

  real(kind=r_def),           intent(in) :: chi(3)
  real(kind=r_def),           intent(in) :: time
  integer,                    intent(in) :: choice
  integer,                    intent(in) :: num_options
  real(kind=r_def),           intent(in) :: domain_max_x
  real(kind=r_def),           intent(in) :: domain_max_y
  real(kind=r_def), optional, intent(in) :: option_arg(num_options)

  ! Local variables
  real(kind=r_def)             :: option(num_options)
  real(kind=r_def)             :: u(3)
  real(kind=r_def)             :: s, r_on_a
  real(kind=r_def)             :: pressure, temperature, density, mr_v
  real(kind=r_def)             :: lat_pole, lon_pole
  real(kind=r_def)             :: k_z                 ! vertical wavenumber
  real(kind=r_def)             :: z                   ! height in spherical coordinates
  real(kind=r_def)             :: long_dash           ! translated longitude
  real(kind=r_def)             :: end_time            ! length of simulation
  real(kind=r_def)             :: u0, w0              ! hardwired velocity values
  real(kind=r_def)             :: u_d                 ! divergence correcting velocity
  real(kind=r_def)             :: p, ptop, b, rhoz    ! pressure variables for DCMIP 1.1
  real(kind=r_def), parameter  :: p0 = 100000.0_r_def ! surface pressure for DCMIP 1.1
  real(kind=r_def), parameter  :: Rd = 287.0_r_def    ! gas constant for DCMIP 1.1
  real(kind=r_def), parameter  :: T0 = 300.0_r_def    ! ref temperature for DCMIP 1.1
  real(kind=r_def), parameter  :: g  = 9.80616_r_def  ! gravity for DCMIP 1.1
  integer(kind=i_def)          :: num_periods

  if ( .not. present(option_arg) ) then
    option(:) = 0.0_r_def
  else
    option(:) = option_arg(:)
  end if

  select case ( choice )

    case ( profile_none )
      u(:) = 0.0_r_def
    case ( profile_solid_body_rotation,                           &
           profile_dcmip301)
      s = 0.5_r_def*(chi(3)/scaled_radius + 1.0_r_def)
      ! No height variation for the dcmip test
      if ( choice == profile_dcmip301) then
        s = 1.0_r_def
      else
        s = 0.5_r_def * (chi(3) / scaled_radius + 1.0_r_def)
      endif
      lat_pole = pi / 2.0_r_def - option(2) * pi
      lon_pole = pi / 2.0_r_def + option(3) * pi

      u(1) = s * option(1) *                                                   &
             (  sin(lat_pole) * cos(chi(2))                                    &
              - cos(lat_pole) * cos(chi(1)-lon_pole) * sin(chi(2)) )
      u(2) = s * option(1) * cos(lat_pole) * sin(chi(1)-lon_pole)
      u(3) = 0.0_r_def
    case ( profile_solid_body_rotation_alt)
      ! Rotated pole version of equation (74) of Staniforth & White (2007)
      ! with m=1, A=0, n therefore arbitrary.
      ! In shallow geometry the r/a factor is replaced by 1 (see section 6
      ! of the reference)
      if (shallow) then
        r_on_a = 1.0_r_def
      else
        r_on_a = chi(3) / scaled_radius
      endif
      lat_pole = pi / 2.0_r_def - option(2) * pi
      lon_pole = pi / 2.0_r_def + option(3) * pi
      u(1) = r_on_a * option(1) *                                              &
             (  sin(lat_pole) * cos(chi(2))                                    &
              - cos(lat_pole) * cos(chi(1)-lon_pole) * sin(chi(2)) )
      u(2) = r_on_a * option(1) * cos(lat_pole) * sin(chi(1)-lon_pole)
      u(3) = 0.0_r_def
    case ( profile_four_part_sbr )
      ! Test from Bendall & Wimmer (2022) involving rotation of a tracer around
      ! a sphere by composing four half-complete solid body rotations. When
      ! transporting a vector-valued field, the final state should equal the
      ! initial state.
      num_periods = FLOOR(time / wind_time_period)
      if ( MOD(num_periods, 2) == 0 ) then
        lat_pole = pi / 2.0_r_def
        lon_pole = 0.0_r_def
      else
        lat_pole = 0.0_r_def
        lon_pole = -pi / 2.0_r_def
      end if

      u(1) = option(1) *                                                       &
             (  sin(lat_pole) * cos(chi(2))                                    &
              - cos(lat_pole) * cos(chi(1)-lon_pole) * sin(chi(2)) )
      u(2) = option(1) * cos(lat_pole) * sin(chi(1)-lon_pole)
      u(3) = 0.0_r_def
    case ( profile_constant_uv )
      u(1) = option(1)
      u(2) = option(2)
      u(3) = 0.0_r_def
    case ( profile_constant_shear_uv )
      u(1) = option(1)*chi(3)/option(3)
      u(2) = option(2)*chi(3)/option(3)
      u(3) = 0.0_r_def
    case ( profile_sin_uv )
      k_z = 2.0_r_def*pi/option(3)
      u(1) = option(1)*sin(chi(3)*k_z)
      u(2) = option(2)*sin(chi(3)*k_z)
      u(3) = 0.0_r_def
    case ( profile_deep_baroclinic_steady, &
           profile_deep_baroclinic_perturbed)
      call deep_baroclinic_wave(chi(1), chi(2), chi(3)-scaled_radius, &
                                pressure, temperature, density, &
                                u(1), u(2), u(3), mr_v)
    case ( profile_vortex )
      u = vortex_wind(chi(2),chi(1),chi(3))
    case ( profile_NL_case_1 )
      u = NL_wind_case_1(chi(1),chi(2),time)
    case ( profile_NL_case_2 )
      u = NL_wind_case_2(chi(1),chi(2),time)
    case ( profile_NL_case_3 )
      u = NL_wind_case_3(chi(1),chi(2),time)
    case ( profile_NL_case_4 )
      u = NL_wind_case_4(chi(1),chi(2),time)
    case ( profile_xy_NL_case_1 )
      u = xy_NL_wind_case_1( chi(1), chi(2), time, &
                             domain_max_x, domain_max_y )
    case ( profile_yz_NL_case_1 )
      u = yz_NL_wind_case_1( chi(2), chi(3), time, &
                             domain_max_y )
    case ( profile_hadley_like_dcmip )
      u = hadley_like_dcmip(chi(2),chi(3),time)
    case ( profile_curl_free_reversible )
      u = curl_free_reversible( chi(1), chi(3), time, &
                                domain_max_x )
    case (profile_sbr_with_vertical)
      ! Solid body rotation around equator with vertical motion transport test.
      ! Set test case constants
      end_time  = 1036800.0_r_def
      u0        = 2.0_r_def * pi * scaled_radius / end_time
      w0        = 0.025_r_def
      ! Get height
      z         = chi(3) - scaled_radius
      ! Translated longitude
      long_dash = chi(1) - 2.0_r_def * pi * time / end_time
      ! Set divergent correction velocity
      u_d = - scaled_radius * cos(chi(2)) * w0 * pi / domain_height &
            * cos(pi * z / domain_height) * sin(long_dash)          &
            * cos(2.0_r_def * pi * time / end_time)
      ! Set velocity
      u(1) = u0 *  cos(chi(2)) + u_d
      u(2) = 0.0_r_def
      u(3) = w0 * sin(pi * z / domain_height) * cos(long_dash) &
             * cos(2.0_r_def * pi * time / end_time)
    case (profile_dcmip_101)
      ! DCMIP 1.1 transport test.
      ! Set test case constants
      end_time  = 1036800.0_r_def
      u0        = 10.0_r_def * scaled_radius / end_time
      w0        = 23000.0_r_def * pi / end_time
      b         = 0.2_r_def
      ! Get height
      z         = chi(3) - scaled_radius
      ! Translated longitude
      long_dash = chi(1) - 2.0_r_def * pi * time / end_time
      ! Pressure related variables
      p         = p0*exp(-g * z / ( Rd * T0 ))
      ptop      = p0*exp(-g * domain_height / ( Rd * T0 ))
      rhoz      = p / ( Rd * T0 )
      ! Set divergent correction velocity
      u_d = scaled_radius * w0 / (b * ptop) * cos(chi(2)) * cos(chi(2))      &
            * cos(long_dash) * cos(2.0_r_def * pi * time / end_time)         &
            * ( -exp( (p - p0)/(b * ptop) ) + exp( (ptop - p)/(b * ptop) ) )
      ! Set velocity
      u(1) = u0 * sin(long_dash) * sin(long_dash) * sin(2.0_r_def * chi(2))  &
                * cos(pi * time / end_time) + u_d                            &
                + 2.0_r_def * pi * scaled_radius / end_time * cos(chi(2))
      u(2) = u0 * sin(2.0_r_def * long_dash) * cos(chi(2)) * cos(pi * time / end_time)
      u(3) = - w0 / ( g*rhoz) * cos(chi(2)) * sin(long_dash)                           &
                * cos(2.0_r_def * pi * time / end_time)                                &
                * ( 1.0_r_def + exp( (ptop-p0) / (b*ptop) ) - exp( (p-p0) / (b*ptop) ) &
                - exp( (ptop-p) / (b*ptop) ) )
    case (profile_vertical_deformation)
      ! Vertical deformation transport test.
      ! Set test case constants
      w0        = 0.02_r_def
      end_time  = 1036800.0_r_def
      ! Get height
      z         = chi(3) - scaled_radius
      ! Translated longitude
      long_dash = chi(1) - 2.0_r_def * pi * time / end_time
      ! Set divergent correction velocity
      u_d = 2.0_r_def * sin(2.0_r_def * pi *( z / domain_height - 0.5_r_def )) *              &
            sin(pi * z / domain_height) - cos(2.0_r_def * pi *( z / domain_height - 0.5_r_def )) &
            * cos(pi * z / domain_height)
      ! Set velocity
      u(1) = scaled_radius * w0 * pi / domain_height                     &
             * sin(long_dash) * sin(long_dash)                        &
             * cos(pi * time / end_time) * cos(chi(2))  * cos(chi(2)) &
             * u_d + 2.0_r_def * pi * scaled_radius / end_time * cos(chi(2))
      u(2) = 0.0_r_def
      u(3) = w0 * sin(2.0_r_def * long_dash)                       &
                * cos(2.0_r_def * pi *( z / domain_height - 0.5_r_def )) &
                * sin(pi * z / domain_height) * cos(pi * time / end_time) * cos(chi(2))
    case default
      write( log_scratch_space, '(A)' )  'Invalid velocity profile choice, stopping'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end select

end function analytic_wind

end module analytic_wind_profiles_mod
