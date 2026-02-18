!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

!> @brief Functions to compute a variety of analytic profiles
!> @details Collection of functions to return the value of a field at a given
!!          point based upon a specified analytic formula
module analytic_temperature_profiles_mod

use constants_mod,                only : r_def, pi, i_def, l_def
use log_mod,                      only : log_event,                &
                                         log_scratch_space,        &
                                         LOG_LEVEL_ERROR
use coord_transform_mod,          only : xyz2llr, central_angle
use idealised_config_mod,         only : test_cold_bubble_x,           &
                                         test_cold_bubble_y,           &
                                         test_gravity_wave,            &
                                         test_warm_bubble,             &
                                         test_warm_bubble_3d,          &
                                         test_solid_body_rotation,     &
                                         test_solid_body_rotation_alt, &
                                         test_deep_baroclinic_wave,    &
                                         test_cos_phi,                 &
                                         test_cosine_bubble,           &
                                         test_bryan_fritsch,           &
                                         test_grabowski_clark,         &
                                         test_squall_line,             &
                                         test_supercell,               &
                                         test_horizontal_mountain,     &
                                         test_breaking_gw
use initial_density_config_mod,    only : r1, x1, y1, r2, x2, y2,      &
                                          density_max, density_background
use initial_pressure_config_mod,   only : surface_pressure
use base_mesh_config_mod,          only : geometry,                    &
                                          geometry_spherical
use planet_config_mod,             only : p_zero, Rd, kappa, scaled_radius, &
                                          scaled_omega, gravity, cp
use reference_profile_mod,         only : reference_profile
use generate_global_gw_fields_mod, only : generate_global_gw_pert
use initial_wind_config_mod,       only : u0, sbr_angle_lat
use deep_baroclinic_wave_mod,      only : deep_baroclinic_wave
use breaking_gravity_wave_mod,     only : breaking_gravity_wave
use formulation_config_mod,        only : shallow
use extrusion_config_mod,          only : domain_height

implicit none

private

public :: analytic_temperature
public :: analytic_theta_pert
public :: integrate_theta_v
public :: integrate_exner_surf

contains

!> @brief Compute an analytic temperature field
!> @param[in] chi Position in physical coordinates
!> @param[in] choice Integer defining which specified formula to use
!> @result temperature The result temperature field
function analytic_temperature(chi, choice, surface_height) result(temperature)

  implicit none
  real(kind=r_def), intent(in) :: chi(3)
  integer,          intent(in) :: choice
  real(kind=r_def), intent(in) :: surface_height
  real(kind=r_def)             :: temperature

  real(kind=r_def)             :: l, dt
  real(kind=r_def)             :: theta0
  real(kind=r_def), parameter  :: XC     = 0.0_r_def
  real(kind=r_def), parameter  :: YC     = 0.0_r_def
  real(kind=r_def), parameter  :: A      = 5000.0_r_def
  real(kind=r_def), parameter  :: H      = 10000.0_r_def
  real(kind=r_def), parameter  :: XR = 4000.0_r_def, &
                                  ZC_cold = 3000.0_r_def, &
                                  ZC_hot = 260.0_r_def, &
                                  ZC_3d  = 350.0_r_def, &
                                  ZR = 2000.0_r_def
  real(kind=r_def)             :: long, lat, radius, z
  real(kind=r_def)             :: l1, l2
  real(kind=r_def)             :: pressure, density, mr_v, exner
  real(kind=r_def)             :: s, u00, f_sb, t0
  real(kind=r_def)             :: r_on_a
  real(kind=r_def)             :: u, v, w
  real(kind=r_def)             :: theta_trop, T_trop, theta_eq
  real(kind=r_def)             :: z_trop
  real(kind=r_def)             :: p_surf, Phi_s, Nsq
  real(kind=r_def)             :: psp, u0

  if ( geometry == geometry_spherical ) then
    call xyz2llr(chi(1),chi(2),chi(3),long,lat,radius)
    call central_angle(long,lat,x1,y1,l1)
    call central_angle(long,lat,x2,y2,l2)
  else
    long = chi(1)
    lat  = chi(2)
    l1 = sqrt((long-x1)**2 + (lat-y1)**2)
    l2 = sqrt((long-x2)**2 + (lat-y2)**2)
    z = chi(3)
  end if

  temperature = 0.0_r_def

  select case( choice )

  case ( test_gravity_wave )
    ! Set default value
    call reference_profile(pressure, density, temperature, chi, choice)
    if ( geometry == geometry_spherical ) then
      temperature = temperature &
                  +  generate_global_gw_pert(long,lat,radius-scaled_radius)
    else
      theta0 = 0.01_r_def
      temperature = temperature + THETA0 * sin ( PI * chi(3) / H ) &
                            / ( 1.0_r_def + ( chi(1) - XC )**2/A**2 )
    end if

  case ( test_cold_bubble_x )
    ! Set default value
    call reference_profile(pressure, density, temperature, chi, choice)
    l = sqrt( ((chi(1)-XC)/XR)**2 + ((chi(3)-ZC_cold)/ZR)**2 )
    if ( l <= 1.0_r_def ) then
      dt =  15.0_r_def/2.0_r_def*(cos(PI*l)+1.0_r_def)
      temperature = temperature - dt/pressure
    end if

  case ( test_cold_bubble_y )
    ! Set default value
    call reference_profile(pressure, density, temperature, chi, choice)
    l = sqrt( ((chi(2)-XC)/XR)**2 + ((chi(3)-ZC_cold)/ZR)**2 )
    if ( l <= 1.0_r_def ) then
      dt =  15.0_r_def/2.0_r_def*(cos(PI*l)+1.0_r_def)
      temperature = temperature - dt/pressure
    end if

  case( test_warm_bubble )
    ! Set default value
    call reference_profile(pressure, density, temperature, chi, choice)
    l = sqrt( ((chi(1)-XC))**2 + ((chi(3)-ZC_hot))**2 )
    if ( l <= 50.0_r_def ) then
      dt = 0.5_r_def
    else
      dt = 0.5_r_def*exp(-(l-50.0_r_def)**2/(100.0_r_def)**2)
    end if
    temperature = temperature + dt

  ! Test from Bryan and Fritsch (2002)
  case( test_bryan_fritsch )

    l = sqrt(((chi(1)-XC)/2000.0_r_def)**2.0_r_def + ((chi(3)-2000.0_r_def)/2000.0_r_def)**2.0_r_def)
    if ( l <= 1.0_r_def ) then
      dt = 2.0_r_def * (cos(PI*l/2.0_r_def))**2.0_r_def
      temperature = 1.0_r_def + dt / 300.0_r_def
    else
      temperature = 1.0_r_def
    end if

  ! Test from Grabowski and Clark (1991)
  case( test_grabowski_clark )

    ! Value of theta at surface with temperature = 283 K
    t0 = 283.0_r_def * (p_zero / surface_pressure) ** kappa
    s = 1.3e-5_r_def  ! stability, in m^{-1}
    temperature = t0*exp(s*chi(3))

  ! Test from DCMIP 2025
  case ( test_squall_line, test_supercell )

    ! Specify potential temperature at equator
    z = MAX(radius - scaled_radius, 0.0_r_def)
    theta_trop = 343.0_r_def
    theta0 = 300.0_r_def
    T_trop = 213.0_r_def
    z_trop = 12000.0_r_def
    if (z < z_trop) then
      theta_eq = theta0 + (theta_trop - theta0)*(z/z_trop)**1.25_r_def
    else
      theta_eq = theta_trop*EXP(gravity*(z-z_trop)/(cp*T_trop))
    end if

    temperature = theta_eq

  case ( test_horizontal_mountain )
    T0 = 288.0_r_def
    Nsq = gravity**2/(cp*T0)
    u0 = 10.0_r_def
    psp = 100000.0_r_def
    Phi_s = gravity*surface_height

    ! Specify surface pressure
    p_surf = psp * EXP(                                                        &
        -scaled_radius*Nsq*u0/(2.0_r_def*gravity**2*kappa)                     &
        * (u0 / scaled_radius + 2.0_r_def*scaled_omega)                        &
        * ((SIN(lat))**2 - 1.0_r_def)                                          &
        - Nsq*Phi_s/(gravity**2*kappa)                                         &
    )

    ! Obtain vertical pressure given isothermal atmosphere
    pressure = p_surf * EXP( -gravity*(radius-scaled_radius) / (Rd*T0) )

    ! Obtain potential temperature
    exner = (pressure / p_zero) ** kappa
    temperature = T0 / exner

  !> Test from Kelly & Giraldo
  case( test_warm_bubble_3d )
    ! Set default value
    call reference_profile(pressure, density, temperature, chi, choice)
    l = sqrt( (chi(1)-XC)**2 + (chi(2)-YC)**2 + (chi(3)-ZC_3d)**2 )

    if ( abs(l) <= 250.0_r_def ) then
      dt = 0.25_r_def*(1.0_r_def + cos(PI*l/250.0_r_def))
    else
      dt = 0.0_r_def
    end if
    temperature = temperature + dt

  case ( test_solid_body_rotation )
    t0   = 280.0_r_def
    s    = (radius / scaled_radius) *                                          &
           ( cos(lat) * cos(sbr_angle_lat * pi) +                              &
             sin(lat) * sin(sbr_angle_lat * pi) * sin(long) )
    u00  = u0 * (u0 + 2.0_r_def * scaled_omega * scaled_radius) / (t0 * Rd)
    f_sb = 0.5_r_def * u00*s**2
    temperature = t0 * exp(gravity * (radius - scaled_radius) / ( cp * t0 ) )  &
                     * exp(-kappa * f_sb)

  case ( test_solid_body_rotation_alt )
    ! See Staniforth & White (2007) (rotated pole version of example of
    ! Section 5.3 with m = 1, A = 0, n therefore arbitrary, Phi0 = 0).
    ! In shallow geometry the r/a factor is replaced by 1, see Section 6 of
    ! Staniforth & White (2007).
    t0 = 280_r_def
    if (shallow) then
      r_on_a = 1.0_r_def
    else
      r_on_a = radius / scaled_radius
    end if
    s    = r_on_a * ( cos(lat) * cos(sbr_angle_lat * pi) +                     &
                      sin(lat) * sin(sbr_angle_lat * pi) * sin(long) )
    u00  = u0 * (u0 + 2.0_r_def * scaled_omega * scaled_radius) / (t0 * Rd)
    ! f_sb is the Q of (69) of Staniforth & White (2007)
    f_sb = 0.5_r_def * u00*s**2
    ! The first exponential factor is the integral on the RHS of the first
    ! equation of (69), the factor r_on_a in the denominator makes the same
    ! form work for both deep and shallow atmospheres.
    ! kappa = R / cp has been used
    ! Note: temperature is potential temperature
    temperature = t0 * exp(  gravity * (radius - scaled_radius)                &
                           / (cp * t0 * r_on_a) )                              &
                     * exp(-kappa * f_sb)
  case( test_deep_baroclinic_wave )
    call deep_baroclinic_wave(long, lat, radius-scaled_radius, &
                              pressure, temperature, density, &
                              u, v, w, mr_v)

  case( test_breaking_gw )
    call breaking_gravity_wave(long, lat, radius-scaled_radius, &
                              pressure, temperature, density, &
                              u, v, w, mr_v)

  case( test_cos_phi )
    temperature = density_max*cos(lat)**4

  case( test_cosine_bubble )
    l1 = sqrt( ((chi(1) - x1)/r1)**2 + ((chi(3) - y1)/r2)**2 )
    if ( l1 < 1.0_r_def ) then
      temperature = density_background + density_max*cos(0.5_r_def*l1*PI)**2
    else
      temperature = density_background
    end if

  case default
    ! Set default value
    call reference_profile(pressure, density, temperature, chi, choice)

  end select

end function analytic_temperature

function analytic_theta_pert(chi, choice) result(theta_pert)

  use extrusion_config_mod, only: planet_radius

  implicit none
  real(kind=r_def),    intent(in) :: chi(3)
  integer(kind=i_def), intent(in) :: choice
  real(kind=r_def)                :: theta_pert
  real(kind=r_def)                :: long, lat, radius, l, z
  real(kind=r_def)                :: lon_p, lat_p, Rtheta, ds, r_p
  real(kind=r_def),     parameter :: lon_c = 2.0_r_def*PI/3.0_r_def
  real(kind=r_def),     parameter :: lat_c = 0.0_r_def
  real(kind=r_def),     parameter :: z_c = 1500.0_r_def
  real(kind=r_def),     parameter :: z_p = 1500.0_r_def
  real(kind=r_def),     parameter :: dtheta = 3.0_r_def
  integer(kind=i_def),  parameter :: num_bubbles = 9

  integer(kind=i_def) :: i

  if ( geometry == geometry_spherical ) then
    call xyz2llr(chi(1), chi(2), chi(3), long, lat, radius)
  else
    long = chi(1)
    lat  = chi(2)
  end if

  select case( choice )

  case ( test_squall_line )

    z = radius - scaled_radius
    theta_pert = 0.0_r_def
    ds = 10000.0_r_def
    r_p = 5000.0_r_def

    do i = 1, num_bubbles
      lon_p = lon_c
      lat_p = lat_c + (real(i, r_def) - 0.5_r_def * real(num_bubbles + 1, r_def)) * ds / scaled_radius
      call central_angle(long, lat, lon_p, lat_p, l)
      Rtheta = SQRT((l * scaled_radius / r_p)**2 + ((z - z_c) / z_p)**2)

      if (Rtheta <= 1.0_r_def) then
        theta_pert = theta_pert + dtheta * (COS(PI/2.0_r_def*Rtheta))**2
      end if
    end do

  case ( test_supercell )

    r_p = 10000.0_r_def

    z = radius - scaled_radius
    theta_pert = 0.0_r_def

    call central_angle(long, lat, lon_c, lat_c, l)
    Rtheta = SQRT((l * scaled_radius / r_p)**2 + ((z - z_c) / z_p)**2)

    if (Rtheta <= 1.0_r_def) then
      theta_pert = theta_pert + dtheta * (COS(PI/2.0_r_def*Rtheta))**2
    end if

  case default
    ! In other cases, relative humidity is set to zero
    theta_pert = 0.0_r_def

  end select

end function analytic_theta_pert


subroutine integrate_theta_v(theta_v, theta_v_prev, dtheta_v_dz, theta_eq,     &
                             dusq_eq_dz, height_wth, lat_points, nl, num_quad_points, test)

  use planet_config_mod, only: cp, gravity

  implicit none

  integer(kind=i_def), intent(in)    :: nl, test
  integer(kind=i_def), intent(in)    :: num_quad_points
  real(kind=r_def),    intent(in)    :: theta_v_prev(nl, num_quad_points)
  real(kind=r_def),    intent(in)    :: theta_eq(nl)
  real(kind=r_def),    intent(in)    :: dtheta_v_dz(nl, num_quad_points)
  real(kind=r_def),    intent(in)    :: height_wth(nl)
  real(kind=r_def),    intent(in)    :: lat_points(num_quad_points+1)
  real(kind=r_def),    intent(inout) :: theta_v(nl, num_quad_points)
  real(kind=r_def),    intent(in)    :: dusq_eq_dz(nl)


  real(kind=r_def)                :: u_eq(nl)
  real(kind=r_def)                :: dusq_eq_dz_2(nl)
  real(kind=r_def)                :: z, zS, dzU, Us, Uc
  real(kind=r_def)                :: A, B, C
  real(kind=r_def)                :: integrand(nl, num_quad_points)
  real(kind=r_def)                :: integral(nl, num_quad_points)
  integer(kind=i_def)             :: k, j

  ! Squall line test case from DCMIP 2025
  if ( test == test_squall_line ) then
    zS = 2500.0_r_def
    dzU = 1000.0_r_def
    Us = 12.0_r_def
    Uc = 5.0_r_def

  ! Supercell test case from DCMIP 2016
  else if ( test == test_supercell ) then
    zS = 5000.0_r_def
    dzU = 1000.0_r_def
    Us = 30.0_r_def
    Uc = 15.0_r_def
  end if

  ! Zonal solid body rotation but with wind shear
  A = 0.25_r_def * (-zs**2 + 2.0_r_def*zS*dzU - dzU**2) / (zS * dzU)
  B = 0.5_r_def * (zS + dzU) / dzU
  C = - 0.25_r_def * (zS / dzU)

  do j = 1, num_quad_points
    do k = 1, nl
      ! Evaluate analytic wind
      z = height_wth(k)
      if (z < zS - dzU + 1.0E-6_r_def) then
        u_eq(k) = Us*z/zS - Uc
        dusq_eq_dz_2(k) = 2.0_r_def*u_eq(k)*Us/zs
      else if (ABS(z - zS) < dzU + 1.0E-6_r_def) then
        u_eq(k) = (A + B*z/zS + C*(z/zS)**2)*Us - Uc
        dusq_eq_dz_2(k) = 2.0_r_def*u_eq(k)*Us*(B/zs + 2.0_r_def*C*z/(zs**2))
      else
        u_eq(k) = Us - Uc
        dusq_eq_dz_2(k) = 0.0_r_def
      end if
      ! Evaluate integrand
      integrand(k,j) = (                                                       &
        0.5_r_def * SIN(2.0_r_def*lat_points(j)) / gravity                     &
        * (u_eq(k)**2 * dtheta_v_dz(k,j) - theta_v_prev(k,j) * dusq_eq_dz_2(k))  &
      )
    end do
  end do

  ! Integrate the integrand over the latitude points
  integral(:,:) = 0.0_r_def
  do j = 2, num_quad_points
    do k = 1, nl
      ! Use trapezoidal rule for integration
      integral(k,j) = integral(k,j-1)                                          &
        + 0.5_r_def * (integrand(k,j) + integrand(k,j-1))                      &
        * (lat_points(j) - lat_points(j-1))
    end do
  end do

  do j = 1, num_quad_points
    theta_v(:,j) = theta_eq(:) + integral(:,j)
  end do

end subroutine integrate_theta_v

subroutine integrate_exner_surf(exner_wt, theta_v, exner_eq,     &
                                height_wth, lat_points, nl, num_quad_points, test)

  use planet_config_mod, only: cp, gravity

  implicit none

  integer(kind=i_def), intent(in)    :: nl, test
  integer(kind=i_def), intent(in)    :: num_quad_points
  real(kind=r_def),    intent(in)    :: exner_eq(nl)
  real(kind=r_def),    intent(in)    :: lat_points(num_quad_points+1)
  real(kind=r_def),    intent(in)    :: theta_v(nl, num_quad_points)
  real(kind=r_def),    intent(in)    :: height_wth(nl)
  real(kind=r_def),    intent(inout) :: exner_wt(nl, num_quad_points)

  real(kind=r_def)                :: u_eq(nl)
  real(kind=r_def)                :: z, zS, dzU, Us, Uc
  real(kind=r_def)                :: A, B, C
  real(kind=r_def)                :: integrand(nl, num_quad_points)
  real(kind=r_def)                :: integral(nl, num_quad_points)
  integer(kind=i_def)             :: j, k

  ! Squall line test case from DCMIP 2025
  if ( test == test_squall_line ) then
    zS = 2500.0_r_def
    dzU = 1000.0_r_def
    Us = 12.0_r_def
    Uc = 5.0_r_def

  ! Supercell test case from DCMIP 2016
  else if ( test == test_supercell ) then
    zS = 5000.0_r_def
    dzU = 1000.0_r_def
    Us = 30.0_r_def
    Uc = 15.0_r_def
  end if

  ! Zonal solid body rotation but with wind shear
  A = 0.25_r_def * (-zs**2 + 2.0_r_def*zS*dzU - dzU**2) / (zS * dzU)
  B = 0.5_r_def * (zS + dzU) / dzU
  C = - 0.25_r_def * (zS / dzU)

  do j = 1, num_quad_points
    do k = 1, nl
      ! Evaluate analytic wind
      z = height_wth(k)
      if (z < zS - dzU) then
        u_eq(k) = Us*z/zS - Uc
      else if (ABS(z - zS) < dzU) then
        u_eq(k) = (A + B*z/zS + C*(z/zS)**2)*Us - Uc
      else
        u_eq(k) = Us - Uc
      end if

      ! Evaluate integrand
      integrand(k,j) = (                                                       &
        u_eq(k)**2 * SIN(lat_points(j)) * COS(lat_points(j))                   &
        / (cp * theta_v(k,j))                                                  &
      )
    end do
  end do

  ! Integrate the integrand over the latitude points
  integral(:,:) = 0.0_r_def
  do j = 2, num_quad_points
    do k = 1, nl
      ! Use trapezoidal rule for integration
      integral(k,j) = integral(k,j-1)                                          &
        + 0.5_r_def * (integrand(k,j) + integrand(k,j-1))                      &
        * (lat_points(j) - lat_points(j-1))
    end do
  end do

  do j = 1, num_quad_points
    exner_wt(:,j) = exner_eq(:) - integral(:,j)
  end do

end subroutine integrate_exner_surf

end module analytic_temperature_profiles_mod
