!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

!> @brief Functions to compute the analytic profiles of the deep atmosphere
!>        baroclinic wave of Ullrich et.al. QJRMS2013
module deep_baroclinic_wave_mod

  use constants_mod,                only : r_def, i_def, pi
  use planet_config_mod,            only : scaled_radius, gravity, Rd, cp, &
                                           scaled_omega, kappa, p_zero
  use formulation_config_mod,       only : shallow
  use initial_wind_config_mod,      only : profile, &
                                           profile_deep_baroclinic_perturbed
  implicit none

  private

  public :: deep_baroclinic_wave, T0E, T0P, B, KK, lapse, pertu0,         &
            pertr, pertlon, pertlat, pertz, dxepsilon
  ! evaluate_streamfunction() is a public subroutine in order to facilitate
  ! unit testing.
  public :: evaluate_streamfunction

!=======================================================================
! Deep baroclinic wave parameters
!=======================================================================
  real(kind=r_def), parameter ::     &
       T0E        = 310.0_r_def,     & ! temperature at equatorial surface (K)
       T0P        = 240.0_r_def,     & ! temperature at polar surface (K)
       B          = 2.0_r_def ,      & ! jet half-width parameter
       KK         = 3.0_r_def ,      & ! jet width parameter
       lapse      = 0.005_r_def        ! lapse rate parameter

  real(kind=r_def), parameter ::           &
       pertu0     = 0.5_r_def,             & ! Perturbation wind velocity (m/s)
       pertr      = 1.0_r_def/6.0_r_def,         & ! Perturbation radius (Earth radii)
       pertlon    = pi/9.0_r_def,          & ! Perturbation longitude
       pertlat    = 2.0_r_def*pi/9.0_r_def,& ! Perturbation latitude
       pertz      = 15000.0_r_def   ,      & ! Perturbation height cap
       dxepsilon  = 1.0e-5_r_def             ! Small value for numerical derivatives


contains

  subroutine deep_baroclinic_wave(lon, lat, z,       &
                                  exner, theta, rho, &
                                  u, v, w, mr_v)

    implicit none

!-----------------------------------------------------------------------
!     input/output params parameters at given location
!-----------------------------------------------------------------------
    real(kind=r_def), intent(in)  :: &
                lon,        & ! Longitude (radians)
                lat,        & ! Latitude (radians)
                z             ! Altitude (m)


    real(kind=r_def), intent(out) :: &
                exner,      & ! Exner pressure
                theta,      & ! Potential Temperature (K)
                rho,        & ! density (kg m^-3)
                u,          & ! Zonal wind (m s^-1)
                v,          & ! Meridional wind (m s^-1)
                w,          & ! Vertical Velocity (m s^-1)
                mr_v          ! Water vapour mixing ratio (kg kg^-1)


    real(kind=r_def) :: p, t, q, eta
    real(kind=r_def) :: T0, constA, constB, constC, constH, scaledZ
    real(kind=r_def) :: tau1, tau2, inttau1, inttau2
    real(kind=r_def) :: rratio, inttermT, inttermU, bigU, rcoslat, omegarcoslat

!-----------------------------------------------------------------------
!    constants / parameters
!-----------------------------------------------------------------------

    T0      = 0.5_r_def * (T0E + T0P)
    constA  = 1.0_r_def / lapse
    constB  = (T0 - T0P) / (T0 * T0P)
    constC  = 0.5_r_def * (KK + 2.0_r_def) * (T0E - T0P) / (T0E * T0P)
    constH  = Rd * T0 / gravity
    scaledZ = z / (B * constH)

!-----------------------------------------------------------------------
!    tau values
!-----------------------------------------------------------------------
    tau1 = constA * lapse / T0 * exp(lapse * z / T0) &
         + constB * (1.0_r_def - 2.0_r_def * scaledZ**2) * exp(- scaledZ**2)
    tau2 = constC * (1.0_r_def - 2.0_r_def * scaledZ**2) * exp(- scaledZ**2)

    inttau1 = constA * (exp(lapse * z / T0) - 1.0_r_def) &
            + constB * z * exp(- scaledZ**2)
    inttau2 = constC * z * exp(- scaledZ**2)

!-----------------------------------------------------------------------
!    radius ratio
!-----------------------------------------------------------------------
    if ( shallow ) then
      rratio = 1.0_r_def
    else
      rratio = (z + scaled_radius) / scaled_radius;
    end if

!-----------------------------------------------------------------------
!    interior term on temperature expression
!-----------------------------------------------------------------------
    inttermT = (rratio * cos(lat))**int(KK,i_def) &
             - KK / (KK + 2.0_r_def) &
             * (rratio * cos(lat))**int(KK + 2.0_r_def,i_def)

!-----------------------------------------------------------------------
!    temperature
!-----------------------------------------------------------------------
    t = 1.0_r_def / (rratio**2 * (tau1 - tau2 * inttermT))

!-----------------------------------------------------------------------
!    hydrostatic pressure
!-----------------------------------------------------------------------
    p = p_zero * exp(- gravity / Rd * (inttau1 - inttau2 * inttermT))

!-----------------------------------------------------------------------
!    density (via ideal gas law)
!-----------------------------------------------------------------------
    rho = p / (Rd * t)
!-----------------------------------------------------------------------
!    Compute exner pressure
!-----------------------------------------------------------------------
    exner = (p/p_zero)**kappa

!-----------------------------------------------------------------------
!    Compute potential temperaure
!-----------------------------------------------------------------------
    theta = t/exner

!-----------------------------------------------------------------------
!    velocity field
!-----------------------------------------------------------------------
    inttermU = (rratio * cos(lat))**int(KK - 1.0_r_def,i_def) &
             - (rratio * cos(lat))**int(KK + 1.0_r_def,i_def)
    bigU = gravity / scaled_radius * KK * inttau2 * inttermU * t
    if ( shallow ) then
      rcoslat = scaled_radius * cos(lat)
    else
      rcoslat = (z + scaled_radius) * cos(lat)
    end if

    omegarcoslat = scaled_omega * rcoslat

    u = - omegarcoslat + sqrt(omegarcoslat**2 + rcoslat * bigU)
    v = 0.0_r_def
    w = 0.0_r_def

!-----------------------------------------------------------------------
!    perturbation on velocity field (streamfunction type)
!-----------------------------------------------------------------------
    if ( profile == profile_deep_baroclinic_perturbed ) then
      u = u - 1.0_r_def / (2.0_r_def * dxepsilon) *             &
          ( evaluate_streamfunction(lon, lat + dxepsilon, z)    &
          - evaluate_streamfunction(lon, lat - dxepsilon, z))

      v = v + 1.0_r_def / (2.0_r_def * dxepsilon * cos(lat)) *  &
          ( evaluate_streamfunction(lon + dxepsilon, lat, z)    &
          - evaluate_streamfunction(lon - dxepsilon, lat, z))
    end if

    eta = p/p_zero

    if ( eta >  1.0e4_r_def/p_zero) then
      q = 1.8e-2_r_def*exp(-(lat/(2.0_r_def*pi/9.0_r_def))**4.0_r_def) &
          *exp(-((eta-1.0_r_def)*p_zero/340.0e2_r_def)**2.0_r_def )
    else
      q = 1.0e-12_r_def
    end if

    mr_v = q/(1.0_r_def - q)


  end subroutine deep_baroclinic_wave

!>@brief stream function perturbation function
  real(kind=r_def) function evaluate_streamfunction(lon, lat, z)

    implicit none

    real(kind=r_def), intent(in)  :: &
                lon,        & ! Longitude (radians)
                lat,        & ! Latitude (radians)
                z             ! Altitude (meters)

    real(kind=r_def) :: greatcircler, perttaper, cospert

    ! Great circle distance
    greatcircler = 1.0_r_def / pertr &
      * acos(sin(pertlat) * sin(lat) &
           + cos(pertlat) * cos(lat) * cos(lon - pertlon))

    ! Vertical tapering of stream function
    if (z < pertz) then
      perttaper = 1.0_r_def - 3.0_r_def * z**2 / pertz**2 &
                            + 2.0_r_def * z**3 / pertz**3
    else
      perttaper = 0.0_r_def
    end if

    ! Horizontal tapering of stream function
    if (greatcircler <= 1.0_r_def) then
      cospert = cos(0.5_r_def * pi * greatcircler)
    else
      cospert = 0.0_r_def
    end if

    evaluate_streamfunction = &
        (- pertu0 * pertr * perttaper * cospert**4)

  end function evaluate_streamfunction

end module deep_baroclinic_wave_mod
