!----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief LFRic interface module for UM code (planet_constants_mod)
!----------------------------------------------------------------------------

module planet_constants_mod

  use, intrinsic :: iso_fortran_env, only: real32
  ! Universal constants
  use constants_mod, only: l_def, i_def, r_def, i_um, r_um, pi, rmdi, imdi, r_bl
  use driver_water_constants_mod, only: gas_constant_h2o
  use conversions_mod, only: rsec_per_day

  implicit none

  private
  public :: set_planet_constants
  public :: c_virtual, cp, cv, etar, g, grcp, kappa, lcrcp, lfrcp, ls, lsrcp, &
            one_minus_epsilon, one_minus_epsilon_32b, p_zero, planet_radius,  &
            pref, r, recip_a2, recip_kappa, repsilon, repsilon_32b, rv, vkman,&
            recip_epsilon, omega, two_omega, s2r, lapse, lcrcp_bl, lsrcp_bl,  &
            g_bl, grcp_bl, vkman_bl, pref_bl, kappa_bl, cp_bl, rd_bl,         &
            c_virtual_bl, etar_bl, repsilon_bl, ls_bl, r_32b, c_virtual_32b,  &
            etar_32b, lcrcp_32b, ls_32b, lsrcp_32b, planet_radius_bl,         &
            recip_kappa_bl, power, ex_power, lcrcp_def, lsrcp_def,            &
            recip_kappa_def

  ! The following variables have been hidden as they are not currently
  ! required to build the extracted UM code. They have been left in
  ! in case they are required as more UM code is drawn into the lfric_atm
  ! build. Should they be required at a later date, they should simply be
  ! added to the public statement above.

  ! Disabled variables:
  !   sclht, omega, two_omega, recip_p_zero,
  !   g_over_r


!----------------------------------------------------------------------
! Primary planet constants
!----------------------------------------------------------------------

  ! Planet radius in metres
  real(r_um), protected :: planet_radius = real(rmdi, r_um)

  ! Mean acceleration due to gravity at the planet surface
  real(r_um), protected :: g = real(rmdi, r_um)

  ! Gas constant for dry air
  real(r_um), protected :: r = real(rmdi, r_um)

  ! Specific heat of dry air at constant pressure
  real(r_um), protected :: cp = real(rmdi, r_um)

  ! Reference surface pressure
  real(r_um), protected :: pref = real(rmdi, r_um)

  ! Mean scale height for pressure
  real(r_um), protected :: sclht = real(rmdi, r_um)

  ! Near surface environmental lapse rate
  real(r_um), protected :: lapse = real(rmdi, r_um)

  ! Powers used in PMSL and pressure level calculations
  real(r_um), protected :: power = real(rmdi, r_um)
  real(r_def), protected :: ex_power = real(rmdi, r_um)

  ! Angular speed of planet rotation
  real(r_um), protected :: omega = real(rmdi, r_um)

  ! Gas constant for water vapour
  real(r_um), parameter :: rv = real(gas_constant_h2o, r_um )

  ! Increment to Earth's hour angle per day number from epoch
  real(r_um), parameter :: earth_dha = 2.0_r_um*pi

  ! Von Karman's constant
  real(r_um), parameter :: vkman = 0.4_r_um
  real(r_bl), parameter :: vkman_bl = real(vkman, r_bl)

!----------------------------------------------------------------------
! Derived planet constants
!----------------------------------------------------------------------

  ! Angular speed of planet rotation x2
  real(r_um), protected :: two_omega

  ! Seconds-to-radians converter
  real(r_um), protected :: s2r                 ! planet_dha/rsec_per_day

  ! Ratio of molecular weights of water and dry air
  real(r_um), protected :: repsilon            ! r/rv

  real(r_um), protected :: p_zero              ! pref
  real(r_um), protected :: recip_p_zero        ! 1.0/pref
  real(r_um), protected :: kappa               ! r/cp
  real(r_um), protected :: recip_kappa         ! 1.0/kappa
  real(r_um), protected :: recip_epsilon       ! 1.0/repsilon
  real(r_um), protected :: c_virtual           ! 1.0/repsilon-1.0
  real(r_um), protected :: one_minus_epsilon   ! 1.0-repsilon
  real(r_um), protected :: etar                ! 1.0/(1.0-repsilon)
  real(r_um), protected :: grcp                ! g/cp
  real(r_um), protected :: lcrcp               ! lc/cp
  real(r_um), protected :: lfrcp               ! lf/cp
  real(r_um), protected :: ls                  ! lc+lf
  real(r_um), protected :: lsrcp               ! (lc+lf)/cp
  real(r_um), protected :: cv                  ! cp-r
  real(r_um), protected :: recip_a2            ! 1.0/(planet_radius*planet_radius)
  real(r_um), protected :: g_over_r            ! g/r

  ! 32-bit versions of variables
  real(real32), protected :: repsilon_32b
  real(real32), protected :: one_minus_epsilon_32b
  real(real32), protected :: r_32b
  real(real32), protected :: c_virtual_32b
  real(real32), protected :: etar_32b
  real(real32), protected :: lcrcp_32b
  real(real32), protected :: ls_32b
  real(real32), protected :: lsrcp_32b

  ! BL precision versions
  real(r_bl), protected :: lcrcp_bl               ! lc/cp
  real(r_bl), protected :: lsrcp_bl               ! (lc+lf)/cp
  real(r_bl), protected :: g_bl
  real(r_bl), protected :: grcp_bl                ! g/cp
  real(r_bl), protected :: pref_bl
  real(r_bl), protected :: kappa_bl               ! r/cp
  real(r_bl), protected :: recip_kappa_bl
  real(r_bl), protected :: cp_bl
  real(r_bl), protected :: c_virtual_bl           ! 1.0/repsilon-1.0
  real(r_bl), protected :: etar_bl                ! 1.0/(1.0-repsilon)
  real(r_bl), protected :: rd_bl
  real(r_bl), protected :: repsilon_bl
  real(r_bl), protected :: ls_bl                  ! lc+lf
  real(r_bl), protected :: planet_radius_bl

  ! Versions at native r_def precision
  real(r_def), protected :: lcrcp_def
  real(r_def), protected :: lsrcp_def
  real(r_def), protected :: recip_kappa_def

contains

subroutine set_planet_constants()

  use extrusion_config_mod, only: radius => planet_radius
  use planet_config_mod, only: gravity, rd, epsilon,                 &
                               lfric_recip_epsilon => recip_epsilon, &
                               lfric_omega => omega,                 &
                               lfric_cp => cp,                       &
                               lfric_p_zero => p_zero

  use orbit_config_mod, only: spin, spin_user, hour_angle_inc

  use driver_water_constants_mod, only: latent_heat_h2o_condensation, &
                                        latent_heat_h2o_fusion,       &
                                        gas_constant_h2o

  implicit none

  omega         = real(lfric_omega, r_um)
  r             = real(rd, r_um)
  cp            = real(lfric_cp, r_um)
  g             = real(gravity, r_um)
  planet_radius = real(radius, r_um)
  pref          = real(lfric_p_zero, r_um)
  repsilon      = real(epsilon, r_um)
  recip_epsilon = real(lfric_recip_epsilon, r_um)
  kappa         = real(rd/lfric_cp, r_um)

  ! These variables left in hardwired to earth values as LFRic does not
  ! currently read in any data for these variables.
  sclht          = 6.8e+03_r_um
  lapse          = 0.0065_r_um
  power          = g / (r * lapse)
  ex_power       = (cp * lapse) / g

  ! Set derived constants
  two_omega         = 2.0_r_um * omega

  select case (spin)
  case (spin_user)
    s2r = hour_angle_inc/rsec_per_day
  case default
    s2r = earth_dha/rsec_per_day
  end select

  p_zero            = pref
  recip_p_zero      = 1.0_r_um / pref
  recip_kappa       = 1.0_r_um / kappa
  one_minus_epsilon = real(1.0_r_def - epsilon, r_um )
  c_virtual         = real(1.0_r_def/epsilon - 1.0_r_def, r_um )

  etar     = real( 1.0_r_def/(1.0_r_def - epsilon), r_um )
  grcp     = real( gravity/lfric_cp, r_um )
  lcrcp    = real( latent_heat_h2o_condensation / lfric_cp, r_um )
  lfrcp    = real( latent_heat_h2o_fusion / lfric_cp, r_um )
  ls       = real( latent_heat_h2o_condensation + latent_heat_h2o_fusion, r_um )
  lsrcp    = real( (latent_heat_h2o_condensation + latent_heat_h2o_fusion) / &
                   lfric_cp, r_um )

  cv       = real( lfric_cp - rd, r_um )

  recip_a2 = real( 1.0_r_def/(radius**2_i_def), r_um )
  g_over_r = real( gravity/rd, r_um )


  ! Set 32-bit versions as required, eg in qsat_mod
  repsilon_32b          = real( epsilon, real32 )
  one_minus_epsilon_32b = real( 1.0_r_def - epsilon, real32 )
  r_32b = real(r, real32)
  c_virtual_32b = real(c_virtual, real32)
  etar_32b = real(etar, real32)
  lcrcp_32b = real(lcrcp, real32)
  ls_32b = real(ls, real32)
  lsrcp_32b = real(lsrcp, real32)

  ! Set BL precision versions
  lcrcp_bl = real(lcrcp, r_bl)
  lsrcp_bl = real(lsrcp, r_bl)
  g_bl = real(g, r_bl)
  grcp_bl = real(grcp, r_bl)
  pref_bl = real(pref, r_bl)
  kappa_bl = real(kappa, r_bl)
  recip_kappa_bl = real(recip_kappa, r_bl)
  cp_bl = real(cp, r_bl)
  rd_bl = real(r, r_bl)
  c_virtual_bl = real(c_virtual, r_bl)
  etar_bl = real(etar, r_bl)
  repsilon_bl = real(repsilon, r_bl)
  ls_bl = real(ls, r_bl)
  planet_radius_bl = real(planet_radius, r_bl)

  ! Set BL precision versions
  lcrcp_def = real(lcrcp, r_def)
  lsrcp_def = real(lsrcp, r_def)
  recip_kappa_def = real(recip_kappa, r_def)
  
end subroutine set_planet_constants

end module planet_constants_mod
