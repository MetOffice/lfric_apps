! *********************************COPYRIGHT************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *********************************COPYRIGHT************************************
! Some of the content of this file has been produced with the assistance of
! Anthropic Claude Opus 5 (Claude Code).
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: atmos_service_vera
!
! Part of the Vera scheme for diagnosing visual range.
!
! This module contains the global data used by Vera.
!

module vera_global_mod

! grab the data types to use
use vera_kind_mod, only: wp => vera_real, wi => vera_integer

! grab a value for Pi
use conversions_mod, only: um_pi => pi

! grab water physical data, tm is [K], rho_water is [kg / m3]
use water_constants_mod, only: um_water_freeze => tm        ,                  &
                               um_rho_water    => rho_water

! grab relative molecular mass data - note the UM values are [kg / mol]
use rel_mol_mass_mod, only: um_rmm_air      => rmm_air  ,                      &
                            um_rmm_water    => rmm_w

! grab VISBTY - air density [kg / m3]
!               base particle radius [m]
!               base particle number density [m-3]
!               power law scaling
!               particle hygroscopy for Kohler curve
!               Kelvin term for Kohler curves [m]
use visbty_constants_mod,  only: visbty_rho_air  => rho_air   ,                &
                                 visbty_radius0  => radius0   ,                &
                                 visbty_n0       => n0        ,                &
                                 visbty_power    => power     ,                &
                                 visbty_b0       => b0        ,                &
                                 visbty_a0       => a0


implicit none

private

! Description:
!   This module contains the global data used by Vera. Also included are two
!   subroutines:
!
!     vera_update_population - to update the definition of the
!                              aerosol population.
!
!     vera_update_covariance - to update the definition of the
!                              noise covariance matrix.
!
! Method:
!   Most of this module consists of definitions of global data for the Vera
!   scheme. There are two subroutines:
!
!     vera_update_population
!
!     This subroutine updates the definition of the aerosol
!     population, i.e. (Nc, Rd, B0) where Nc is the number concentration,
!     Rd is the aerosol dry particle radius and B0 is the hygroscopy. Each
!     aerosol species is described by such a triplet.
!
!     vera_update_covariance
!
!     This subroutine updates the global noise covariance matrix using the
!     covariances described in the derived type vera_noise_covariance.
!     For example, the covariance between noise on the input temperature and
!     specific humidity is specified by vera_noise_covariance%cov_t_q and in
!     the (5,5) square noise covariance matrix, this is specific covariance
!     is entry (2,3), i.e.
!
!     vera_noise_covariance%cov_matrix(2,3) = vera_noise_covariance%cov_t_q
!
!   For more detail, please refer to the Vera user guide.
!
! Code description:
!   Language: Fortran 2003
!   This code is written to UMDP3 standards.
! name of this module

character (len=*), parameter, private :: ModuleName='VERA_GLOBAL_MOD'

public  :: vera_update_population    ! this is the subroutine that updates
                                     ! the aerosol population definition

!=============================================================================
! numerical constants
!=============================================================================
! value of zero
real    (wp), parameter, public :: vera_zero   = 0.0_wp

! value of one half
real    (wp), parameter, public :: vera_half   = 0.5_wp

! value of one third
real    (wp), parameter, public :: vera_third  = 1.0_wp / 3.0_wp

! value of one
real    (wp), parameter, public :: vera_one    = 1.0_wp

! value of two
real    (wp), parameter, public :: vera_two    = 2.0_wp

! value of two
real    (wp), parameter, public :: vera_minus_two = -2.0_wp

! value of three
real    (wp), parameter, public :: vera_three   = 3.0_wp

! value of four
real    (wp), parameter, public :: vera_four    = 4.0_wp

! value of six
real    (wp), parameter, public :: vera_six     = 6.0_wp

! value of one hundred
real    (wp), parameter, public :: vera_hundred = 100.0_wp

! value of one thousand
real    (wp), parameter, public :: vera_thousand = 1000.0_wp

! value of one thousand
real    (wp), parameter, public :: vera_nano     = 1.0e-9_wp

! value of Pi, taken from conversions_mod
real    (wp), parameter, public :: vera_pi       = um_pi

! value of four thirds, i.e. 4/3
real    (wp), parameter, public :: vera_four_thirds = vera_four / vera_three

! value of three halves, i.e. 3/2
real    (wp), parameter, public :: vera_three_halves = vera_three / vera_two

! value of four thirds pi, i.e. 4pi/3
real    (wp), parameter, public :: vera_four_thirds_pi =                       &
                                     vera_four_thirds * vera_pi

! value of two pi, i.e. 2pi
real    (wp), parameter, public :: vera_two_pi      = vera_two * vera_pi

! value of four pi, i.e. 4pi
real    (wp), parameter, public :: vera_four_pi     = vera_four * vera_pi

! value of square root of two * pi, i.e. sqrt(2pi)
real    (wp), parameter, public :: vera_root_two_pi = sqrt( vera_two * vera_pi )

! value of -1 as an integer
integer (wi), parameter, public :: vera_minus_one_i    = -1_wi

! value of 0 as an integer
integer (wi), parameter, public :: vera_zero_i         = 0_wi

! value of 1 as an integer
integer (wi), parameter, public :: vera_one_i          = 1_wi

! value of 2 as an integer
integer (wi), parameter, public :: vera_two_i          = 2_wi

! value of 3 as an integer
integer (wi), parameter, public :: vera_three_i        = 3_wi

! value of 100 as an integer
integer (wi), parameter, public :: vera_one_hundred_i  = 100_wi

! value of 1000 as an integer
integer (wi), parameter, public :: vera_one_thousand_i = 1000_wi

! value of the complex number (0 + i)
complex (wp), parameter, public :: vera_complex_i      = ( 0_wp, 1_wp )

! gather the numerical constants into a handy structure
type :: constants
  real    (wp) :: pi             = vera_pi
  real    (wp) :: three_halves   = vera_three_halves
  real    (wp) :: four_thirds    = vera_four_thirds
  real    (wp) :: four_thirds_pi = vera_four_thirds_pi
  real    (wp) :: two_pi         = vera_two_pi
  real    (wp) :: four_pi        = vera_four_pi
  real    (wp) :: root_two_pi    = vera_root_two_pi
  real    (wp) :: zero           = vera_zero
  real    (wp) :: one            = vera_one
  real    (wp) :: two            = vera_two
  real    (wp) :: three          = vera_three
  real    (wp) :: four           = vera_four
  real    (wp) :: six            = vera_six
  real    (wp) :: minus_two      = vera_minus_two
  real    (wp) :: half           = vera_half
  real    (wp) :: third          = vera_third
  real    (wp) :: hundred        = vera_hundred
  real    (wp) :: thousand       = vera_thousand
  real    (wp) :: nano           = vera_nano
  integer (wi) :: minus_one_i    = vera_minus_one_i
  integer (wi) :: zero_i         = vera_zero_i
  integer (wi) :: one_i          = vera_one_i
  integer (wi) :: two_i          = vera_two_i
  integer (wi) :: three_i        = vera_three_i
  integer (wi) :: one_hundred_i  = vera_one_hundred_i
  integer (wi) :: one_thousand_i = vera_one_thousand_i
  complex (wp) :: complex_i      = vera_complex_i
end type constants

! specify the constants as Parameters, and initialise the object
type(constants), public, parameter :: vera_constants = constants(              &
  pi             = vera_pi                                        ,            &
  four_thirds    = vera_four_thirds                               ,            &
  four_thirds_pi = vera_four_thirds_pi                            ,            &
  two_pi         = vera_two_pi                                    ,            &
  four_pi        = vera_four_pi                                   ,            &
  root_two_pi    = vera_root_two_pi                               ,            &
  zero           = vera_zero                                      ,            &
  one            = vera_one                                       ,            &
  two            = vera_two                                       ,            &
  minus_two      = vera_minus_two                                 ,            &
  half           = vera_half                                      ,            &
  third          = vera_third                                     ,            &
  three          = vera_three                                     ,            &
  hundred        = vera_hundred                                   ,            &
  thousand       = vera_thousand                                  ,            &
  nano           = vera_nano                                      ,            &
  minus_one_i    = vera_minus_one_i                               ,            &
  zero_i         = vera_zero_i                                    ,            &
  one_i          = vera_one_i                                     ,            &
  two_i          = vera_two_i                                       )

!=============================================================================
! Vera top level switch
!=============================================================================

! define the states of the switch vera_scheme_flag:
!
! vera_scheme_flag_off - don't use the Vera scheme, this is the default
! vera_scheme_flag_on  - use the Vera scheme
integer (wi), parameter :: vera_scheme_flag_on      = 1_wi
integer (wi), parameter :: vera_scheme_flag_off     = 0_wi
integer (wi), parameter :: vera_scheme_flag_default = vera_scheme_flag_off

! define the Vera configurations used in vera_phantom_list_mod.F90
type :: scheme_flag
  integer (wi) :: vera_scheme_flag     = vera_scheme_flag_default
  integer (wi) :: vera_scheme_flag_on  = vera_scheme_flag_on
  integer (wi) :: vera_scheme_flag_off = vera_scheme_flag_off
end type scheme_flag

type(scheme_flag), public, save :: vera_flag

!=============================================================================
! Vera configurations
!=============================================================================
! define the Vera configurations used in vera_phantom_list_mod.F90
type :: configs
  integer (wi) :: visbty_emulation_flexible          = 1_wi
  integer (wi) :: visbty_emulation_mie_scattering    = 2_wi
  integer (wi) :: log_normal_small_particles_300_mie = 3_wi
  integer (wi) :: log_normal_1018_rd_mie_scattering  = 4_wi
  integer (wi) :: log_normal_generic                 = 9_wi
  integer (wi) :: log_normal_generic_mie_scattering  = 10_wi
  integer (wi) :: visbty_emulation_noisy             = 101_wi
  integer (wi) :: visbty_emulation                   = 102_wi
  integer (wi) :: log_normal_expensive_noisy         = 103_wi
  integer (wi) :: log_normal_expensive               = 104_wi
  integer (wi) :: log_normal_cheap_noisy             = 105_wi
  integer (wi) :: log_normal_cheap                   = 106_wi
  integer (wi) :: log_normal_generic_constant_b0     = 500_wi
end type configs

type(configs), public, save :: vera_configs

!=============================================================================
! properties of water
!=============================================================================
! freezing temperature of water, in Kelvin
real (wp), parameter :: water_freeze = um_water_freeze

!=============================================================================
! densities of aerosol, water and air [kg / m3]
!=============================================================================
real (wp), parameter :: rho_aerosol  = 1700.0_wp
real (wp), parameter :: rho_air      = visbty_rho_air
real (wp), parameter :: rho_water    = um_rho_water

!=============================================================================
! relative moelcular weights, water and dry air [g / mol]
! Note - the UM values are [kg / mol]
!=============================================================================
real (wp), parameter :: mw_water     = um_rmm_water * vera_thousand
real (wp), parameter :: mw_dry_air   = um_rmm_air   * vera_thousand

! gather the physical constants into a handy structure
type :: physics
  real (wp) :: water_freeze = water_freeze
  real (wp) :: rho_aerosol  = rho_aerosol
  real (wp) :: rho_air      = rho_air
  real (wp) :: rho_water    = rho_water
  real (wp) :: mw_water     = mw_water
  real (wp) :: mw_dry_air   = mw_dry_air
end type physics

type(physics), public, save :: vera_physics

!=============================================================================
! water in the scheme, i.e. the saturation water vapour pressure
! and the total water ql + qv from the UM
!=============================================================================
! saturation vapour pressure [Pa]
real (wp), parameter :: vera_e_sat             = vera_zero

! saturation specific humidity [kg / kg]
real (wp), parameter :: vera_q_sat             = vera_zero

! total water mixing ratio in the UM [kg / kg]
real (wp), parameter :: vera_q_total_perturbed = vera_zero

! total water mixing ratio to use in the hydration scheme, [kg / kg]
real (wp), parameter :: vera_q_total_use       = vera_zero

! switch to decide whther to use the UM routines to comute qsat, or the
! equivalent routine in vera_water_mod.F90 - by default, use the UM routines
integer (wi), parameter :: switch_qsat_um_on  = vera_one_i
integer (wi), parameter :: switch_qsat_um_off = vera_zero_i
integer (wi), parameter :: switch_qsat_um     = switch_qsat_um_on

! gather the water variables into a handy structure
type :: water
  real    (wp) :: e_sat              = vera_e_sat
  real    (wp) :: q_sat              = vera_q_sat
  real    (wp) :: q_total_use        = vera_q_total_use
  real    (wp) :: q_total_perturbed  = vera_q_total_perturbed
  integer (wi) :: switch_qsat_um_on  = switch_qsat_um_on
  integer (wi) :: switch_qsat_um_off = switch_qsat_um_off
  integer (wi) :: switch_qsat_um     = switch_qsat_um
end type water

type(water), public, save :: vera_water


!=============================================================================
! aerosol source
!=============================================================================

! aerosol source, possible choices are:
! 1 - background aerosol mass mixing ratio - this is the default
! 2 - MURK field
integer (wi), parameter :: aerosol_source_background = vera_one_i
integer (wi), parameter :: aerosol_source_murk       = vera_two_i
integer (wi), parameter :: default_aerosol_source    = aerosol_source_background

! backgroud aerosol mass mixing ratio, in micrograms / kg
real (wp), parameter :: am_background = 18.956_wp

! gather the MURK aerosol scaling constants into a handy structure
type :: aerosol
  integer (wi) :: aerosol_source             = default_aerosol_source
  integer (wi) :: aerosol_source_background  = aerosol_source_background
  integer (wi) :: aerosol_source_murk        = aerosol_source_murk
  real    (wp) :: am_background              = am_background
end type aerosol

type(aerosol), public, save :: vera_aerosol

!=============================================================================
! power law for scaling the MURK aerosol size and mass mixing ratio
!=============================================================================

! base state of dry monodisperse aerosol -
! particle radius [m], number concentration [m-3] and power law scaling
real (wp), parameter :: radius0       = visbty_radius0
real (wp), parameter :: n0            = visbty_n0
real (wp), parameter :: power         = visbty_power

! the cube of the base radius [m3]
real (wp), parameter :: radius0_cubed = radius0 * radius0 * radius0

! base state aerosol mass mixing ratio, am0  [kg / kg]
real (wp), parameter :: am0 = vera_four_thirds_pi * radius0_cubed *            &
                              ( rho_aerosol / rho_air ) * n0

! gather the MURK aerosol scaling constants into a handy structure
type :: scaling
  real (wp) :: radius0    = radius0
  real (wp) :: n0         = n0
  real (wp) :: power      = power
  real (wp) :: am0        = am0
end type scaling

type(scaling), public, save :: vera_scaling

! define an object to hold the scaled values
type :: scaled
  real (wp) :: rd          = radius0
  real (wp) :: nc          = n0
  real (wp) :: am_over_am0 = am0
  real (wp) :: nc_rd_cubed = radius0_cubed * n0
end type scaled

type(scaled), public, save :: vera_scaled
!$OMP THREADPRIVATE(vera_scaled)

!=============================================================================
! Kohler curve parameters
!=============================================================================
! constant for Kelvin term in the Kohler equation
real (wp), parameter :: a0           = visbty_a0

! default hygroscopy parameter used in Rault's term in VISBTY
real (wp), parameter :: visbty_in_b0 = visbty_b0

! gather the Kohler curve constants into a handy structure
type :: kohler
  real (wp) :: a0        = a0
  real (wp) :: visbty_b0 = visbty_in_b0
end type kohler

type(kohler), public, save :: vera_kohler

!=============================================================================
! VISBTY parameters
!=============================================================================
! default critical Rh, rhcrit
real    (wp), parameter :: rhcrit_default           = 0.92_wp

! default probability of fog, calc_prob_of_fog
real    (wp), parameter :: calc_prob_of_fog_default = 0.40_wp

logical     , parameter :: l_murk_default           = .true.

! default shape parameter, eta, the normalised ratio of the second to
! the third moments of the size distribution, i.e. eta = m2 / m3^(2/3),
! nominally, VISBTY assumes a log-normal size distribution
real    (wp), parameter :: eta            = 0.75_wp

! switch to select the VISBTY scheme for perturbing the total water q_tot
! as a function of rhcrit and prob.
! If this switch is = 1 then the perurbation is switched on, otherwise
! the computed q_total = q + qcl is left unaltered if the switch is set to 0.
!
! The default is for this switch to be set to 0, i.e. the q_tot perturbation
! scheme is turned off. To turn this switch on or off use
!
!  vera_visbty%switch_q_total = vera_visbty%switch_q_total_on
!
!  vera_visbty%switch_q_total = vera_visbty%switch_q_total_off
!
integer (wi), parameter :: switch_q_total_on  = vera_one_i
integer (wi), parameter :: switch_q_total_off = vera_zero_i
integer (wi), parameter :: switch_q_total     = switch_q_total_off

! switch to dtermine if the old UM VISBTY scheme is used rather than Vera:
!   switch_visbty = 0 use Vera
!   switch_visbty = 1 use VISBTY
!
! The default is for this switch to be set to 0, i.e. switched off so that
! the Vera scheme is used. To turn this switch on or off use
!
!  vera_visbty%switch_visbty = vera_visbty%switch_visbty_on
!
!  vera_visbty%switch_visbty = vera_visbty%switch_visbty_off
!
integer (wi), parameter :: switch_visbty_on  = vera_one_i
integer (wi), parameter :: switch_visbty_off = vera_zero_i
integer (wi), parameter :: switch_visbty     = switch_visbty_off

! gather the old VISBTY parameters into a handy structure
type :: visbty
  real    (wp) :: rhcrit             = rhcrit_default
  real    (wp) :: prob               = calc_prob_of_fog_default
  logical      :: l_murk             = l_murk_default
  real    (wp) :: eta                = eta
  integer (wi) :: switch_visbty      = switch_visbty
  integer (wi) :: switch_visbty_on   = switch_visbty_on
  integer (wi) :: switch_visbty_off  = switch_visbty_off
  integer (wi) :: switch_q_total     = switch_q_total
  integer (wi) :: switch_q_total_on  = switch_q_total_on
  integer (wi) :: switch_q_total_off = switch_q_total_off
end type visbty

type(visbty), public, save :: vera_visbty

!=============================================================================
! Koschmeider parameters, for implementing Koschmeider's Law to compute
! the visible range from the scattering coefficients of the aerosol particles
!=============================================================================
! liminal contrast - a parametrisation of the accuity of human vision
real (wp), parameter :: liminal_contrast   = 0.05_wp

! clear air visibility, 100km expressed in metres
real (wp), parameter :: clear_air_vis      = 100.0_wp * 1000.0_wp

! minimum visibility to consider when computing the geometric mean
! visibility from a noisy visibility population [m]
real (wp), parameter :: minimum_visibility = 0.01_wp

! Rayleigh scattering coeffiecient
real (wp), parameter :: beta_rayleigh      = -log(liminal_contrast) /          &
                                              clear_air_vis

! gather the Koschmeider parameters into a handy structure
type :: koschmeider

  real (wp) :: liminal_contrast     = liminal_contrast
  real (wp) :: log_liminal_contrast = -log(liminal_contrast)
  real (wp) :: clear_air_vis        = clear_air_vis

  real (wp) :: minimum_visibility   = minimum_visibility
  real (wp) :: maximum_visibility   = clear_air_vis

  ! maximum aerosol scattering coeffiecient
  real (wp) :: max_beta_aerosol     = ( -log(liminal_contrast) /               &
                                        minimum_visibility ) - beta_rayleigh

  ! minimum aerosol scattering coeffiecient set to zero
  real (wp) :: min_beta_aerosol     = vera_zero

  ! Rayleigh scattering coeffiecient
  real (wp) :: beta_rayleigh       = beta_rayleigh

end type koschmeider

type(koschmeider), public, save :: vera_koschmeider


!=============================================================================
! equation solver parameters
!=============================================================================
! default tolerance to use for the equation solver in MINPACK
real (wp), parameter :: tolerance            = 0.000001_wp

! initial guess to use for the aerosol growth factor
! this must be > 1.0 for Rh > 0
real (wp), parameter :: initial_growth_guess_small =   1.1_wp
real (wp), parameter :: initial_growth_guess_large = 100.0_wp

! growth factor to use if the MINPACK solver fails
real (wp), parameter :: no_growth = 1.0_wp

! flag that specifies which solver to use in MINPACK
!
!    solver_flag = 1   solver uses a numerical approximation to the
!                      Jacobian J(x)
!
!    solver_flag = 2   solver uses an analytical computation of the
!                      Jacobian J(x)
!
! the default is to use solver_flag = 2
!
! for details of the MINPACK solver routines, look at
!
! "User Guide for MINPACK-1"
!  Jorge J. More, Burton S. Garbow, Kenneth E. Hillstrom, 1980
!
integer (wi), parameter :: solver_flag_numeric  = vera_one_i
integer (wi), parameter :: solver_flag_analytic = vera_two_i
integer (wi), parameter :: solver_flag_default  = solver_flag_analytic

! gather the MINPACK solver parameters into a handy structure
type :: minpack
  real    (wp) :: tolerance                  = tolerance
  real    (wp) :: initial_growth_guess_small = initial_growth_guess_small
  real    (wp) :: initial_growth_guess_large = initial_growth_guess_large
  real    (wp) :: initial_growth_guess       = initial_growth_guess_small
  real    (wp) :: no_growth                  = no_growth
  integer (wi) :: solver_flag                = solver_flag_default
  integer (wi) :: solver_flag_numeric        = solver_flag_numeric
  integer (wi) :: solver_flag_analytic       = solver_flag_analytic
end type minpack

type(minpack), public, save :: vera_minpack

!=============================================================================
! phantom aerosol populations parameters
!=============================================================================
! default configuration of Vera
integer (wi), parameter :: vera_default_config       = 10_wi

! switch to determine if a phantom population needs to be cast
! and set the default to ON
integer (wi), parameter :: cast_switch_on            = vera_one_i
integer (wi), parameter :: cast_switch_off           = vera_zero_i
integer (wi), parameter :: cast_switch_default       = cast_switch_on

! default number of distinct aerosol particle sizes in the population
integer (wi), parameter :: n_rd_default              = 16_wi

! default number of distinct aerosol hygroscopies in the population
integer (wi), parameter :: n_b0_default              = 4_wi

! default number of distinct aerosol species in the population
integer (wi), parameter :: n_aerosol_species_default =                         &
                           n_rd_default * n_b0_default

! scaling factor for the mode of the particle size distribution
real    (wp), parameter :: rd_mode_scale             = 1.0_wp

! width of the Rd log-normal distribution as a multiple of the variance
real    (wp), parameter :: rd_spread                 = 1.0_wp

! variance of the Rd log-normal distribution
real    (wp), parameter :: rd_sigma                  = sqrt( 10.0_wp )

! maximum value of the Rd log-normal distribution
real    (wp), parameter :: rd_max                    = rd_sigma**rd_spread

! minimum value of the Rd log-normal distribution
real    (wp), parameter :: rd_min                    = 1.0_wp / rd_max

! maximum value of the B0 distribution
real    (wp), parameter :: b0_max                    = 1.00_wp

! minimum value of the B0 distribution
real    (wp), parameter :: b0_min                    = 0.00_wp

! peak value of a triangular B0 distribution
real    (wp), parameter :: b0_peak                   = 0.14_wp

! tiny nudge to use when computing log-normal distributions,
! used to stretch the log range to avoid rounding error. Value is not
! critical, but must be of the form 1+delta with 0<delta<<1
real (wp), parameter :: delta = 0.1_wp ** ( precision(0.0_wp) - 2_wi )
real (wp), parameter :: nudge = 1.0_wp + delta

! parameters that can be used to construct the phantom population
!
! these parameters can be used to construct a log-normal size
! distribution and a triangular hygroscopy distribution
!
! the range of rd values is given by
! rd = [rd_min, rd_max]*rd_mode_scale*rd_mode
!
! the range of B0 values is given by
! B0 = [b0_min, b0_max]
!
type :: phantom_parameters
  integer (wi) :: vera_config         = vera_default_config
  integer (wi) :: vera_default_config = vera_default_config
  integer (wi) :: n_aerosol_species   = n_aerosol_species_default
  integer (wi) :: n_rd                = n_rd_default
  integer (wi) :: cast_switch         = cast_switch_default
  integer (wi) :: cast_switch_on      = cast_switch_on
  integer (wi) :: cast_switch_off     = cast_switch_off
  real    (wp) :: rd_mode_scale       = rd_mode_scale
  real    (wp) :: rd_sigma            = rd_sigma
  real    (wp) :: rd_spread           = rd_spread
  real    (wp) :: rd_max              = rd_max
  real    (wp) :: rd_min              = rd_min
  integer (wi) :: n_b0                = n_b0_default
  real    (wp) :: b0_max              = b0_max
  real    (wp) :: b0_min              = b0_min
  real    (wp) :: b0_peak             = b0_peak
  real    (wp) :: default_b0          = visbty_b0
  real    (wp) :: nudge               = nudge
end type phantom_parameters

type(phantom_parameters), public, save :: vera_phantom

! A thread private version of vera_phantom%n_aerosol_species
integer (wi), public :: n_aerosol_species_thread
!$OMP THREADPRIVATE(n_aerosol_species_thread)

! define the aerosol population type, popping in default values
type, public :: vera_population_type

  ! aerosol number concentration
  real (wp)  :: nc = n0

  ! dry aerosol radius
  real (wp)  :: rd = radius0

  ! hygroscopy to use in Kohler equation Rault's term
  real (wp)  :: b0 = visbty_b0

end type vera_population_type

! define the polydisperse aerosol population
type(vera_population_type), allocatable, public, save ::                       &
                                                    vera_aerosol_population(:)
!$OMP THREADPRIVATE(vera_aerosol_population)

!=============================================================================
! scattering parameters
!=============================================================================
! large particle upper limit for extinction efficiency Qext
real    (wp), parameter :: qext_upper_limit     = 2.0_wp

! small particle lower limit for extinction efficiency Qext
real    (wp), parameter :: qext_lower_limit     = 0.0_wp

! switch to select geometric scattering, i.e. Mie large particle limit
! with Qext = 2.
! If this switch is = 1 then geometric scattering is switched on, otherwise
! the default Mie scattering scheme is used if the switch is set to 0.
!
! The default is for this switch to be set to 0, i.e. turn off geometric
! scattering. To turn this switch on or off use
!
!  vera_mie%geometric_scattering = vera_mie%geometric_scattering_on
!
!  vera_mie%geometric_scattering = vera_mie%geometric_scattering_off
!
integer (wi), parameter :: geometric_scattering_on  = vera_one_i
integer (wi), parameter :: geometric_scattering_off = vera_zero_i
integer (wi), parameter :: geometric_scattering = geometric_scattering_off

! switch to select the scattering scheme used in VISBTY,
! i.e. an effective extinction efficiency, Qeff = eta . Qext,
! where Qext=2, the Mie large particle limit and eta is the shape parameter
! derived from the assumed log-normal aerosol size distribution,
! with eta = m2 / m3(2/3) the normalised ratio of the second to third moments
! and eta is set to 0.75
!
! So what this means is that Qeff is set to
! Qeff = eta . Qext(large particle limit) = 0.75 x 2 = 1.5
!
! If this switch is = 1 then VISBTY scattering is switched on, otherwise
! the default Mie scattering scheme is used if the switch is set to 0.
!
! The default is for this switch to be set to 0, i.e. turn off the
! VISBTY scattering scheme. To turn this switch on or off use
!
!  vera_mie%visbty_scattering = vera_mie%visbty_scattering_on
!
!  vera_mie%visbty_scattering = vera_mie%visbty_scattering_off
!
integer (wi), parameter :: visbty_scattering_on  = vera_one_i
integer (wi), parameter :: visbty_scattering_off = vera_zero_i
integer (wi), parameter :: visbty_scattering     = visbty_scattering_off

! switch to select using Blumel's approximation to compute Qext for very
! small aerosol particles. If this switch is = 1 then Blumel's approximation
! is used, otherwise very small particles are assigned a value of Qext = 0.
!
! The default is for this switch to be set to 1, i.e. Blumel's approximation
! is turned on. To turn this switch on or off use
!
!  vera_mie%switch_blumel = vera_mie%switch_blumel_on
!
!  vera_mie%switch_blumel = vera_mie%switch_blumel_off
!
integer (wi), parameter :: switch_blumel_on  = vera_one_i
integer (wi), parameter :: switch_blumel_off = vera_zero_i
integer (wi), parameter :: switch_blumel     = switch_blumel_on

! light wavelength to use if computing Qext using Blumel's approximation
real    (wp), parameter :: wavelength        = 550.0_wp * vera_nano

! refractive index to use if computing Qext using Blumel's approximation
complex (wp), parameter :: refractive_index  = (1.53_wp, -0.007_wp)

! gather the mie constants into a handy structure
type :: mie
  real    (wp) :: qext_upper_limit         = qext_upper_limit
  real    (wp) :: qext_lower_limit         = qext_lower_limit
  integer (wi) :: switch_blumel            = switch_blumel
  integer (wi) :: switch_blumel_on         = switch_blumel_on
  integer (wi) :: switch_blumel_off        = switch_blumel_off
  real    (wp) :: wavelength               = wavelength
  complex (wp) :: refractive_index         = refractive_index
  integer (wi) :: geometric_scattering     = geometric_scattering
  integer (wi) :: geometric_scattering_on  = geometric_scattering_on
  integer (wi) :: geometric_scattering_off = geometric_scattering_off
  integer (wi) :: visbty_scattering        = visbty_scattering
  integer (wi) :: visbty_scattering_on     = visbty_scattering_on
  integer (wi) :: visbty_scattering_off    = visbty_scattering_off
end type mie

type(mie), public, save :: vera_mie

!=============================================================================
! synthetic noise settings
!=============================================================================
! default number of instances of synthetic noise to use is set to zero,
! i.e. don't use the synthetic noise scheme
integer   (wi), parameter :: n_noise_default    = vera_zero_i

! define a default number of noisy instances to use when switching on the
! noise scheme from a trial configuration
integer   (wi), parameter :: n_noise_trial      = vera_one_hundred_i

! number of variables in a Vera input set i.e. (P, T, q, qcl, am)
integer   (wi), parameter :: n_inputs           = 5_wi

! switch to indicate no noise is required on the inputs to Vera
integer   (wi), parameter :: noise_switch_off   = vera_zero_i

! switch to indicate whether noise is required on the aerosol
! input to Vera
integer  (wi), parameter :: aerosol_noise_on    = vera_one_i
integer  (wi), parameter :: aerosol_noise_off   = vera_zero_i

! switch to indicate whether to sort a list of noisy visibilities
integer  (wi), parameter :: sort_on             = vera_one_i
integer  (wi), parameter :: sort_off            = vera_zero_i

! maximum number of centiles to compute from a noisy visibility population
integer   (wi), parameter :: n_centiles_max     = 20_wi

! maximum number of range probablities to compute from a noisy
! visibility population
integer   (wi), parameter :: n_ranges_max       = 20_wi

! character length of a Vera input variable
integer   (wi), parameter :: cl                 = 4_wi

! variables Vera input set i.e. (P, T, q, qcl, am)
character (cl), parameter :: vera_inputs(n_inputs) =                           &
                             ['P  ', 'T  ', 'q  ', 'qcl', 'am ']

! definition of the median, i.e. the 50th centile
real    (wp), parameter :: fifty_cent(1)        = 50.0_wp

! gather the noisy output options into a handy structure
type :: outputs
  real    (wp) :: fifty_cent(1)                 = fifty_cent
  real    (wp) :: median_centile(1)             = fifty_cent
end type outputs

type(outputs), public, save :: vera_outputs

! define some numbers to populate the noise covariance matrix
real    (wp), parameter, private   :: zero         = vera_zero
real    (wp), parameter, private   :: one          = vera_one

real    (wp), parameter, private   :: cov_t_q      = 0.80_wp
real    (wp), parameter, private   :: cov_t_qcl    = zero
real    (wp), parameter, private   :: cov_t_am     = zero
real    (wp), parameter, private   :: cov_p_t      = zero
real    (wp), parameter, private   :: cov_p_q      = zero
real    (wp), parameter, private   :: cov_p_qcl    = zero
real    (wp), parameter, private   :: cov_p_am     = zero
real    (wp), parameter, private   :: cov_q_qcl    = zero
real    (wp), parameter, private   :: cov_q_am     = zero
real    (wp), parameter, private   :: cov_qcl_am   = zero

! covariance matrix for the noise on an input set (P, T, q, qcl, am)
! NOTE - Fortran uses column-major order, so in this covariance matrix
!        definition, the five sub-arrays specify the columns, so the matrix
!        specified is as shown here:
!
!              P        T        q        qcl      am
!
!        P   | 1        cov_p_t   cov_p_q   cov_p_qcl cov_p_am   |
!        T   | 0        1         cov_t_q   cov_t_qcl cov_t_am   |
!        q   | 0        0         1         cov_q_qcl cov_q_am   |
!        qcl | 0        0         0         1         cov_q_clam |
!        am  | 0        0         0         0         1          |
!
real     (wp), parameter :: default_cov_matrix5(5, 5) =                        &
  reshape( [one       , zero      , zero      , zero       , zero      ,       &
             cov_p_t   , one       , zero      , zero       , zero      ,      &
             cov_p_q   , cov_t_q   , one       , zero       , zero      ,      &
             cov_p_qcl , cov_t_qcl , cov_q_qcl , one        , zero      ,      &
             cov_p_am  , cov_t_am  , cov_q_am  , cov_qcl_am , one     ],       &
           [n_inputs, n_inputs] )

! switch to determine whether to use correlated noise, i.e. use the
! error covariance matrix
!
! default is to switch this scheme off, i.e. use uncorrelated
! gaussian noise
!
integer (wi), parameter :: switch_covariance_on  = vera_one_i
integer (wi), parameter :: switch_covariance_off = vera_zero_i
integer (wi), parameter :: switch_covariance     = switch_covariance_off


! gather the noise covariance data into a handy structure
type :: noise_covariance

  real     (wp) :: cov_t_q     =   cov_t_q
  real     (wp) :: cov_t_qcl   =   cov_t_qcl
  real     (wp) :: cov_t_am    =   cov_t_am
  real     (wp) :: cov_p_t     =   cov_p_t
  real     (wp) :: cov_p_q     =   cov_p_q
  real     (wp) :: cov_p_qcl   =   cov_p_qcl
  real     (wp) :: cov_p_am    =   cov_p_am
  real     (wp) :: cov_q_qcl   =   cov_q_qcl
  real     (wp) :: cov_q_am    =   cov_q_am
  real     (wp) :: cov_qcl_am  =   cov_qcl_am

  ! switch to turn the covariance matrix scheme on and off, default
  ! is off
  integer  (wi) :: switch_covariance     = switch_covariance
  integer  (wi) :: switch_covariance_on  = switch_covariance_on
  integer  (wi) :: switch_covariance_off = switch_covariance_off

  ! covariance matrix for the noise perturbations this is
  ! contructed with the variables ordered (P, T, q, qcl, am)
  real     (wp) :: cov_matrix(n_inputs, n_inputs) = default_cov_matrix5

  ! number of variables in a Vera input set i.e. (P, T, q, qcl, am)
  integer  (wi) :: n_variables           = n_inputs

  ! variables in a Vera input set i.e. (P, T, q, qcl, am)
  character (cl):: vera_inputs(n_inputs) = vera_inputs

end type noise_covariance

type(noise_covariance), public, save :: vera_noise_covariance

! flag to indicate whether the random number generator seed has been
! planted, this is initialised to zero
!   seed_flag = 0 not set yet, so reset the seed
!   seed_flag > 0 has been set, so leave the seed alone
integer   (wi), parameter :: seed_flag_reset = vera_zero_i
integer   (wi), parameter :: seed_flag_leave = vera_one_i

! seed to use for the random number generator random_number
integer   (wi), parameter :: seed_default    = 19701219_wi

! method used to set the seed, options are:
!
!   seed_method_constant - constant value defined by seed_default
!   seed_method_model    - model run date/time data
!   seed_method_computer - computer clock date/time data computed at runtime
integer   (wi), parameter :: seed_method_constant = vera_zero_i
integer   (wi), parameter :: seed_method_model    = vera_one_i
integer   (wi), parameter :: seed_method_computer = vera_two_i

! Post-Processing list of visible threshold probabilities to compute
real      (wp), parameter :: ranges_pp(13)   =                                 &
 [ 50.0_wp, 100.0_wp, 200.0_wp, 400.0_wp, 600.0_wp, 800.0_wp, 1000.0_wp,       &
   1500.0_wp, 2000.0_wp, 5000.0_wp, 10000.0_wp, 20000.0_wp, 40000.0_wp   ]

! Post-processing list of centiles to compute
real      (wp), parameter :: centiles_pp(5)   =                                &
                            [ 1.0_wp, 10.0_wp, 50.0_wp, 90.0_wp, 99.0_wp ]

! gather the noise control settings into a handy structure
type :: noise_control

  ! switch to indicate whether noise is required on the aerosol
  ! input to Vera - default is no aerosol noise
  integer  (wi) :: aerosol_noise     = aerosol_noise_off
  integer  (wi) :: aerosol_noise_on  = aerosol_noise_on
  integer  (wi) :: aerosol_noise_off = aerosol_noise_off

  ! switch to indicate whether to sort a list of noisy visibilities
  integer  (wi) :: sort_switch  = sort_on
  integer  (wi) :: sort_on      = sort_on
  integer  (wi) :: sort_off     = sort_off

  ! number of variables in a Vera input set i.e. (P, T, q, qcl, am)
  integer  (wi) :: n_variables       = n_inputs

  ! number of noise perturbations to generate, per variable
  integer  (wi) :: n_noise           = n_noise_default

  ! set the default number of noise perturbations to generate
  integer  (wi) :: n_noise_default  = n_noise_default

  ! number of noise perturbations to generate from a trial configuration
  integer  (wi) :: n_noise_trial    = n_noise_trial

  ! switch to indicate if synthetic noise is required
  integer  (wi) :: noise_switch_off  = noise_switch_off

  ! maximum number of centiles to compute from a noisy
  ! visibility population
  integer  (wi) :: n_centiles_max    = 20_wi

  ! number of centiles to compute from a noisy
  ! visibility population
  integer  (wi) :: n_centiles        = 5_wi

  ! maximum number of range probablities to compute from a noisy
  ! visibility population
  integer  (wi) :: n_ranges_max      = 20_wi

  ! number of range probablities to compute from a noisy
  ! visibility population
  integer  (wi) :: n_ranges          = 13_wi

  ! ensemble number to seed the random number generator random_number
  integer  (wi) :: ensemble_member   = 0_wi

  ! default seed to use for the random number generator random_number
  integer  (wi) :: seed_default      = seed_default

  ! seed to use for the random number generator random_number
  integer  (wi) :: seed              = seed_default

  ! flag to indicate whether the random number generator seed has been
  ! planted, this is initialised to zero
  !   seed_flag = 0 not set
  !   seed_flag > 0 set
  integer  (wi) :: seed_flag_reset   = seed_flag_reset
  integer  (wi) :: seed_flag_leave   = seed_flag_leave
  integer  (wi) :: seed_flag         = seed_flag_reset

  ! method used to seed the random number generator random_number
  integer  (wi) :: seed_method_constant = seed_method_constant
  integer  (wi) :: seed_method_model    = seed_method_model
  integer  (wi) :: seed_method_computer = seed_method_computer
  integer  (wi) :: seed_method          = seed_method_model

  ! maximum number of loops to attempt to generate bounded gaussian noise
  integer  (wi) :: max_loops         = 10000000_wi

  ! centiles to compute from a population of noisy visibilities
  real     (wp) :: centiles(n_centiles_max)            = vera_zero
  real     (wp) :: centiles_pp(size(centiles_pp))      = centiles_pp
  real     (wp) :: centiles_default(size(centiles_pp)) = centiles_pp

  ! visible ranges [m] for which to estimate probabilities from a
  ! population of noisy visibilities
  real     (wp) :: ranges(n_ranges_max)            = vera_zero
  real     (wp) :: ranges_fog(1)                   = [ 1000.0_wp ]
  real     (wp) :: ranges_fog_mist(2)              = [ 1000.0_wp, 5000.0_wp ]
  real     (wp) :: ranges_pp(size(ranges_pp))      = ranges_pp
  real     (wp) :: ranges_default(size(ranges_pp)) = ranges_pp

end type noise_control

type(noise_control), public, save :: vera_noise_control

! gather the noise scalings into a handy structure - this defines
! standard deviation and bias, of the gaussian noise
type :: noise_scaling

  ! bounds for un-scaled gaussian noise, the range [lower, upper]sigma
  real     (wp) :: upper_bound_normal =   1.0_wp
  real     (wp) :: lower_bound_normal =  -1.0_wp

  ! sd and bias for pressure noise [Pa]
  ! both of these are absolute values, with the noise scaling as
  ! delta_p = N(gaussian) * width_p + bias_p
  real     (wp) :: width_p            =   vera_zero
  real     (wp) :: bias_p             =   vera_zero

  ! sd and bias for temperature noise [K]
  ! both of these are absolute values, with the noise scaling as
  ! delta_t = N(gaussian) * width_t + bias_t
  real     (wp) :: width_t            =   0.1_wp
  real     (wp) :: bias_t             =   vera_zero

  ! sd and bias for specific humidity noise,
  ! the bias is an absolute value in [kg/kg],
  ! whereas the width the a relative value, i.e. width_q=0.01
  ! will set sigma_q to be 0.01 of the unperturbed q value, i.e.
  ! the noise scales as:
  ! delta_q = N(gaussian) * width_q * q + bias_q
  real     (wp) :: width_q            =   0.01_wp
  real     (wp) :: bias_q             =   vera_zero

  ! sd and bias for liquid water noise,
  ! the bias is an absolute value in [kg/kg],
  ! whereas the width the a relative value, i.e. width_qcl=0.01
  ! will set sigma_qcl to be 0.01 of the unperturbed qcl value, i.e.
  ! the noise scales as:
  ! delta_qcl = N(gaussian) * width_qcl * qcl + bias_qcl
  real     (wp) :: width_qcl          =   vera_zero
  real     (wp) :: bias_qcl           =   vera_zero

  ! sd and bias for aerosol mass mixing ratio noise,
  ! the bias is an absolute value in [microgram/kg],
  ! whereas the width the a relative value, i.e. width_am=0.01
  ! will set sigma_q to be 0.01 of the unperturbed am value, i.e.
  ! the noise scales as:
  ! delta_am = N(gaussian) * width_am * am + bias_am
  real     (wp) :: width_am           =   vera_zero
  real     (wp) :: bias_am            =   vera_zero

end type noise_scaling

type(noise_scaling), public, save :: vera_ns

! A workspace for MINPACK-solve related variables
type :: thread_minpack
  real    (wp) :: q_sat              = vera_q_sat
  real    (wp) :: q_total_use        = vera_q_total_use
end type thread_minpack

type(thread_minpack), public, save :: vera_thread_minpack
!$OMP THREADPRIVATE(vera_thread_minpack)

contains

  !=============================================================================
  !
  ! vera_update_population - Provides a method to update the
  ! global aerosol popluation definition, ie (Nc, Rd, B0).
  !
  !=============================================================================

subroutine vera_update_population( population_in )

  ! use the DrHook stuff
use parkind1,   only: jpim, jprb
use yomhook,    only: lhook, dr_hook

implicit none

!===========================================================================
! argument for vera_update_population
!===========================================================================

! new definition for the global aerosol population
type(vera_population_type), intent(in) :: population_in(:)

!===========================================================================
! local variables for vera_update_population
!===========================================================================

! number of species in new aerosol population
integer (wi) :: n_species
! and in the old aerosol population
integer (wi) :: n_species_old

! DrHook variables
integer (kind=jpim), parameter :: zhook_in  = 0 ! DrHook tracing entry
integer (kind=jpim), parameter :: zhook_out = 1 ! DrHook tracing exit
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName = 'VERA_UPDATE_POPULATION'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!===========================================================================
! start of the executable code for vera_update_population
!===========================================================================

n_species = size( population_in )
! if the aerosol population structure already exists, then we will check
! whether it is the right size or whether we should get rid of it.
if ( allocated( vera_aerosol_population ) ) then
  ! Must not call size on an unallocated array
  n_species_old = size( vera_aerosol_population )
  if (n_species /= n_species_old) then
    deallocate( vera_aerosol_population )
  end if
end if
! set up a new structure for the aerosol population definition if required
if (.not. allocated( vera_aerosol_population ) ) then
  allocate( vera_aerosol_population( n_species ) )
end if
! Update n_aerosol_species_thread
n_aerosol_species_thread = n_species

! fill in the global aerosol population structure
vera_aerosol_population = population_in

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine vera_update_population

!=============================================================================
!
! vera_update_covariance - Provides a method to update the
! global noise covariance matrix.
!
!=============================================================================

subroutine vera_update_covariance( )

  ! use the DrHook stuff
use parkind1,   only: jpim, jprb
use yomhook,    only: lhook, dr_hook

implicit none

!===========================================================================
! argument for vera_update_covariance
!===========================================================================

! this routine does not have any arguments

!===========================================================================
! local variables for vera_update_covariance
!===========================================================================

! DrHook variables
integer (kind=jpim), parameter :: zhook_in  = 0 ! DrHook tracing entry
integer (kind=jpim), parameter :: zhook_out = 1 ! DrHook tracing exit
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName = 'VERA_UPDATE_COVARIANCE'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!===========================================================================
! start of the executable code for vera_update_covariance
!===========================================================================

! update the covariance matrix in the global data module

! set covariance(P,T)
vera_noise_covariance%cov_matrix(1,2) = vera_noise_covariance%cov_p_t

! set covariance(P,q)
vera_noise_covariance%cov_matrix(1,3) = vera_noise_covariance%cov_p_q

! set covariance(P,qcl)
vera_noise_covariance%cov_matrix(1,4) = vera_noise_covariance%cov_p_qcl

! set covariance(P,am)
vera_noise_covariance%cov_matrix(1,5) = vera_noise_covariance%cov_p_am

! set covariance(T,q)
vera_noise_covariance%cov_matrix(2,3) = vera_noise_covariance%cov_t_q

! set covariance(T,qcl)
vera_noise_covariance%cov_matrix(2,4) = vera_noise_covariance%cov_t_qcl

! set covariance(T,am)
vera_noise_covariance%cov_matrix(2,5) = vera_noise_covariance%cov_t_am

! set covariance(q,qcl)
vera_noise_covariance%cov_matrix(3,4) = vera_noise_covariance%cov_q_qcl

! set covariance(q,am)
vera_noise_covariance%cov_matrix(3,5) = vera_noise_covariance%cov_q_am

! set covariance(qcl,am)
vera_noise_covariance%cov_matrix(4,5) = vera_noise_covariance%cov_qcl_am

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine vera_update_covariance

end module vera_global_mod
