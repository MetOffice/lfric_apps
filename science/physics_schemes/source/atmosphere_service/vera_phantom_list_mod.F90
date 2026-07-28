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
! This module contains definitions for casting phantom aerosol populations.
!

module vera_phantom_list_mod

use vera_kind_mod,          only: wp => vera_real, wi => vera_integer

use vera_global_mod,        only: vera_constants          ,                    &
                                  vera_configs            ,                    &
                                  vera_visbty             ,                    &
                                  vera_mie                ,                    &
                                  vera_phantom            ,                    &
                                  vera_kohler             ,                    &
                                  vera_scaled             ,                    &
                                  vera_noise_control      ,                    &
                                  vera_population_type    ,                    &
                                  vera_aerosol_population ,                    &
                                  vera_update_population

use vera_phantom_cast_mod,  only: vera_cast_rd_log_normal_b0_triangle,         &
                                  vera_cast_rd_log_normal_b0         ,         &
                                  vera_volume_scale

implicit none

! Description:
!   This module contains definitions for casting phantom aerosol populations.
!
! Method:
!
!   This module comprises just the single subroutine:
!
!   vera_phantom_list
!     This subroutine consists of a series of executable blocks of code
!     that contain a definition of a phantom aerosol species and any
!     associated Vera settings, e.g. whether or not to employ the Mie
!     scattering scheme to compute the scattering coefficient of the
!     hydrated aerosol particles.
!
!     The form of these blocks is:
!
!     if ( vera_config == vera_configs % choice_of_configuration ) then
!       .
!       .
!       .
!       code defining the aerosol population and Vera settings
!       .
!       .
!       .
!     end if
!
!   There are some pre-defined configurations, that are used for testing:
!
!   vera_config ==  vera_configs % visbty_emulation
!     Replicates the old VISBTY scheme used in the UM since 1992. The MURK
!     aerosol is scaled using a power law to produce monodisperse aerosol.
!     Simple geometric scattering is used, i.e. Qext=2, together with a
!     shape parameter, eta=0.75, to give an effective extinction coefficient
!     Qeff = eta . Q_ext = 0.75 . 2 = 1.5
!
!   vera_config == vera_configs % visbty_emulation_mie_scattering
!     Monodisperse aerosol scaled off the MURK field using the VISBTY
!     power law. Scattering is computed using simple Mie scattering.
!
!   vera_config == log_normal_small_particles_300_mie
!     Log-normal size distribution with 100 particle sizes, scaled off the
!     MURK aerosol field, but biased towards small particles. Each particle
!     size has 3 hygroscopies giving a total of 300 particle species.
!     Scattering is computed using simple Mie scattering.
!
!   vera_config == log_normal_1018_rd_mie_scattering
!     Log-normal size distribution with 1018 particle sizes, scaled off the
!     MURK aerosol field. Each particle has identical hygroscopy.
!     Scattering is computed using simple Mie scattering.
!
!   vera_config == vera_configs % log_normal_generic_mie_scattering
!     Definition for an aerosol population with log-normal size
!     distribution and triangular hygroscopy distribution. The scaling of
!     both the particle sizes and hygroscopy are determined by the
!     variables stored in the object vera_phantom in the Vera
!     global data module vera_global_module.F90. The variables in the
!     object vera_phantom can be assigned during runtime, so this
!     configuration is very flexible.
!
!     The log-normal size distribution isscaled off the MURK field using
!     the VISBTY power law scaling.
!
!     The hygroscopy assumes a triangular b0 distribution.
!
!     If there are n_b0 hygroscopy values and n_rd size values, then
!     the aerosol population consists of n_b0.n_rd species of particles,
!     i.e. there is an instance of each hygroscopy for each particle
!     size.
!
!     Mie scattering is used to compute the scattering coefficient.
!
!
!   For more detail, please refer to the Vera user guide.
!
! Code description:
!   Language: Fortran 2003
!   This code is written to UMDP3 standards.

! name of this module
character (len=*), parameter, private :: ModuleName='VERA_PHANTOM_LIST_MOD'

private

! make all the subroutines in this module Public,
! they are all called by vera_scheme_mod.F90
public :: vera_phantom_list

contains

  !=============================================================================
  !
  ! vera_phantom_list
  !
  ! Casts phantom aerosol populations, selected by the input vera_config, e.g.
  !
  ! vera_config == 1  monodisperse aerosol, scaled off MURK using the VISBTY
  !                   scheme power laws, and with geometric scattering,
  !                   i.e. Qext=2
  !
  !=============================================================================

subroutine vera_phantom_list( vera_config, initialise_vera_parameters )

  ! use the DrHook stuff
use parkind1,   only: jpim, jprb
use yomhook,    only: lhook, dr_hook

implicit none

!===========================================================================
! arguments for vera_phantom_list
!
!===========================================================================

! configuration number of the phantom aerosol population to cast
integer (wi), intent(in)  :: vera_config

! Switch to determine whether we initialise parameters.
! We should only need to do this on the first call to this routine.
logical, intent(in) :: initialise_vera_parameters

!===========================================================================
! local variables for vera_phantom_list
!===========================================================================

! the aerosol population to cast - remember to deallocate this array before
!                                  returning from this routine
type (vera_population_type), allocatable :: cast_population(:)

! number of rd values to use
integer (wi)              :: n_rd

! maxmimum value of rd
real    (wp)              :: rd_max

! minmimum value of rd
real    (wp)              :: rd_min

! mode of the rd distribution
real    (wp)              :: rd_mode

! geometric standard deviation of the rd distribution,
! i.e. the width of the distribution
real    (wp)              :: rd_sigma

! number of b0 values to use
integer (wi)              :: n_b0

! values of b0 in the distribution
real    (wp), allocatable :: b0(:)

! maxmimum value of b0 in the distribution
real    (wp)              :: b0_max

! minmimum value of b0 in the distribution
real    (wp)              :: b0_min

! b0 value corresponding to the peak of the triangular b0 distribution
! NOTE: x_peak must be in the range ]b0_min, b0_max[
real    (wp)              :: b0_peak

! DrHook variables
integer (kind=jpim), parameter :: zhook_in  = 0 ! DrHook tracing entry
integer (kind=jpim), parameter :: zhook_out = 1 ! DrHook tracing exit
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName = 'VERA_PHANTOM_LIST'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!===========================================================================
! executable code for vera_phantom_list
!===========================================================================

if ( vera_config == vera_configs % visbty_emulation_flexible ) then
  ! This aerosol definition replicates the old VISBTY scheme,
  ! i.e. monodisperse aerosol scaled off the MURK field and
  ! a scattering scheme using an effective extinction efficiency computed
  ! from eta the shape parameter and the large particle limit for geometric
  ! scattering, Qext=2, i.e.
  !
  ! Qeff = eta . Qext = 0.75 x 2 = 1.5

  if ( initialise_vera_parameters ) then
!$OMP CRITICAL(vera_phantom_list_switches)
    ! switch on the VISBTY scheme to perturb the total water q_tot,
    vera_visbty%switch_q_total = vera_visbty%switch_q_total_on

    ! switch on the VISBTY scattering scheme,
    vera_mie%visbty_scattering = vera_mie%visbty_scattering_on
!$OMP end CRITICAL(vera_phantom_list_switches)
  end if

  ! just the one aerosol species
  allocate( cast_population(vera_constants%one_i) )

  ! simple scaled monodisperse aerosol population, scaled off the MURK field
  cast_population%nc = vera_scaled%nc
  cast_population%rd = vera_scaled%rd
  cast_population%b0 = vera_kohler%visbty_b0

  ! update the aerosol definition in the Vera global data module
  call vera_update_population( cast_population )

  ! tidy up - get rid of the temporary definition of the
  ! cast aerosol population
  if ( allocated( cast_population ) ) deallocate( cast_population )

  ! end of vera_config = vera_configs % visbty_emulation_flexible
end if

!===========================================================================

if ( vera_config == vera_configs % visbty_emulation_mie_scattering ) then
  ! Monodisperse aerosol, based on the VISBTY scaling from the MURK field,
  ! but using Mie scattering.

  if ( initialise_vera_parameters ) then
!$OMP CRITICAL(vera_phantom_list_switches)
    ! switch on the VISBTY scheme to perturb the total water q_tot,
    vera_visbty%switch_q_total = vera_visbty%switch_q_total_on

    ! ensure that Mie scattering is used, so switch off both geometric
    ! scattering and the old VISBTY scattering scheme
    vera_mie%geometric_scattering = vera_mie%geometric_scattering_off
    vera_mie%visbty_scattering    = vera_mie%visbty_scattering_off
!$OMP end CRITICAL(vera_phantom_list_switches)
  end if

  ! just the one aerosol species
  allocate( cast_population(vera_constants%one_i) )

  ! simple scaled monodisperse aerosol population, scaled off the MURK field
  cast_population%nc = vera_scaled%nc
  cast_population%rd = vera_scaled%rd
  cast_population%b0 = vera_kohler%visbty_b0

  ! update the aerosol definition in the Vera global data module
  call vera_update_population( cast_population )

  ! tidy up - get rid of the temporary definition of the
  ! cast aerosol population
  if ( allocated( cast_population ) ) deallocate( cast_population )

  ! end of vera_config = vera_configs % visbty_emulation_mie_scattering
end if

!===========================================================================

if ( vera_config == vera_configs % log_normal_small_particles_300_mie ) then
  !
  ! MURK scale particles, broad range of rd, biased towards small particles.
  !
  ! Log-normal rd distribution scaled off the VISBTY MURK scaling
  ! with mode(rd) = MURK rd / 2, broad gaussian rd_sigma = 5.0
  ! all of distribution is small particles < mode rd
  !
  ! 100 particle sizes, just 3 hygroscopy values
  !

  if ( initialise_vera_parameters ) then
!$OMP CRITICAL(vera_phantom_list_switches)
    ! switch off the VISBTY scheme to perturb the total water q_tot,
    vera_visbty%switch_q_total = vera_visbty%switch_q_total_off

    ! ensure that Mie scattering is used, so switch off both geometric
    ! scattering and the old VISBTY scattering scheme
    vera_mie%geometric_scattering = vera_mie%geometric_scattering_off
    vera_mie%visbty_scattering    = vera_mie%visbty_scattering_off
!$OMP end CRITICAL(vera_phantom_list_switches)
  end if

  ! specify the particle size distribution
  n_rd     = 100_wi
  rd_mode  = vera_scaled%rd / 1.50_wp
  rd_max   = rd_mode * 2.0_wp
  rd_min   = rd_mode / 2.0_wp
  rd_sigma = 5.0_wp

  ! specify the hygroscopy values
  n_b0     = 3_wi
  allocate( b0(n_b0) )
  b0       = [ 0.05_wp, 0.14_wp, 0.50_wp ]

  ! cast the phantom aerosol population
  ! NOTE - this will populate the aerosol population
  !        vera_aerosol_population in the Vera global data module
  call vera_cast_rd_log_normal_b0( n_rd, rd_max, rd_min ,                      &
                                   rd_mode, rd_sigma    ,                      &
                                   b0                     )

  ! tidy up
  if ( allocated( b0 ) ) deallocate( b0 )

  ! scale the aerosol population to retain the same aerosol
  ! total volume as MURK
  call vera_volume_scale( vera_aerosol_population )

  ! end of vera_config = vera_configs % log_normal_small_particles_300_mie
end if

!===========================================================================

if ( vera_config == vera_configs % log_normal_1018_rd_mie_scattering ) then
  !
  ! MURK scale particles, broad range of rd.
  !
  ! Log-normal rd distribution scaled off the VISBTY MURK scaling,
  ! broad gaussian rd_sigma = 5.0
  !
  ! 1018 particle sizes, just 1 hygroscopy value - this is tha maximum
  !                                                number of particle
  !                                                species Vera can cope
  !                                                with. Any more gives
  !                                                a memory fault.
  !

  if ( initialise_vera_parameters ) then
!$OMP CRITICAL(vera_phantom_list_switches)
    ! switch off the VISBTY scheme to perturb the total water q_tot,
    vera_visbty%switch_q_total = vera_visbty%switch_q_total_off

    ! ensure that Mie scattering is used, so switch off both geometric
    ! scattering and the old VISBTY scattering scheme
    vera_mie%geometric_scattering = vera_mie%geometric_scattering_off
    vera_mie%visbty_scattering    = vera_mie%visbty_scattering_off
!$OMP end CRITICAL(vera_phantom_list_switches)
  end if

  ! specify the particle size distribution
  n_rd     = 1018_wi
  rd_mode  = vera_scaled%rd * vera_phantom%rd_mode_scale
  rd_max   = rd_mode * 2.0_wp
  rd_min   = rd_mode / 2.0_wp
  rd_sigma = 5.0_wp

  ! specify the hygroscopy value
  n_b0     = 1_wi
  allocate( b0(n_b0) )
  b0       = [ 0.14_wp ]

  ! cast the phantom aerosol population
  ! NOTE - this will populate the aerosol population
  !        vera_aerosol_population in the Vera global data module
  call vera_cast_rd_log_normal_b0( n_rd, rd_max, rd_min ,                      &
                                   rd_mode, rd_sigma    ,                      &
                                   b0                     )

  ! tidy up
  if ( allocated( b0 ) ) deallocate( b0 )

  ! scale the aerosol population to retain the same aerosol
  ! total volume as MURK
  call vera_volume_scale( vera_aerosol_population )

  ! end of vera_config = vera_configs % log_normal_1018_rd_mie_scattering
end if

!===========================================================================

if ( vera_config == vera_configs % log_normal_generic ) then
  !
  ! Definition for an aerosol population with log-normal size
  ! distribution and triangular hygroscopy distribution. The scaling of
  ! both the particle sizes and hygroscopy are scaled from the parameters
  ! stored in the object vera_phantom in the global data module
  ! vera_global_module.F90
  !
  ! log-normal rd distribution scaled off the VISBTY MURK scaling
  !
  ! triangular b0 distribution
  !
  ! this configuration does not set the scattering properties, this
  ! should be done via the namelist RUN_Vera

  ! specify the particle size distribution
  n_rd     = vera_phantom%n_rd
  rd_mode  = vera_scaled%rd * vera_phantom%rd_mode_scale
  rd_max   = rd_mode        * vera_phantom%rd_max
  rd_min   = rd_mode        * vera_phantom%rd_min
  rd_sigma = vera_phantom%rd_sigma

  ! specify the hygroscopy distribution
  n_b0     = vera_phantom%n_b0
  b0_max   = vera_phantom%b0_max
  b0_min   = vera_phantom%b0_min
  b0_peak  = vera_phantom%b0_peak

  ! cast the aerosol population
  ! NOTE - this will populate the aerosol population
  !        vera_aerosol_population in the Vera global data module
  call vera_cast_rd_log_normal_b0_triangle( n_rd, rd_max, rd_min ,             &
                                            rd_mode, rd_sigma    ,             &
                                            n_b0, b0_max, b0_min ,             &
                                            b0_peak                )

  ! scale the aerosol population to retain the same aerosol
  ! total volume as MURK
  call vera_volume_scale( vera_aerosol_population )

  ! end of vera_config = vera_configs % log_normal_generic
end if

!===========================================================================

if ( vera_config == vera_configs % log_normal_generic_mie_scattering ) then
  !
  ! Definition for an aerosol population with log-normal size
  ! distribution and triangular hygroscopy distribution. The scaling of
  ! both the particle sizes and hygroscopy are scaled from the parameters
  ! stored in the object vera_phantom in the global data module
  ! vera_global_module.F90
  !
  ! log-normal rd distribution scaled off the VISBTY MURK scaling
  !
  ! triangular b0 distribution
  !

  if ( initialise_vera_parameters ) then
!$OMP CRITICAL(vera_phantom_list_switches)
    ! switch off the VISBTY scheme to perturb the total water q_tot,
    vera_visbty%switch_q_total = vera_visbty%switch_q_total_off

    ! ensure that Mie scattering is used, so switch off both geometric
    ! scattering and the old VISBTY scattering scheme
    vera_mie%geometric_scattering = vera_mie%geometric_scattering_off
    vera_mie%visbty_scattering    = vera_mie%visbty_scattering_off
!$OMP end CRITICAL(vera_phantom_list_switches)
  end if

  ! specify the particle size distribution
  n_rd     = vera_phantom%n_rd
  rd_mode  = vera_scaled%rd * vera_phantom%rd_mode_scale
  rd_max   = rd_mode        * vera_phantom%rd_max
  rd_min   = rd_mode        * vera_phantom%rd_min
  rd_sigma = vera_phantom%rd_sigma

  ! specify the hygroscopy distribution
  n_b0     = vera_phantom%n_b0
  b0_max   = vera_phantom%b0_max
  b0_min   = vera_phantom%b0_min
  b0_peak  = vera_phantom%b0_peak

  ! cast the aerosol population
  ! NOTE - this will populate the aerosol population
  !        vera_aerosol_population in the Vera global data module
  call vera_cast_rd_log_normal_b0_triangle( n_rd, rd_max, rd_min ,             &
                                            rd_mode, rd_sigma    ,             &
                                            n_b0, b0_max, b0_min ,             &
                                            b0_peak                )

  ! scale the aerosol population to retain the same aerosol
  ! total volume as MURK
  call vera_volume_scale( vera_aerosol_population )

  ! end of vera_config = vera_configs % log_normal_generic_mie_scattering
end if

!===========================================================================

if ( vera_config == vera_configs % visbty_emulation ) then
  ! This aerosol definition replicates the old VISBTY scheme,
  ! i.e. monodisperse aerosol scaled off the MURK field and
  ! a scattering scheme using an effective extinction efficiency computed
  ! from eta the shape parameter and the large particle limit for geometric
  ! scattering, Qext=2, i.e.
  !
  ! Qeff = eta . Qext = 0.75 x 2 = 1.5
  !
  ! This configuration generates synthetic noise on the T and q inputs.
  !
  ! The output is the median of the computed set of noisy visibilties.

  if ( initialise_vera_parameters ) then
!$OMP CRITICAL(vera_phantom_list_switches)
    ! switch on the VISBTY scheme to perturb the total water q_tot,
    vera_visbty%switch_q_total = vera_visbty%switch_q_total_on

    ! switch on the VISBTY scattering scheme,
    vera_mie%visbty_scattering = vera_mie%visbty_scattering_on

    ! switch off the noise scheme
    vera_noise_control % n_noise = vera_constants % zero_i
!$OMP end CRITICAL(vera_phantom_list_switches)
  end if

  ! just the one aerosol species
  allocate( cast_population(vera_constants%one_i) )

  ! simple scaled monodisperse aerosol population, scaled off the MURK field
  cast_population%nc = vera_scaled%nc
  cast_population%rd = vera_scaled%rd
  cast_population%b0 = vera_kohler%visbty_b0

  ! update the aerosol definition in the Vera global data module
  call vera_update_population( cast_population )

  ! tidy up - get rid of the temporary definition of the
  ! cast aerosol population
  if ( allocated( cast_population ) ) deallocate( cast_population )

  ! end of vera_config = vera_configs % visbty_emulation
end if

!===========================================================================

if ( vera_config == vera_configs % log_normal_expensive ) then
  !
  ! Definition for an aerosol population with log-normal size
  ! distribution and triangular hygroscopy distribution. The scaling of
  ! both the particle sizes and hygroscopy are scaled from the parameters
  ! stored in the object vera_phantom in the global data module
  ! vera_global_module.F90
  !
  ! 128 particle log-normal rd distribution scaled off the
  ! VISBTY MURK scaling
  !
  ! 6 particle triangular b0 distribution
  !

  if ( initialise_vera_parameters ) then
!$OMP CRITICAL(vera_phantom_list_switches)
    ! switch off the VISBTY scheme to perturb the total water q_tot,
    vera_visbty%switch_q_total = vera_visbty%switch_q_total_off

    ! ensure that Mie scattering is used, so switch off both geometric
    ! scattering and the old VISBTY scattering scheme
    vera_mie%geometric_scattering = vera_mie%geometric_scattering_off
    vera_mie%visbty_scattering    = vera_mie%visbty_scattering_off
!$OMP end CRITICAL(vera_phantom_list_switches)
  end if

  ! specify the particle size distribution
  n_rd     = 32_wi
  rd_mode  = vera_scaled%rd * vera_phantom%rd_mode_scale
  rd_max   = rd_mode        * vera_phantom%rd_max
  rd_min   = rd_mode        * vera_phantom%rd_min
  rd_sigma = vera_phantom%rd_sigma

  ! specify the hygroscopy distribution
  n_b0     = 8_wi
  b0_max   = vera_phantom%b0_max
  b0_min   = vera_phantom%b0_min
  b0_peak  = vera_phantom%b0_peak

  ! cast the aerosol population
  ! NOTE - this will populate the aerosol population
  !        vera_aerosol_population in the Vera global data module
  call vera_cast_rd_log_normal_b0_triangle( n_rd, rd_max, rd_min ,             &
                                            rd_mode, rd_sigma    ,             &
                                            n_b0, b0_max, b0_min ,             &
                                            b0_peak                )

  ! scale the aerosol population to retain the same aerosol
  ! total volume as MURK
  call vera_volume_scale( vera_aerosol_population )

  ! end of vera_config = vera_configs % log_normal_expensive
end if

!===========================================================================

if ( vera_config == vera_configs % log_normal_cheap ) then
  !
  ! Definition for an aerosol population with log-normal size
  ! distribution and triangular hygroscopy distribution. The scaling of
  ! both the particle sizes and hygroscopy are scaled from the parameters
  ! stored in the object vera_phantom in the global data module
  ! vera_global_module.F90
  !
  ! 4 particle log-normal rd distribution scaled off the
  ! VISBTY MURK scaling
  !
  ! 2 particle triangular b0 distribution
  !

  if ( initialise_vera_parameters ) then
!$OMP CRITICAL(vera_phantom_list_switches)
    ! switch off the VISBTY scheme to perturb the total water q_tot,
    vera_visbty%switch_q_total = vera_visbty%switch_q_total_off

    ! ensure that Mie scattering is used, so switch off both geometric
    ! scattering and the old VISBTY scattering scheme
    vera_mie%geometric_scattering = vera_mie%geometric_scattering_off
    vera_mie%visbty_scattering    = vera_mie%visbty_scattering_off
!$OMP end CRITICAL(vera_phantom_list_switches)
  end if

  ! specify the particle size distribution
  n_rd     = 4_wi
  rd_mode  = vera_scaled%rd * vera_phantom%rd_mode_scale
  rd_max   = rd_mode        * vera_phantom%rd_max
  rd_min   = rd_mode        * vera_phantom%rd_min
  rd_sigma = vera_phantom%rd_sigma

  ! specify the hygroscopy distribution
  n_b0     = 2_wi
  b0_max   = vera_phantom%b0_max
  b0_min   = vera_phantom%b0_min
  b0_peak  = vera_phantom%b0_peak

  ! cast the aerosol population
  ! NOTE - this will populate the aerosol population
  !        vera_aerosol_population in the Vera global data module
  call vera_cast_rd_log_normal_b0_triangle( n_rd, rd_max, rd_min ,             &
                                            rd_mode, rd_sigma    ,             &
                                            n_b0, b0_max, b0_min ,             &
                                            b0_peak                )

  ! scale the aerosol population to retain the same aerosol
  ! total volume as MURK
  call vera_volume_scale( vera_aerosol_population )

  ! end of vera_config = vera_configs % log_normal_cheap
end if

!===========================================================================

if ( vera_config == vera_configs % log_normal_generic_constant_b0 ) then
  !
  ! MURK scale particles.
  !
  ! Log-normal rd distribution scaled off the VISBTY MURK scaling.
  !
  ! Single value of b0.
  !
  !

  ! specify the particle size distribution
  n_rd     = vera_phantom%n_rd
  rd_mode  = vera_scaled%rd * vera_phantom%rd_mode_scale
  rd_max   = rd_mode        * vera_phantom%rd_max
  rd_min   = rd_mode        * vera_phantom%rd_min
  rd_sigma = vera_phantom%rd_sigma

  ! specify the hygroscopy value
  n_b0     = 1_wi
  allocate( b0(n_b0) )
  b0       = [ 0.14_wp ]

  ! cast the phantom aerosol population
  ! NOTE - this will populate the aerosol population
  !        vera_aerosol_population in the Vera global data module
  call vera_cast_rd_log_normal_b0( n_rd, rd_max, rd_min ,                      &
                                   rd_mode, rd_sigma    ,                      &
                                   b0                     )

  ! tidy up
  if ( allocated( b0 ) ) deallocate( b0 )

  ! scale the aerosol population to retain the same aerosol
  ! total volume as MURK
  call vera_volume_scale( vera_aerosol_population )

  ! end of vera_config = 500
end if

!===========================================================================

! tidy up - get rid of the temporary definition of the
! cast aerosol population
if ( allocated( cast_population ) ) deallocate( cast_population )

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine vera_phantom_list

end module vera_phantom_list_mod
