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
! Part of the Vera scheme for diagnosing visual range
!
! This is the main work-horse for the Vera visibility scheme.
!
! This module implements the Vera visibility scheme for sets of inputs of
! (P, T, q, qcl, am).
!

module vera_scheme_mod

  ! specify the precision to use
use vera_kind_mod, only: wp => vera_real, wi => vera_integer

implicit none

! Description:
!   This module implements the Vera scheme for sets of inputs of
!   (P, T, q, qcl, am).
!
! Method:
!   This module comprises just the the single subroutine:
!
!   vera_scheme
!     Scales the input MURK aerosol mass mixing ratio, then casts
!     a phantom dry aerosol field. This aerosol is then hydrated
!     using Kohler curves. The scattering coefficient of the hydrated
!     aerosol is computed, which us used to estimate the visibility via
!     Koschmeider's Law.
!
!   For more detail, please refer to the Vera user guide.
!
! Code description:
!   Language: Fortran 2003
!   This code is written to UMDP3 standards.

! name of this module
character (len=*), parameter, private :: ModuleName='VERA_SCHEME_MOD'

private

! make the Vera scheme Public, together with the routine to cast aerosol
! populations from the MURK aerosol mass mixing ratio
public :: vera_murk_cast,                                                      &
          vera_scheme

contains

subroutine vera_murk_cast( am, vera_config )

! grab modules required for vera_murk_cast
use vera_murk_scaling_mod, only: vera_murk_scaling

use vera_phantom_list_mod, only: vera_phantom_list

! use the DrHook stuff
use parkind1,   only: jpim, jprb
use yomhook,    only: lhook, dr_hook

implicit none

!===========================================================================
! arguments for vera_murk_cast
!===========================================================================

! aerosol mass mixing ratio, [micrograms / kilograms]
real    (wp), intent(in) :: am

! Vera configuration to use
integer (wi), intent(in) :: vera_config

!===========================================================================
! local variables for vera_murk_cast
!===========================================================================
! Switch to determine whether we initialise parameters in vera_phantom_list().
! We should only need to do this on the first call to vera_phantom_list().
logical, parameter :: initialise_vera_parameters_off = .false.

! DrHook variables
integer (kind=jpim), parameter :: zhook_in  = 0 ! DrHook tracing entry
integer (kind=jpim), parameter :: zhook_out = 1 ! DrHook tracing exit
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName = 'VERA_MURK_CAST'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!===========================================================================
! scale the aerosol using the MURK power law scaling
!===========================================================================

! scale the input MURK aerosol mass mixing ratio using a power law, and
! compute the monodisperse (Nc, rd), scaled off the base state (Nc0, rd0)
call vera_murk_scaling( am )

! cast the phantom aerosol population
call vera_phantom_list( vera_config, initialise_vera_parameters_off )

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine vera_murk_cast

subroutine vera_scheme( p, t, q, qcl, am ,                                     &
                        vera_config      ,                                     &
                        scattering_ls    ,                                     &
                        scattering_c     ,                                     &
                        fractional_ls    ,                                     &
                        fractional_c     ,                                     &
                        vis_no_precip    ,                                     &
                        vis_precip         )

! grab modules required for Vera
use vera_global_mod,        only: vera_constants          ,                    &
                                  vera_visbty             ,                    &
                                  vera_koschmeider        ,                    &
                                  vera_phantom            ,                    &
                                  vera_minpack            ,                    &
                                  vera_aerosol_population ,                    &
                                  vera_thread_minpack     ,                    &
                                  n_aerosol_species_thread

use vera_water_mod,         only: vera_temperature_pressure_to_q_sat,          &
                                  vera_perturb_q_total

use vera_hydration_mod,     only: vera_hydrate

use vera_function_mod,      only: vera_kohler_function

use vera_mie_mod,           only: vera_mie_lookup


! use the DrHook stuff
use parkind1,   only: jpim, jprb
use yomhook,    only: lhook, dr_hook

implicit none

!===========================================================================
! arguments for vera_scheme
!
! NOTE: the atmospheric inputs (P, T, q, qcl, am) are all vectors
!       i.e. 1-d arrays.
!
!       The output visibility is also a 1-d array, with the number of
!       elements matching the number of input sets of (P, T, q, qcl, am).
!
!===========================================================================

! atmospheric pressure [Pa]
real    (wp),           intent(in)     :: p(:)

! atmospheric temperature [K]
real    (wp),           intent(in)     :: t(:)

! atmospheric specific humidity [kg/kg]
real    (wp),           intent(in)     :: q(:)

! atmospheric liquid water [kg/kg]
real    (wp),           intent(in)     :: qcl(:)

! aerosol mass mixing ratio [micrograms/kg]
real    (wp),           intent(in)     :: am(:)

! population configuration to use
integer (wi), optional, intent(in)     :: vera_config

! large scale precipitation scattering coefficients used to compute total
! visibility in precipitation
real    (wp), optional, intent (in)    :: scattering_ls(:)

! convective precipitation scattering coefficients used to compute total
! visibility in precipitation
real    (wp), optional, intent (in)    :: scattering_c(:)

! total large scale fractional cloud coverage used to compute total
! visibility in precipitation
real    (wp), optional, intent (in)    :: fractional_ls(:)

! convective fractional cloud coverage used to compute total
! visibility in precipitation
real    (wp), optional, intent (in)    :: fractional_c(:)

! computed visibility [m] using just aerosol and Rayleigh scattering,
! i.e. no contribution from precipitation
real    (wp), optional, intent(in out) :: vis_no_precip(:)

! computed visibility [m] including the contribution from precipitation
real    (wp), optional, intent(in out) :: vis_precip(:)

!===========================================================================
! local variables for vera_scheme
!===========================================================================

! configuration to use - could be the default value in the global data
! module, or could use a value input to this routine
integer (wi)                        :: vera_config_use

! number of input points to the Vera scheme, i.e. the size of
! the set (P, T, q, qcl, am)
integer (wi)                        :: n_input_points

! counter to loop over the input (P, T, q, qcl, am)
integer (wi)                        :: ii_point

! saturation specific humidity [kg/kg]
real    (wp)                        :: q_sat

! total atmospheric specific humidity [kg/kg], i.e. q + qcl
real    (wp)                        :: q_total

! perturbed atmospheric specific humidity [kg/kg],
! i.e. f(q + qcl, prob, rhcrit)
real    (wp)                        :: q_total_perturbed

! large scale precipitation scattering coefficients used to compute total
! visibility in precipitation
real    (wp)                        :: precip_scattering_ls(size(p))

! convective precipitation scattering coefficients used to compute total
! visibility in precipitation
real    (wp)                        :: precip_scattering_c(size(p))

! fractional large scale cloud cover used to compute total
! visibility in precipitation
real    (wp)                        :: precip_fractional_ls(size(p))

! fractional large scale cloud cover modified to remove the convective
! cloud cover
real    (wp)                        :: precip_fractional_ls_no_c(size(p))

! fractional convective cloud cover used to compute total
! visibility in precipitation
real    (wp)                        :: precip_fractional_c(size(p))

! growth factor of the aerosol particles, i.e. hydrated radius / dry radius
real    (wp), allocatable           :: growth_factor(:)

! extinction efficiency Qext, computed using simple Mie scattering
real    (wp), allocatable           :: mie_qext(:)

! scattering coefficients of the aerosol species
real    (wp), allocatable           :: beta_aerosol(:)

! total aerosol scattering coefficient - sum of the scattering  from each
! aerosol species
real    (wp)                        :: beta_aerosol_total

! visibility excluding precipitation
real    (wp)                        :: vis_excluding_precip(size(p))

! visibility including scattering from large scale precipitation
real    (wp)                        :: vis_precip_ls(size(p))

! visibility including scattering from convective precipitation
real    (wp)                        :: vis_precip_c(size(p))

! DrHook variables
integer (kind=jpim), parameter :: zhook_in  = 0 ! DrHook tracing entry
integer (kind=jpim), parameter :: zhook_out = 1 ! DrHook tracing exit
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName = 'VERA_SCHEME'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!===========================================================================
! start of the executable code for vera_scheme
!===========================================================================

! how many input points are there?
n_input_points = size(p)

! if a configuration was passed to this routine, then use that
if (present(vera_config)) then
  ! use the input configuration
  vera_config_use = vera_config
else
  ! use the default configuration for Vera
  vera_config_use = vera_phantom%vera_default_config
end if

! check for large scale precipitation scattering coefficients
if (present(scattering_ls)) then
  ! use the input precipitation scattering coefficient
  precip_scattering_ls = scattering_ls
else
  ! set the precipitation scattering coefficient to zero
  precip_scattering_ls = vera_constants%zero
end if

! check for convetive precipitation scattering coefficients
if (present(scattering_c)) then
  ! use the input precipitation scattering coefficient
  precip_scattering_c = scattering_c
else
  ! set the precipitation scattering coefficient to zero
  precip_scattering_c = vera_constants%zero
end if

! check for fractional large scale cloud coverage
if (present(fractional_ls)) then
  ! use the input fractional large scale cloud coverage
  precip_fractional_ls = fractional_ls
else
  ! set the large scale cloud coverage to zero
  precip_fractional_ls = vera_constants%zero
end if

! check for fractional convective cloud coverage
if (present(fractional_c)) then
  ! use the input fractional convective cloud coverage
  precip_fractional_c = fractional_c
else
  ! set the convective cloud coverage to zero
  precip_fractional_c = vera_constants%zero
end if

! compute the modified large scale cloud cover by removing the
! contribution of convective cloud
precip_fractional_ls_no_c = ( vera_constants % one - precip_fractional_c )*    &
                            precip_fractional_ls

! loop over the inputs, ie the input sets of ( P, T, q, qcl, am )
do ii_point = vera_constants%one_i, n_input_points

  !=========================================================================
  ! cast an aerosol population from the MURK aerosol mass mising ratio
  !=========================================================================

  if ( vera_phantom % cast_switch == vera_phantom % cast_switch_on ) then
    call vera_murk_cast( am(ii_point), vera_config_use )
  end if

  !=========================================================================
  ! compute total water content and saturation vapour pressure, esat(P, T)
  !=========================================================================

  ! compute the total atmospheric water content from the input
  ! specific humidity and liquid water
  q_total = q(ii_point) + qcl(ii_point)

  ! compute the saturation specific humidity to use in the hydration scheme
  call vera_temperature_pressure_to_q_sat( t(ii_point), p(ii_point), q_sat )

  ! set the global saturation specific humidity
  vera_thread_minpack%q_sat = q_sat

  ! set the global total q
  vera_thread_minpack%q_total_use = q_total

  ! specify the total q to use in the hydration scheme
  ! i.e. whether to apply the perturbation driven by RHcrit and Prob_Fog
  ! as used in VISBTY, or to ignore this scheme. The switch
  ! switch_q_total controls the use of this pertubation scheme.
  ! If switch_q_total > 0 then the perturbation is applied.
  if ( vera_visbty%switch_q_total == vera_visbty%switch_q_total_on ) then

    call vera_perturb_q_total( q_total, q_sat, q_total_perturbed )

    ! update the global total q with the pertubed q_total_perturbed
    vera_thread_minpack%q_total_use = q_total_perturbed

  end if

  !=======================================================================
  ! hydrate the aerosol particles
  !=======================================================================

  ! form a first guess as to the aerosol particles' growth factors
  allocate( growth_factor(n_aerosol_species_thread) )
  growth_factor = vera_minpack%initial_growth_guess

  ! hydrate the aerosol particles using the Kohler functions
  call vera_hydrate( growth_factor, vera_kohler_function )

  !=========================================================================
  ! compute the aerosol scattering coefficients
  !=========================================================================

  ! compute the extinction efficiency, Qext
  allocate( mie_qext(n_aerosol_species_thread) )

  call vera_mie_lookup(growth_factor * vera_aerosol_population%rd, mie_qext)

  ! compute the aerosol scattering coefficients
  allocate( beta_aerosol(n_aerosol_species_thread) )

  beta_aerosol = mie_qext * vera_constants%pi *                                &
                 vera_aerosol_population%nc   *                                &
        ( (growth_factor * vera_aerosol_population%rd)**vera_constants%two )

  ! form the sum of the aerosol scattering coefficients
  beta_aerosol_total = sum( beta_aerosol )

  ! check that the aerosol scattering coefficient is not set to NaN,
  ! i.e. Not a Number. If this is set to NaN, then reset it to the
  ! the maximum allowed scattering coefficient, i.e. for minimum
  ! visibility. This test works because NaN is not equal to anything,
  ! not even itself.
  if ( beta_aerosol_total /= beta_aerosol_total ) then
    beta_aerosol_total = vera_koschmeider%max_beta_aerosol
  end if

  ! check that the aerosol scattering coefficient falls within the range
  ! [ vera_koschmeider%min_beta_aerosol, vera_koschmeider%max_beta_aerosol ]
  if ( beta_aerosol_total > vera_koschmeider%max_beta_aerosol ) then
    beta_aerosol_total = vera_koschmeider%max_beta_aerosol
  end if
  if ( beta_aerosol_total < vera_koschmeider%min_beta_aerosol ) then
    beta_aerosol_total = vera_koschmeider%min_beta_aerosol
  end if

  !=========================================================================
  ! compute visibile range from the scattering coefficients using
  ! Koschmeider's law
  !=========================================================================

  ! compute the visibility using just aerosol and Rayleigh scattering, i.e.
  ! excluding scattering from precipitation
  vis_excluding_precip(ii_point) = vera_koschmeider%log_liminal_contrast /     &
                                  ( beta_aerosol_total                   +     &
                                    vera_koschmeider%beta_rayleigh         )

  ! output visibility using just aerosol and Rayleigh scattering, i.e.
  ! no contribution from precipitation
  if (present(vis_no_precip)) then
    vis_no_precip(ii_point) = vis_excluding_precip(ii_point)
  end if

  ! visibility including the contribution from precipitation
  if (present(vis_precip)) then

    ! compute visibility in convective precipitation
    if ( precip_fractional_c(ii_point) > vera_constants%zero ) then
      vis_precip_c(ii_point) = vera_koschmeider%log_liminal_contrast /         &
                               ( beta_aerosol_total                  +         &
                                 vera_koschmeider%beta_rayleigh      +         &
                                 precip_scattering_c(ii_point)         )
    else
      vis_precip_c(ii_point) = vis_excluding_precip(ii_point)
    end if

    ! compute visibility in large scale precipitation
    if ( precip_fractional_ls_no_c(ii_point) > vera_constants%zero) then
      vis_precip_ls(ii_point) = vera_koschmeider%log_liminal_contrast /        &
                               ( beta_aerosol_total                  +         &
                                 vera_koschmeider%beta_rayleigh      +         &
                                 precip_scattering_ls(ii_point)        )
    else
      vis_precip_ls(ii_point) = vis_excluding_precip(ii_point)
    end if

    ! compute the visibility including contributions to scattering
    ! from both large scale and convective precipitation
    vis_precip(ii_point) = min( ( (vera_constants % one                 -      &
                                   precip_fractional_c(ii_point)        -      &
                                   precip_fractional_ls_no_c(ii_point)) *      &
                                  vis_excluding_precip(ii_point) )      +      &
                                ( precip_fractional_ls_no_c(ii_point)   *      &
                                  vis_precip_ls(ii_point) )             +      &
                                ( precip_fractional_c(ii_point)         *      &
                                  vis_precip_c(ii_point) )              ,      &
                                vis_excluding_precip(ii_point) )

  end if

  !=========================================================================
  ! tidy up - deallocate the assumed-shape, semi-dynamic arrays,
  !           do this in the reverse order that they were created,
  !           i.e. get rid of the most recent stuff first
  !=========================================================================

  ! do away with the aerosol scattering coefficients, the extinction
  ! efficiencies and the aerosol growth factors
  if ( allocated( beta_aerosol  ) ) deallocate( beta_aerosol  )
  if ( allocated( mie_qext      ) ) deallocate( mie_qext      )
  if ( allocated( growth_factor ) ) deallocate( growth_factor )

  ! finally, flush out the vera_aerosol_population
  if ( vera_phantom % cast_switch == vera_phantom % cast_switch_on ) then
    if ( allocated( vera_aerosol_population ) )                                &
      deallocate( vera_aerosol_population )
  end if

  ! end of looping over the inputs
end do

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine vera_scheme

end module vera_scheme_mod
