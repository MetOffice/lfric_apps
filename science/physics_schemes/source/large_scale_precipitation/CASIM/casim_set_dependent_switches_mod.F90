! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Cloud Aerosol Interacting Microphysics (CASIM)
! Set variables in module casim_switches dependent on values in
! run_precip namelist

module casim_set_dependent_switches_mod

implicit none

logical ::  l_process             ! True if aerosol processing
logical ::  l_passivenumbers      ! True if aerosol processing and passive
                                  ! numbers in liquid
logical ::  l_passivenumbers_ice  ! True if aerosol processing and passive



character(len=*), parameter, private ::                                        &
  ModuleName='CASIM_SET_DEPENDENT_SWITCHES_MOD'

contains

subroutine casim_set_dependent_switches

use mphys_inputs_mod, only: l_mcr_qrain, l_mcr_qcf2, l_mcr_qgraup,             &
                            l_psd, casim_moments_choice,                       &
                            casim_aerosol_process_level, casim_aerosol_option, &
                            casim_aerosol_couple_choice,                       &
                            l_separate_process_rain, graupel_option,           &
                            no_graupel, casim_iopt_act, fixed_number,          &
                            abdul_razzak_ghan, casim_iopt_inuc

use cloud_inputs_mod,  only: i_cld_vn
use pc2_constants_mod, only: i_cld_off

use casim_switches, only: l_mp_cloudnumber, l_mp_rainnumber, l_mp_rain3mom,    &
                          l_mp_icenumber, l_mp_snownumber, l_mp_snow3mom,      &
                          l_mp_graupnumber, l_mp_graup3mom,                    &
                          l_mp_activesolliquid, l_mp_activesolrain,            &
                          l_mp_activeinsolice, l_mp_activesolice,              &
                          l_mp_activeinsolliquid, l_mp_activesolnumber,        &
                          l_mp_activeinsolnumber, l_casim_warm_only,           &
                          n_casim_progs, casim_moments_option,                 &
                          l_fix_aerosol, l_tracer_aerosol, l_ukca_aerosol,     &
                          l_ukca_feeding_in, l_ukca_feeding_out,               &
                          cloud_mom, rain_mom, ice_mom, snow_mom, graup_mom,   &
                          no_moments, single_moment, double_moment,            &
                          triple_moment, no_processing, passive_processing,    &
                          full_processing, passive_ice_only,                   &
                          passive_liquid_only, no_aerosol_modes,               &
                          soluble_insoluble_modes, soluble_all_modes,          &
                          l_kfsm, l_set_casim_lbc_number

use ukca_option_mod,      only: l_ukca

use model_domain_mod,       only: model_type, mt_lam


use ereport_mod,            only: ereport
use errormessagelength_mod, only: errormessagelength

use yomhook,    only: lhook, dr_hook
use parkind1,   only: jprb, jpim
use umprintmgr, only: newline

implicit none

! ------------------------------------------------------------------------------
! Description:
!   This routine overrides default values held in casim_switches
!   with values depending on the run_precip namelist input. Called from readlsta
!   reconfiguration or scm_shell after run_precip namelist read in.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.1.
!
! ------------------------------------------------------------------------------

! Local parameters for the different settings are contained below

! Parameters for casim_moments_choice
integer, parameter :: all_single_moment = 0
integer, parameter :: all_double_moment = 1
integer, parameter :: warm_cloud2_rain3 = 2
integer, parameter :: warm_cloud1_rain2 = 3
integer, parameter :: warm_cloud1_rain1 = 4
integer, parameter :: warm_cloud1_rain3 = 5
integer, parameter :: warm_cloud2_rain2 = 6
integer, parameter :: double_no_graupel = 7
integer, parameter :: all_moments_on    = 8

! Parameters for casim_aerosol_couple_choice
integer, parameter :: fixed_aerosol      = 0
integer, parameter :: tracer_aerosol     = 1
integer, parameter :: ukca_aerosol_in    = 2
integer, parameter :: ukca_aerosol_inout = 3

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='CASIM_SET_DEPENDENT_SWITCHES'

integer :: errorstatus            ! Return code : 0 Normal Exit : >0 Error
character(len=errormessagelength) :: cmessage
                                  ! Error message if Errorstatus /=0


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!-------------------------------------------------------------------------------
! 1. Initialization of variables
!-------------------------------------------------------------------------------
! Initialise error status to zero (normal/no error)
! Initialise number of extra casim prognostics to 0
errorstatus   = 0
n_casim_progs = 0

! Set all moments to zero initially and then modify depending on
! casim_moments_choice
cloud_mom = 0
rain_mom  = 0
ice_mom   = 0
snow_mom  = 0
graup_mom = 0

! Set l_mcr_qgraup based on graupel option (as required later)
if ( graupel_option > no_graupel ) l_mcr_qgraup = .true.

!-------------------------------------------------------------------------------
! 2.1 Select moments for microphysical species dependent on casim_moments_choice

! Note: casim_moments_option is moments for cloud, rain, ice, snow, graupel
! It is required as an integer in this form for the CASIM code itself
!-------------------------------------------------------------------------------

select case (casim_moments_choice)

case (all_single_moment)
  ! Equivalent to current 3D setup in convective-scale models
  casim_moments_option = 11111
  cloud_mom            = 1
  rain_mom             = 1
  ice_mom              = 1
  snow_mom             = 1
  graup_mom            = 1

  ! Set logical switch for single moment code within the UM.
  ! This lets high-level UM routines know that CASIM is running
  ! as a single moment simulation.
  l_kfsm = .true.

case (all_double_moment)
  ! All species active and double-moment
  casim_moments_option = 22222
  cloud_mom            = 2
  rain_mom             = 2
  ice_mom              = 2
  snow_mom             = 2
  graup_mom            = 2

case (warm_cloud2_rain3)
  ! Warm rain only
  ! Double-moment cloud
  ! Triple-moment rain
  casim_moments_option = 23000
  cloud_mom            = 2
  rain_mom             = 3
  l_casim_warm_only    = .true.

case (warm_cloud1_rain2)
  ! Warm rain only
  ! Single-moment cloud
  ! Double-moment rain
  casim_moments_option = 12000
  cloud_mom            = 1
  rain_mom             = 2
  l_casim_warm_only    = .true.

case (warm_cloud1_rain1)
  ! Warm rain only
  ! Single-moment cloud
  ! Single-moment rain
  casim_moments_option = 11000
  cloud_mom            = 1
  rain_mom             = 1
  l_casim_warm_only    = .true.

case (warm_cloud1_rain3)
  ! Warm rain only
  ! Single-moment cloud
  ! Triple-moment rain
  casim_moments_option = 13000
  cloud_mom            = 1
  rain_mom             = 3
  l_casim_warm_only    = .true.

case (warm_cloud2_rain2)
  ! Warm rain only
  ! Double-moment cloud
  ! Double-moment rain
  casim_moments_option = 22000
  cloud_mom            = 2
  rain_mom             = 2
  l_casim_warm_only    = .true.

case (double_no_graupel)
  ! No graupel
  ! Double-moment cloud
  ! Double-moment rain
  ! Double-moment ice
  ! Double-moment snow
  casim_moments_option = 22220
  cloud_mom            = 2
  rain_mom             = 2
  ice_mom              = 2
  snow_mom             = 2

case (all_moments_on)
  ! All microphysics prognostics on
  ! Useful for prognostic variable testing
  casim_moments_option = 23233
  cloud_mom            = 2
  rain_mom             = 3
  ice_mom              = 2
  snow_mom             = 3
  graup_mom            = 3

case DEFAULT
  ! If this occurs then the user has clearly set up the namelist
  ! incorrectly. Throw an error message and stop the UM immediately
  errorstatus = 1

  write(cmessage, '(A,I0)')                                                    &
    'Variable casim_moments_choice in the run_precip namelist has '// newline//&
    'been set incorrectly. Please check namelist input and correct'// newline//&
    'It should take an integer value >= 0 and agree with one of   '// newline//&
    'the available choices: see documentation for details.'        // newline//&
    'casim_moments_choice =', casim_moments_choice

  call ereport('casim_set_dependent_switches', errorstatus, cmessage)

end select

!-------------------------------------------------------------------------------
! 2.2 Select aerosol modes dependent on activation scheme in use. This saves a
!     switch in the run_precip namelist and prevents the user from setting up
!     a silly combination of switches.
!     However, this option can be expanded or replaced in the future to give
!     more choice to the user, once the CASIM team are confident as to what
!     works sensibly.
!-------------------------------------------------------------------------------

! This is a temporary piece of code to allow UKCA and non-UKCA simulations
! work independently of each other, whilst still minimising the number
! of switches added to the run_precip namelist.
if ( casim_iopt_act >= abdul_razzak_ghan ) then

  ! if iopt_act = 4 then using tracers.
  ! Set up switches:
  ! process_level=3 to enable l_passivenumbers_ice.to keep track of
  !                                         dust numbers in droplets.
  ! Using process_level=2 will make the heterogeneous
  !                         nucleation insensitive to number of dust.
  ! droplet activation set to iopt_act=4: Dust numbers in droplets
  !                                   only available with iopt_act=4.
  ! ice nucleation iopt_inuc=4 as well for a demott scheme using
  !                               dust numbers in activated droplets.
  if ( casim_iopt_act == 4) then
    casim_iopt_inuc=4
    casim_aerosol_process_level=passive_ice_only
    casim_aerosol_couple_choice=1
  end if

  if ( casim_iopt_inuc >= 4 ) then
    ! Need soluble and insoluble modes for more complex ice nucleation schemes
    casim_aerosol_option = soluble_insoluble_modes !( = 3)
  else
    ! Revert to just soluble modes for simpler ice nucleation schemes
    casim_aerosol_option = soluble_all_modes !( = 2)
  end if

  ! Can use this for murk or arcl aerosol. But if
  ! this is a UKCA coupled run too, so set this up.
  if ( l_ukca ) then
    casim_aerosol_couple_choice = ukca_aerosol_in
  end if
else if ( casim_iopt_act == fixed_number ) then
  casim_aerosol_option = no_aerosol_modes        !( = 0)
end if

!-------------------------------------------------------------------------------
! 2.3 Select level of aerosol coupling based on casim_aerosol_couple_choice
!-------------------------------------------------------------------------------

! First, define fixed aerosol if using no_aerosol modes
if ( casim_aerosol_option == no_aerosol_modes ) then
  casim_aerosol_couple_choice = fixed_aerosol
end if

select case (casim_aerosol_couple_choice)

  ! N.B. Most defensive programming and checking for casim aerosol options is
  ! performed in init_casim_run due to ensuring other information required
  ! has been loaded in by other namelists (which may be after this point)

case (fixed_aerosol)
  ! Use a fixed value read from mphys_inputs_mod
  l_fix_aerosol = .true. ! Other related switches already are false
  ! N. B. In the case where casim_aerosol_option == no_modes, the
  ! values will not be read in, but this switch should still be true
  ! otherwise it will try and process aerosol.

case (tracer_aerosol)
  ! Use tracer aerosol data
  l_tracer_aerosol = .true. ! Other related switches already are false

case (ukca_aerosol_in)
  ! Use UKCA aerosol, inwards only
  l_ukca_aerosol    = .true.
  l_ukca_feeding_in = .true. ! UKCA set to feed in only. Other switches
                             ! are already false

case (ukca_aerosol_inout)
  ! Use UKCA aerosol in and out
  l_ukca_aerosol     = .true.
  l_ukca_feeding_in  = .true.
  l_ukca_feeding_out = .true. ! UKCA set to feed in and out. Other switches
                              ! are already false

case DEFAULT
  ! If this occurs then the user has clearly set up the namelist
  ! incorrectly. Throw an error message and stop the UM immediately
  errorstatus = 1
  write(cmessage, '(A,I0)')                                                    &
    'Variable casim_aerosol_couple_choice in the run_precip '      // newline//&
    'namelist has been set incorrectly. Please check namelist'     // newline//&
    'input and correct. It should take an integer value >= 0 and'  // newline//&
    'agree with one of the available choices: see documentation '  // newline//&
    'for details. casim_aerosol_couple_choice =', casim_aerosol_couple_choice

  call ereport('casim_set_dependent_switches', errorstatus, cmessage)

end select

if ( l_ukca_feeding_out .or. ( l_ukca_aerosol .and.                            &
          casim_aerosol_process_level > 0 )) then

    ! UKCA has not yet been coupled to accept processed aerosol from CASIM.
    ! The rose settings should have prevented this from happening, but
    ! if not the model needs to throw an error immediately
    ! Aerosol processing with UKCA on is also currently a bad idea.

  errorstatus = 1
  cmessage    = 'CASIM has been set to feed back  to UKCA    ' //newline//     &
                'aerosol or UKCA is feeding in and processing' //newline//     &
                'is on. These options are not yet available, ' //newline//     &
                'so for your model to run, please choose an  ' //newline//     &
                'alternative aerosol input option in the     ' //newline//     &
                'CASIM settings of the run_precip namelist.'

  call ereport('casim_set_dependent_switches', errorstatus, cmessage)

end if ! l_ukca_feeding_out

!-------------------------------------------------------------------------------
! 2.4 Select moments for microphysical species dependent on casim_aerosol_choice
!-------------------------------------------------------------------------------

! First, define fixed aerosol if using no_aerosol modes
if ( casim_aerosol_option        == no_aerosol_modes .or.                      &
     casim_aerosol_couple_choice == fixed_aerosol         ) then
  casim_aerosol_process_level = no_processing
end if

select case (casim_aerosol_process_level)

case (no_processing)
  l_process            = .false.
  l_passivenumbers     = .false.
  l_passivenumbers_ice = .false.

case (passive_processing)
  l_process            = .true.
  l_passivenumbers     = .true.
  l_passivenumbers_ice = .true.

case (full_processing)
  l_process            = .true.
  l_passivenumbers     = .false.
  l_passivenumbers_ice = .false.

case (passive_ice_only)
  l_process            = .true.
  l_passivenumbers     = .false.
  l_passivenumbers_ice = .true.

case (passive_liquid_only)
  l_process            = .true.
  l_passivenumbers     = .true.
  l_passivenumbers_ice = .false.

case DEFAULT
  ! If this occurs then the user has clearly set up the namelist
  ! incorrectly. Throw an error message and stop the UM immediately
  errorstatus = 1
  write(cmessage, '(A,I0)' )                                                   &
    'Variable casim_aerosol_process_level in the run_precip       '// newline//&
    'namelist has been set incorrectly. Please check namelist     '// newline//&
    'input and correct It should take an integer value >= 0 and   '// newline//&
    'agree with one of the available choices: see documentation   '// newline//&
    'for details. casim_aerosol_process_level =', casim_aerosol_process_level


  call ereport('casim_set_dependent_switches', errorstatus, cmessage)

end select

! Set up prognostic switches based on level of aerosol procesing selected
! and the equivalent on the CASIM repository

if (l_process) then

  if (l_casim_warm_only) then

    l_mp_activesolliquid   = .true.
    l_mp_activeinsolliquid = .false.
    l_mp_activesolice      = .false.
    l_mp_activeinsolice    = .false.

    n_casim_progs = n_casim_progs + 1

  else ! l_casim_warm_only

    l_mp_activesolliquid   = .true.
    l_mp_activeinsolliquid = .true.
    l_mp_activesolice      = .true.
    l_mp_activeinsolice    = .true.

    n_casim_progs = n_casim_progs + 4

  end if ! l_casim_warm_only

  if ( l_passivenumbers ) then
    l_mp_activesolnumber = .true.
    n_casim_progs        = n_casim_progs + 1
  else
    l_mp_activesolnumber = .false.
  end if ! l_passivenumbers

  if ( l_passivenumbers_ice ) then
    l_mp_activeinsolnumber = .true.
    n_casim_progs          = n_casim_progs + 1
  else
    l_mp_activeinsolnumber = .false.
  end if ! l_passivenumbers_ice

  if (l_separate_process_rain) then
    ! Rain is in its own category for processing
    l_mp_activesolrain = .true.
    n_casim_progs      = n_casim_progs + 1
  else
    ! No extra prognostic required for rain
    l_mp_activesolrain = .false.
  end if ! l_separate_process_rain

else ! l_process

  l_mp_activesolliquid    = .false.
  l_mp_activesolrain      = .false.
  l_mp_activeinsolice     = .false.
  l_mp_activesolice       = .false.
  l_mp_activeinsolliquid  = .false.
  l_mp_activesolnumber    = .false.
  l_mp_activeinsolnumber  = .false.

  ! ( Number of CASIM prognostics does not need increasing )

end if ! l_process

!-------------------------------------------------------------------------------
! 3. User input error checking and correction

! While the existing (3D) microphysics and CASIM exist side by side in the UM
! code, options set for the existing microphysics will not work for CASIM.

! The GUI should check inputs and produce a warning, but adding some checking
! code here is a useful 'belt and braces' approach to ensure there are no
! unwanted effects of the code, especially if the user just modifies the
! namelist input, without using the GUI.
!-------------------------------------------------------------------------------

if ( rain_mom == no_moments .and. l_mcr_qrain ) then

  ! Correct conflict between no-rain and prognostic rain

  errorstatus = -100
  cmessage    =                                                     newline//  &
  'Prognostic rain has been switched on (variable l_mcr_qrain   '// newline//  &
  'set to .true. in the run_precip namelist), but this conflicts'// newline//  &
  'with your choice of a no-rain simulation in CASIM. This has  '// newline//  &
  'been corrected. To prevent this warning appearing, set       '// newline//  &
  'l_mcr_qrain to .false.'

  call ereport('casim_set_dependent_switches', errorstatus, cmessage)

  l_mcr_qrain = .false. ! Correct the error

else if ( rain_mom > no_moments                         .and.                  &
          casim_moments_choice > all_single_moment      .and.                  &
          .not. l_mcr_qrain ) then

  ! Correct conflict between rain-included and no prognostic rain
  ! Do not correct if casim_moments_choice set to 0.

  errorstatus = -100
  cmessage    =                                                     newline//  &
  'Prognostic rain has been switched off (variable l_mcr_qrain  '// newline//  &
  'set to .false. in the run_precip namelist), but this         '// newline//  &
  'conflicts with your choice of a simulation in CASIM including'// newline//  &
  'rain. This has been corrected. To prevent this warning       '// newline//  &
  'appearing, set l_mcr_qrain to .true.'

  call ereport('casim_set_dependent_switches', errorstatus, cmessage)


  l_mcr_qrain = .true. ! Correct the error

end if ! rain_mom and l_mcr_qrain

if ( graup_mom == no_moments .and. l_mcr_qgraup ) then

  ! Correct conflict between no graupel in CASIM and prognostic graupel

  errorstatus = -100
  cmessage    =                                                     newline//  &
  'Prognostic graupel has been switched on                      '// newline//  &
  '(graupel_option > 0 in the run_precip namelist),             '// newline//  &
  'but this conflicts with your choice of a no-graupel          '// newline//  &
  'simulation in CASIM. This has been corrected.                '// newline//  &
  'To prevent this warning appearing, set graupel_option = 0    '

  call ereport('casim_set_dependent_switches', errorstatus, cmessage)

  l_mcr_qgraup = .false. ! Correct the error

else if ( graup_mom > no_moments                         .and.                 &
          casim_moments_choice > all_single_moment       .and.                 &
          .not. l_mcr_qgraup ) then

  ! Correct conflict between graupel-included and no prognostic graupel
  ! Do not correct if casim_moments_choice set to 0.

  errorstatus = -100
  cmessage    =                                                       newline//&
  'Prognostic graupel has been switched off (graupel_option set to'// newline//&
  'zero in the run_precip namelist), but this conflicts with      '// newline//&
  'your choice of a graupel-included simulation in CASIM. This    '// newline//&
  'has been corrected. To prevent this warning appearing, set     '// newline//&
  'graupel_option >= 1.'

  call ereport('casim_set_dependent_switches', errorstatus, cmessage)

  l_mcr_qgraup = .true. ! Correct the error

end if ! graup_mom and l_mcr_qgraup

! N.B. CASIM uses the second ice cloud prognostic, which is in mphys_inputs_mod
! but not explicitly set in the run_precip namelist as it has not been used in
! any Wilson and Ballard microphysics jobs in the operational model. Rather than
! include it back in the namelist and have an extra option for users to worry
! about getting correct, it is far simpler to just set it on here, as it will
! always be .true. for CASIM simulations where ice_mom > 0.
if ( ice_mom > no_moments ) l_mcr_qcf2 = .true.

! The second ice cloud prognostic which we have just turned on does not work
! with the generic ice particle size distribution (l_psd set to .true.).
! Check whether the user has either of these turned on and if so,
! throw a warning message and correct.

if ( l_psd ) then

  errorstatus = -100
  cmessage    =                                                     newline//  &
  'The generic ice particle size distribution has been turned on'// newline//  &
  'in your run (l_psd set to .true. n the run_precip namelist), '// newline//  &
  'but this does not work with CASIM. This has been corrected.  '// newline//  &
  'To remove this warning please set this logical to .false.'

  call ereport('casim_set_dependent_switches', errorstatus, cmessage)

  ! Now correct the error
  l_psd        = .false.

end if ! l_psd

!------------------------------------------------------------------------------
! 4.1 Determine (liquid) cloud moment and initialise switches
!------------------------------------------------------------------------------

if (cloud_mom == double_moment) then
  l_mp_cloudnumber = .true.
  n_casim_progs    = n_casim_progs + 1
end if

!-------------------------------------------------------------------------------
! 4.2 Determine ice cloud moment and initialise switches
!------------------------------------------------------------------------------

if (ice_mom == double_moment) then
  l_mp_icenumber = .true.
  n_casim_progs  = n_casim_progs + 1
end if

!------------------------------------------------------------------------------
! 4.3 Determine rain moment and initialise switches
!------------------------------------------------------------------------------
select case (rain_mom)

case (single_moment) ! Single-moment rain

  l_mp_rainnumber = .false.
  l_mp_rain3mom   = .false.

case (double_moment) ! Double-moment rain

  l_mp_rainnumber = .true.
  l_mp_rain3mom   = .false.
  n_casim_progs   = n_casim_progs + 1

case (triple_moment) ! Triple-moment rain

  l_mp_rainnumber = .true.
  l_mp_rain3mom   = .true.
  n_casim_progs   = n_casim_progs + 2

end select

!------------------------------------------------------------------------------
! 4.4 Determine snow moment and initialise switches
!------------------------------------------------------------------------------
select case (snow_mom)

case (no_moments) ! No snow

  l_mp_snownumber = .false.
  l_mp_snow3mom   = .false.

case (single_moment) ! Single-moment snow

  l_mp_snownumber = .false.
  l_mp_snow3mom   = .false.

case (double_moment) ! Double-moment snow

  l_mp_snownumber = .true.
  l_mp_snow3mom   = .false.
  n_casim_progs   = n_casim_progs + 1

case (triple_moment) ! Triple-moment snow

  l_mp_snownumber = .true.
  l_mp_snow3mom   = .true.
  n_casim_progs   = n_casim_progs + 2

end select

!------------------------------------------------------------------------------
! 4.5 Determine graupel moment and initialise switches
!------------------------------------------------------------------------------
select case (graup_mom)

case (no_moments) ! No graupel

  l_mp_graupnumber = .false.
  l_mp_graup3mom   = .false.

case (single_moment) ! Single-moment graupel

  l_mp_graupnumber = .false.
  l_mp_graup3mom   = .false.

case (double_moment) ! Double-moment graupel

  l_mp_graupnumber = .true.
  l_mp_graup3mom   = .false.
  n_casim_progs    = n_casim_progs + 1

case (triple_moment) ! Triple-moment graupel

  l_mp_graupnumber = .true.
  l_mp_graup3mom   = .true.
  n_casim_progs    = n_casim_progs + 2

end select

!set l_set_casim_lbc_number based on global or lam.
!If rim size was set to 0 for
!global this would not be necessary.
!Or if global model passes number conc to LAM
if ( ( model_type == mt_lam ) .and. ( i_cld_vn /= i_cld_off ) ) then
  l_set_casim_lbc_number = .true.
else
  l_set_casim_lbc_number = .false.
end if





!------------------------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine casim_set_dependent_switches

!==============================================================================
! Break between subroutines
!==============================================================================

subroutine casim_print_dependent_switches

use casim_switches, only: l_mp_cloudnumber, l_mp_rainnumber, l_mp_rain3mom,    &
                          l_mp_icenumber, l_mp_snownumber, l_mp_snow3mom,      &
                          l_mp_graupnumber, l_mp_graup3mom,                    &
                          l_mp_activesolliquid, l_mp_activesolrain,            &
                          l_mp_activeinsolice, l_mp_activesolice,              &
                          l_mp_activeinsolliquid, l_mp_activesolnumber,        &
                          l_mp_activeinsolnumber, n_casim_progs,               &
                          casim_moments_option, l_fix_aerosol,                 &
                          l_tracer_aerosol, l_ukca_aerosol, l_ukca_feeding_in, &
                          l_ukca_feeding_out

use mphys_inputs_mod, only: casim_aerosol_option

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

use umprintmgr, only: umprint, ummessage, PrNorm

implicit none

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='CASIM_PRINT_DEPENDENT_SWITCHES'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

call umprint('Model run includes CASIM Microphysics',                          &
              src='casim_set_dependent_switches_mod')

call umprint('Switches set for CASIM based on inputs from namelist run_precip',&
              src='casim_set_dependent_switches_mod')

write(ummessage, '(A, L1)')'l_mp_cloudnumber = ', l_mp_cloudnumber
call umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

write(ummessage, '(A, L1)')'l_mp_rainnumber = ', l_mp_rainnumber
call umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

write(ummessage, '(A, L1)')'l_mp_rain3mom = ', l_mp_rain3mom
call umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

write(ummessage, '(A, L1)')'l_mp_icenumber = ', l_mp_icenumber
call umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

write(ummessage, '(A, L1)')'l_mp_snownumber = ', l_mp_snownumber
call umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

write(ummessage, '(A, L1)')'l_mp_snow3mom = ', l_mp_snow3mom
call umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

write(ummessage, '(A, L1)')'l_mp_graupnumber = ', l_mp_graupnumber
call umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

write(ummessage, '(A, L1)')'l_mp_graup3mom = ', l_mp_graup3mom
call umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

call umprint('- - - - - - - - - - - - - - - - - - - - - - - - - - ',           &
             src='casim_set_dependent_switches_mod')

write(ummessage, '(A, L1)')'l_mp_activesolliquid = ', l_mp_activesolliquid
call umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

write(ummessage, '(A, L1)')'l_mp_activesolrain = ', l_mp_activesolrain
call umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

write(ummessage, '(A, L1)')'l_mp_activeinsolice = ', l_mp_activeinsolice
call umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

write(ummessage, '(A, L1)')'l_mp_activesolice = ', l_mp_activesolice
call umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

write(ummessage, '(A, L1)')'l_mp_activeinsolliquid = ', l_mp_activeinsolliquid
call umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

write(ummessage, '(A, L1)')'l_mp_activesolnumber = ', l_mp_activesolnumber
call umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

write(ummessage, '(A, L1)')'l_mp_activeinsolnumber = ', l_mp_activeinsolnumber
call umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

write(ummessage,'(A,I5)')'casim_moments_option = ', casim_moments_option
call umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

call umprint('- - - - - - - - - - - - - - - - - - - - - - - - - - ',           &
             src='casim_set_dependent_switches_mod')

write(ummessage,'(A,I2)')'Number of additional CASIM prognostics = ',          &
                          n_casim_progs
call umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

call umprint('- - - - - - - - - - - - - - - - - - - - - - - - - - ',           &
             src='casim_set_dependent_switches_mod')

write(ummessage,'(A,I2)')'casim_aerosol_option = ',                            &
                          casim_aerosol_option
call umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

write(ummessage, '(A, L1)')'l_fix_aerosol = ', l_fix_aerosol
call umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

write(ummessage, '(A, L1)')'l_tracer_aerosol = ', l_tracer_aerosol
call umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

write(ummessage, '(A, L1)')'l_ukca_aerosol = ', l_ukca_aerosol
call umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

write(ummessage, '(A, L1)')'l_ukca_feeding_in = ', l_ukca_feeding_in
call umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

write(ummessage, '(A, L1)')'l_ukca_feeding_out = ', l_ukca_feeding_out
call umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

call umprint('- - - - - - end of CASIM switches - - - - - -',                  &
    src='casim_set_dependent_switches_mod')

!------------------------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine casim_print_dependent_switches

end module casim_set_dependent_switches_mod
