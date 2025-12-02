! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!  data module for switches/options concerned with the precipitation scheme.
  ! Description:
  !   Module containing runtime options/data used by the precipitation scheme

  ! Method:
  !   Switches and associated data values used by the precipitation scheme
  !   are defined here and assiged default values. These may be overridden
  !   by namelist input.

  !   A description of what each switch or number refers to is provided
  !   with the namelist

  !   Any routine wishing to use these options may do so with the 'use'
  !   statement.
  !
  ! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

  ! Code Description:
  !   Language: FORTRAN 90


module mphys_inputs_mod

use missing_data_mod, only: rmdi, imdi
use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim
use errormessagelength_mod, only: errormessagelength
use ereport_mod, only: ereport
use umprintmgr,  only: newline

use um_types, only: real_umphys

implicit none

!===========================================================================
! integer options set from RUN_PRECIP namelist
!===========================================================================

! Method for iterating microphysics scheme
integer :: i_mcr_iter       = imdi
! Values to which this can be set
! No iterating the microphysics
integer, parameter :: i_mcr_iter_none   = 0
! Iterating on: user sets number of iterations per model timestep
integer, parameter :: i_mcr_iter_niters = 1
! Iterating on: user set microphysics "timestep" and number of
!               iterations is calculated
integer, parameter :: i_mcr_iter_tstep  = 2
!------------------------------------

! Number of iterations of microphysics scheme
!  Input parameter for i_mcr_iter=i_mcr_iter_niters
!  Calculated at run-time for i_mcr_iter=i_mcr_iter_tstep
integer :: niters_mp        = imdi
!------------------------------------

! Requested (input) length of microphysics timestep (s)
!  Input parameter for i_mcr_iter=i_mcr_iter_tstep
integer :: timestep_mp_in    = imdi

!------------------------------------

! Choice for each microphysical species for CASIM
integer :: casim_moments_choice = imdi

!------------------------------------
! Option for the graupel scheme
integer :: graupel_option = imdi

! Integer options for the new graupel scheme (graupel_option)
integer, parameter :: no_graupel   = 0  ! No prognostic graupel
integer, parameter :: gr_orig      = 1  ! Original scheme
integer, parameter :: gr_srcols    = 2  ! Snow-rain collections produce graupel
integer, parameter :: gr_field_psd = 3  ! New (Field et al) scheme

!-------------------------------------
! Option for the location of sedimentation within the microphysics
integer :: sediment_loc = imdi

! Integer options for the sedimentation
integer, parameter :: all_sed_start = 1 ! Historical choice, sedimentation
                                        ! at start of microphysics
integer, parameter :: fall_end      = 2 ! Move lsp_fall call to the end
integer, parameter :: all_sed_end   = 3 ! Move lsp_fall and lsp_settle to the
                                        ! end of the microphysics
integer, parameter :: rain_sed_end  = 4 ! Just move the rain sedimentation
                                        ! to the end of the microphysics
integer, parameter :: warm_sed_end  = 5 ! Have all warm sedimentation
                                        ! (lsp_settle and lsp_fall_rain)
                                        ! at the end of the microphysics

! Activation option for CASIM
!(0- fixed number, 3 -ukca/murk/arcl aerosol, 4- tracers)
integer :: casim_iopt_act = imdi

! Parameter options for casim_iopt_act
integer, parameter :: fixed_number      = 0
integer, parameter :: abdul_razzak_ghan = 3
integer, parameter :: arg_dust_ft       = 4

! Options for how to update the prognostic precip fraction consistent
! with the increments to precip mass:
integer :: i_update_precfrac = imdi
integer, parameter :: i_homog_areas = 1  ! Partition the precip fraction and
                                         ! its increments into homogeous areas
integer, parameter :: i_sg_correl = 2    ! Assume sub-grid spatial correlation
                                         ! between precip frac and increments.

!===========================================================================
! logical options set from RUN_PRECIP namelist
!===========================================================================

! Use improved warm rain microphysics scheme
logical :: l_warm_new       = .false.
!--------------------------------------

! Use generic ice psd (particle size distribution)
! ice psd
logical :: l_psd            = .false.

! Use global version - hard-wired in mphys_constants_mod
! l_psd_global     = .true.
!------------------------------------

! Use murk aerosol to
! calculate the droplet number
logical :: l_autoconv_murk  = .false.
!------------------------------------

! Enable tapering of cloud droplets
! towards surface
logical :: l_droplet_tpr    = .false.
!------------------------------------

! Use Aerosol climatologies to generate drop number
logical :: l_mcr_arcl       = .false.
!-----------------------------------

! Include prognostic rain
logical :: l_mcr_qrain      = .false.
!-----------------------------------

! Prognostic rain lbcs active
logical :: l_mcr_qrain_lbc  = .false.
!-----------------------------------

! Prognostic graupel lbcs active
logical :: l_mcr_qgraup_lbc = .false.
!-----------------------------------

! Turns precipitation code on/off
logical :: l_rain           = .false.
!-----------------------------------

! Turns seeder-feeder code on/off
logical :: l_orograin       = .false.
!-----------------------------------

! In seeder feeder code enhance riming also
logical :: l_orogrime       = .false.
!-----------------------------------

! Account for blocking in seeder feeder
logical :: l_orograin_block = .false.
!-----------------------------------

! Use same lwc FSD for autoconv and accretion as in cloud generator
logical :: l_fsd_generator = .false.

! Use different fallspeed relations for
! crystals and aggregates with the generic psd
logical :: l_diff_icevt = .false.
!-----------------------------------

! Use sulphate aerosol in microphysics
logical :: l_use_sulphate_autoconv = .false.

! Use sea-salt aerosol in microphysics
logical :: l_use_seasalt_autoconv = .false.

! Use biomass aerosol in microphysics
logical :: l_use_bmass_autoconv = .false.

! Use fossil-fuel organic carbon in microphysics
logical :: l_use_ocff_autoconv = .false.

! Use ammonium nitrate aerosol in microphysics
logical :: l_use_nitrate_autoconv = .false.

! Use autoconversion de-biasing scheme in microphysics
logical :: l_auto_debias = .false.

! Produce extra qcl by turbulent processes
logical :: l_subgrid_qcl_mp = .false.

! Apply capacitance of ice crystals changes
logical :: l_ice_shape_parameter = .false.

! Logical for using tnuc_dust
logical :: l_progn_tnuc = .false.

! Apply the shape-dependent riming parametrization
logical :: l_shape_rime = .false.

! CASIM Microphysics- if switch is True, run with CASIM, otherwise
! use Wilson and Ballard Microphysics
logical :: l_casim = .false.

! Use heterogeneous nucleation of rain
logical :: l_het_freezing_rain = .false.

! Include R^2 factor in layer-masses rho dz, to improve conservation of
! moisture in model-runs that don't use the shallow atmosphere approximation.
logical :: l_mphys_nonshallow = .false.

! Use prognostic precipitation fraction for rain (and optionally graupel)
logical :: l_mcr_precfrac = .false.

! Apply the prognostic precipitation fraction for graupel as well as rain
logical :: l_subgrid_graupel_frac = .false.

! Apply microphysics increments to the precip fluxes as well as prognostics
! on level k (fixes a numerical problem where most of the rain or graupel
! falls straight through a model-level in 1 sub-step, so that there's not
! enough left on level k to evaporate; this causes overestimation of the
! amount of precip reaching the surface, vs evaporating on the way down).
logical :: l_proc_fluxes = .false.


!#############################
!to allow revert to ral3.1

! needed for DA cycling (set to T) so that number concs are not reset.
logical :: l_casim_warmstart = .false.  ! ral3.1 default
! needed for microphysics (set to T do micro in rim - preferred option)
logical :: l_micro_in_rim = .false.     ! ral3.1 default

!===========================================================================
! real values set from RUN_PRECIP namelist
!===========================================================================

! Rain particle size distribution
! values
real(kind=real_umphys) :: x1r                 = rmdi
real(kind=real_umphys) :: x2r                 = rmdi
!------------------------------------

! Ice mass-diameter relationship
! values
real(kind=real_umphys) :: ai                  = rmdi
real(kind=real_umphys) :: bi                  = rmdi
!------------------------------------

!-----------------------------------------------------------------------------
! Fallspeed parameters for crystals and aggregates to use with generic psd to
! us when l_diff_fallspeed = .true.
!-----------------------------------------------------------------------------
!-- Linear ice vt ------------------------------------------------------------
! Fallspeed parameters for crystals
real(kind=real_umphys) :: cic_input = rmdi !1042.18
real(kind=real_umphys) :: dic_input = rmdi !1.0
! Fallspeed parameters for aggregates
real(kind=real_umphys) :: ci_input = rmdi !14.2611
real(kind=real_umphys) :: di_input = rmdi !0.416351
!-----------------------------------------------------------------------------

! Droplet number at z_surf and below:

real(kind=real_umphys) :: ndrop_surf          = rmdi

! Height at which droplet number reaches ndrop_surf
real(kind=real_umphys) :: z_surf = rmdi

!------------------------------------

! Cloud-rain correlation coefficient for inhomogeneity parametrization:
! i.e. when l_inhomog=.true.

real(kind=real_umphys) :: c_r_correl = rmdi
!------------------------------------

! Cloud water content exponent used in autoconversion rate calculation
! if using Khairoutdinov & Kogan (2000, MWR) parameterization
! (l_warm_new=.true.)
real(kind=real_umphys) :: aut_qc = rmdi
!------------------------------------

! Aerosol climatology scaling factor, to account for spatial/temporal
! inhomogeneity
real(kind=real_umphys) :: arcl_inhom_sc  = rmdi
!------------------------------------

! Axial ratios for aggregates and crystals
real(kind=real_umphys) :: ar  = rmdi
real(kind=real_umphys) :: arc = rmdi
!-----------------------------------

! Scaling value for grid spacing
real(kind=real_umphys) :: mp_dz_scal   = rmdi

! Shape dependent riming rate parameters
!
real(kind=real_umphys) :: qclrime = rmdi
real(kind=real_umphys) :: a_ratio_fac = rmdi
real(kind=real_umphys) :: a_ratio_exp = rmdi

! Critical Froude number for blocking in Seeder Feeder scheme
real(kind=real_umphys) :: fcrit = rmdi
! Scaling parameter for sub-grid hill peak-to-trough height
real(kind=real_umphys) :: nsigmasf = rmdi
! Scaling parameter for horizontal wavelength of sub-grid sinusoidal
! ridge for use in vertical decay function
real(kind=real_umphys) :: nscalesf = rmdi

! factor for subgrid vertical velocity variance to be added to
! explicit velocity to be used in aerosol activation
real(kind=real_umphys) :: wvarfac = rmdi

! Parameter used in mixed phase scheme
! in the second-order Lagrangian structure function
real(kind=real_umphys) :: mp_czero = rmdi

! Parameter for limit on the phase relaxation timescale
! Used to avoid complications with very non-turbulent situations
! and large phase relaxation timescales
real(kind=real_umphys) :: mp_tau_lim = rmdi

! Factor for enhancing rain evaporation at high rain mixing-ratios,
! to account for droplet break-up
real(kind=real_umphys) :: heavy_rain_evap_fac = rmdi

!----------------------------------------------------------------------

! Define the RUN_PRECIP namelist

namelist/run_precip/                                                           &
       l_warm_new, l_psd, l_autoconv_murk,                                     &
       l_use_sulphate_autoconv, l_use_seasalt_autoconv, l_use_bmass_autoconv,  &
       l_use_ocff_autoconv, l_use_nitrate_autoconv,                            &
       x1r, x2r, c_r_correl, aut_qc, ai, bi, ar, arc,                          &
       cic_input,dic_input,ci_input,di_input,l_diff_icevt,                     &
       ndrop_surf, z_surf, l_droplet_tpr,                                      &
       i_mcr_iter, niters_mp, timestep_mp_in,                                  &
       l_mcr_arcl, arcl_inhom_sc,                                              &
       l_mcr_qrain, graupel_option, l_mcr_qrain_lbc, l_mcr_qgraup_lbc,         &
       l_rain, l_fsd_generator, l_subgrid_qcl_mp, l_ice_shape_parameter,       &
       l_progn_tnuc, mp_dz_scal,l_shape_rime, qclrime, a_ratio_fac,            &
       a_ratio_exp, l_orograin, l_orogrime, l_orograin_block,                  &
       l_het_freezing_rain, sediment_loc,                                      &
       fcrit, nsigmasf, nscalesf, wvarfac, heavy_rain_evap_fac,                &
       l_casim, casim_moments_choice,                                          &
       casim_iopt_act, l_mphys_nonshallow,                                     &
       l_mcr_precfrac, l_subgrid_graupel_frac, l_proc_fluxes,                  &
       i_update_precfrac,                                                      &
       mp_czero, mp_tau_lim, l_casim_warmstart, l_micro_in_rim

!===========================================================================
! logical options not set in namelist
!===========================================================================
! Do not use generic ice particle size distribution in calculations
logical, parameter :: not_generic_size_dist = .false.

! Second ice variable lbcs active
logical :: l_mcr_qcf2_lbc   = .false.
!-----------------------------------

! Include second ice variable
logical :: l_mcr_qcf2       = .false.
!-----------------------------------

! Use soot aerosol in microphysics (not set in namelist)
logical :: l_use_soot_autoconv = .false.
!-----------------------------------

! Include prognosic graupel
logical :: l_mcr_qgraup     = .false.
!-----------------------------------

!==========================================================================
! CASIM Options for future inclusion in namelist
!==========================================================================

! These variables are intended to form part of the run_precip namelist
! in the near future and hence are set and initialised here to values which
! are known to work with CASIM at present.

! Choice for the aerosol modes option for CASIM
integer :: casim_aerosol_option = 0

! Choice for the aerosol coupling method in CASIM
integer :: casim_aerosol_couple_choice = 0

! Choice for the aerosol processing level in CASIM
integer :: casim_aerosol_process_level = 0

! Fixed accumulation model mass and number
real(kind=real_umphys), parameter ::                                           &
  fixed_accum_sol_mass   = 1.5e-9_real_umphys
real(kind=real_umphys), parameter ::                                           &
  fixed_accum_sol_number = 100.0e6_real_umphys

! CASIM sedimentation timestep in seconds
real(kind=real_umphys), parameter :: casim_max_sed_length = 120.0_real_umphys

! Ice Nucleation option for CASIM - currently set as default (1; Cooper ice
! nucleation). N.B. Values set at 4 or above need insoluble modes in CASIM.
integer :: casim_iopt_inuc = 1

! CASIM aerosol processing has separate rain category
logical, parameter :: l_separate_process_rain = .false.

! DrHook-related parameters
integer(kind=jpim), parameter, private :: zhook_in  = 0
integer(kind=jpim), parameter, private :: zhook_out = 1

logical :: l_no_cf = .false.
! A temporary logical for an advanced user to turn off the effect of cloud
! fraction in the radar reflectivity code. This is not intended to be added
! to the namelist (where it may confuse the user) and is just for testing
! and checking of the radar reflectivity code.

character(len=*), parameter, private :: ModuleName='MPHYS_INPUTS_MOD'

contains

subroutine print_nlist_run_precip()
use umPrintMgr, only: umPrint
implicit none
character(len=50000) :: lineBuffer
real(kind=jprb) :: zhook_handle

character(len=*), parameter :: RoutineName='PRINT_NLIST_RUN_PRECIP'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

call umPrint('Contents of namelist run_precip',                                &
    src='mphys_inputs_mod')

write(lineBuffer,'(A,L1)')' l_warm_new = ',l_warm_new
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,L1)')' l_psd = ',l_psd
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,L1)')' l_autoconv_murk = ',l_autoconv_murk
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,L1)') ' l_use_sulphate_autoconv = ',l_use_sulphate_autoconv
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,L1)') ' l_use_seasalt_autoconv = ',l_use_seasalt_autoconv
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,L1)') ' l_use_ocff_autoconv = ',l_use_ocff_autoconv
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,L1)') ' l_use_nitrate_autoconv = ',l_use_nitrate_autoconv
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,L1)') ' l_use_bmass_autoconv = ',l_use_bmass_autoconv
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,F0.4)')' x1r = ',x1r
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,F0.4)')' x2r = ',x2r
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,F0.4)')' c_r_correl = ',c_r_correl
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,F0.4)')' aut_qc = ',aut_qc
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,F0.4)')' ai = ',ai
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,F0.4)')' bi = ',bi
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,F0.4)')' ar = ',ar
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,F0.4)')' arc = ',arc
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,F0.4)')' ndrop_surf = ',ndrop_surf
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,F0.4)')' z_surf = ',z_surf
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,L1)')' l_droplet_tpr = ',l_droplet_tpr
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,I0)')' i_mcr_iter = ',i_mcr_iter
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,I0)')' niters_mp = ',niters_mp
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,I0)')' timestep_mp_in = ',timestep_mp_in
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,L1)')' l_mcr_arcl = ',l_mcr_arcl
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,F0.4)')' arcl_inhom_sc = ',arcl_inhom_sc
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,L1)')' l_mcr_qrain =', l_mcr_qrain
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,L1)')' l_mcr_qrain_lbc =', l_mcr_qrain_lbc
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,L1)')' l_mcr_qgraup_lbc =', l_mcr_qgraup_lbc
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,L1)')' l_rain =',l_rain
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,L1)')' l_fsd_generator =',l_fsd_generator
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,L1)')' l_subgrid_qcl_mp =',l_subgrid_qcl_mp
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,L1)')' l_ice_shape_parameter =',l_ice_shape_parameter
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,L1)')' l_progn_tnuc =',l_progn_tnuc
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,F0.4)')' mp_dz_scal =', mp_dz_scal
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,L1)')' l_shape_rime =', l_shape_rime
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,F0.4)')' qclrime =', qclrime
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,F0.4)')' a_ratio_fac =', a_ratio_fac
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,F0.4)')' a_ratio_exp =', a_ratio_exp
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,L1)')' l_orograin =', l_orograin
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,L1)')' l_orogrime =', l_orogrime
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,L1)')' l_orograin_block =', l_orograin_block
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,F0.6)')' fcrit =', fcrit
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,F0.6)')' nsigmasf =', nsigmasf
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,F0.6)')' nscalesf =', nscalesf
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,F0.6)')' wvarfac =', wvarfac
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,F0.6)')' heavy_rain_evap_fac =', heavy_rain_evap_fac
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,L1)')' l_casim =', l_casim
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,I0)')' casim_moments_choice =', casim_moments_choice
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,L1)')' l_het_freezing_rain =', l_het_freezing_rain
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,L1)')' l_mphys_nonshallow =', l_mphys_nonshallow
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,I0)')' sediment_loc = ',sediment_loc
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,fmt='(A,I0)')'casim_iopt_act =', casim_iopt_act
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,L1)')' l_mcr_precfrac = ',l_mcr_precfrac
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,L1)')' l_subgrid_graupel_frac = ',l_subgrid_graupel_frac
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,L1)')' l_proc_fluxes = ',l_proc_fluxes
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,fmt='(A,I0)')'i_update_precfrac =', i_update_precfrac
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,F0.6)')' mp_czero =', mp_czero
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,F0.6)')' mp_tau_lim =', mp_tau_lim
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,L1)')' l_casim_warmstart = ',l_casim_warmstart
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,L1)')' l_micro_in_rim = ',l_micro_in_rim
call umPrint(lineBuffer,src='mphys_inputs_mod')
call umPrint('- - - - - - end of namelist - - - - - -',                        &
    src='mphys_inputs_mod')

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine print_nlist_run_precip


subroutine check_run_precip

use chk_opts_mod,        only: chk_var, def_src
use murk_inputs_mod,     only: l_murk
use cv_run_mod,          only: i_convection_vn, i_cv_comorph

implicit none

character(len=*), parameter :: RoutineName='CHECK_RUN_PRECIP'

character(len=errormessagelength) :: comments
character(len=100) :: ChkStr

integer :: ErrorStatus

real(kind=jprb) :: zhook_handle

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
def_src = RoutineName

ErrorStatus = 0

! Set graupel option to 'OFF' in the case that the whole
! precipitation scheme is off
if (.not. l_rain) then
  graupel_option = no_graupel
end if

! Option to set l_mcr_qgraup depending on graupel_option
! This code is temporary, pending the removal of l_mcr_qgraup
! from the UM code base.
if ( graupel_option == no_graupel ) then
  l_mcr_qgraup = .false.
else
  l_mcr_qgraup = .true.
end if
! End of temporary code based on graupel_option

! Aerosol indirect effect
if ( l_autoconv_murk .and. l_mcr_arcl ) then
  comments='Cannot set both l_autoconv_murk and l_mcr_arcl to .true.'
  ErrorStatus=1
  call ereport(RoutineName,ErrorStatus,comments)
end if

! Comorph convection scheme requires large-scale precip and
! prognostic precip fraction to be on
if ( i_convection_vn == i_cv_comorph ) then
  if ( .not. ( l_rain .and. l_mcr_precfrac ) ) then
    comments='Running with the comorph convection scheme '         //newline// &
             'requires using the large-scale precipitation '       //newline// &
             'scheme (l_rain=.true.) with the prognostic '         //newline// &
             'precipitation fraction turned on '                   //newline// &
             '(l_mcr_precfrac=.true.).'
    ErrorStatus=100
    call ereport(RoutineName,ErrorStatus,comments)
  end if
end if

! These parameters will always need checking as they are used if the bimodal
! scheme is used
call chk_var(mp_czero,  'mp_czero', '[2.0:12.0]')
call chk_var(mp_tau_lim,  'mp_tau_lim', '[600.0:9600.0]')

if (l_rain) then
  ! Only bother checking the namelist if the precip scheme is in use.

  if (l_casim) then
    call chk_var( casim_moments_choice, 'casim_moments_choice', '[0, 1, 8]')

    if (l_psd) then

      ! This does not work with CASIM and can cause user
      ! jobs to hang without crashing if used together
      ! Best is to Call Ereport to prevent this happening.
      ErrorStatus = 100

      comments = 'Generic ice PSD does not work with CASIM.     '// newline//  &
      'Please open the Rose Gui and switch off l_psd.           '// newline//  &
      'This model run crash is to save you the annoyance of it  '// newline//  &
      'hanging and later timing out without any obvious error'

      call ereport(RoutineName, ErrorStatus, comments)

    end if

    write(ChkStr,'(3(A,I1),A)') '[',fixed_number,',',                          &
                                          abdul_razzak_ghan,',',arg_dust_ft,']'
    call chk_var(casim_iopt_act, 'casim_iopt_act', ChkStr)

    call chk_var(wvarfac, 'wvarfac', '[0.0:5.0]')

  else ! l_casim

    ! Check options for the Wilson and Ballard Microphysics

    write(ChkStr, '(3(A,I1),A)')                                               &
         '[',i_mcr_iter_none,',',i_mcr_iter_niters,',',i_mcr_iter_tstep,']'

    ! First check i_mcr_iter and fall over if it is not set correctly
    call chk_var(i_mcr_iter, 'i_mcr_iter', ChkStr)

    select case (i_mcr_iter)

      ! N.B. Case of i_mcr_iter_none has nothing to check

    case (i_mcr_iter_niters)
      ! Check number of iterations matches metadata
      call chk_var(niters_mp, 'niters_mp', '[1:999]')

    case (i_mcr_iter_tstep)
      ! Check microphysics timestep matches metadata
      call chk_var(timestep_mp_in, 'timestep_mp_in', '[1:3600]')

    end select

    ! Checking of graupel_option
    write(ChkStr, '(4(A,I1),A)')                                               &
         '[',no_graupel,',',gr_orig,',',gr_srcols,',',gr_field_psd,']'

    call chk_var(graupel_option, 'graupel_option', ChkStr)

    ! Checking of sediment_loc
    write(ChkStr, '(5(A,I1),A)')                                               &
         '[',all_sed_start,',',fall_end,',',all_sed_end,',',rain_sed_end,      &
         ',',warm_sed_end,']'

    call chk_var(sediment_loc, 'sediment_loc', ChkStr)

    ! The options below will always need checking in the run
    call chk_var(x1r, 'x1r', '[>=0.0]')
    call chk_var(x2r, 'x2r', '[-10.0:10.0]')
    call chk_var(ai,  'ai', '[>=0.0]')
    call chk_var(bi,  'bi', '[>=0.0]')
    call chk_var(heavy_rain_evap_fac, 'heavy_rain_evap_fac', '[0.0:1.0E9]')

    if (.not. l_ice_shape_parameter) then
        !ar value is used from namelist only with l_ice_shape_parameter = false
      call chk_var(ar,  'ar', '[0.01:10.0]')
    end if

    if (l_psd .and. l_diff_icevt) then

      ! With both options set, need to check the input parameters to the
      ! split fall speed (cic_input, dic_input, ci_input and di_input)

      call chk_var(cic_input, 'cic_input', '[0.0:1.0E7]')
      call chk_var(dic_input, 'dic_input', '[0.0:1.0E7]')
      call chk_var(ci_input,  'ci_input', '[0.0:1.0E7]')
      call chk_var(di_input,  'di_input', '[0.0:1.0E7]')

    else if (.not. l_psd) then
      ! Run is not using the Generic Ice PSD. Must check all crystal options
      ! are correctly set.

      call chk_var(arc, 'arc', '[0.01:10.0]')

    end if ! not l_psd

    if (l_warm_new) call chk_var(c_r_correl, 'c_r_correl', '[-1.0:1.0]')

    if (l_warm_new) call chk_var(aut_qc, 'aut_qc', '[2.0:3.5]')

    if (l_mcr_arcl) call chk_var(arcl_inhom_sc, 'arcl_inhom_sc', '[0.01:10.0]')

    if (l_subgrid_qcl_mp) then

      call chk_var(mp_dz_scal, 'mp_dz_scal', '[0.01:10.0]')

    end if ! l_subgrid_qcl_mp

    if (l_shape_rime) then

      call chk_var(a_ratio_fac, 'a_ratio_fac','[0.0:1.0]')
      call chk_var(a_ratio_exp, 'a_ratio_exp','[-1.0:0.0]')
      call chk_var(qclrime,     'qclrime','[0.0:1.0E-2]')

    end if

    if (l_mcr_qcf2 .and. l_psd) then

      ! These do not work together and have caused jobs
      ! to hang in the past when someone set these logicals
      ! to true. Thus it is probably best to throw an error
      ! right now.

      ErrorStatus = 100

      comments = 'Generic ice PSD does not work with second ice '// newline//  &
      'crystal prognostic. Please open the Rose Gui and switch  '// newline//  &
      'off l_psd.'

      call ereport(RoutineName, ErrorStatus, comments)

    end if ! l_mcr_qcf2 / l_psd

    if (l_autoconv_murk .and. .not. l_murk) then

      ! Someone has set up MURK aerosol for autoconversion without the
      ! Murk aerosol on. This will not work.
      ! Note also that the l_murk variable must be set before this checking
      ! takes place.

      ErrorStatus = 100

      comments = 'l_autoconv_murk has been set without l_murk   '// newline//  &
      'i.e. requesting Murk aerosol to do some autoconversion   '// newline//  &
      'from cloud to rain when it is not available. Please open '// newline//  &
      'the GUI and switch either l_autoconv_murk off or l_murk on.'

      call ereport(RoutineName, ErrorStatus, comments)

    end if
    if (l_orograin) then

      if (l_orograin_block) then
        call chk_var(fcrit, 'fcrit','[>0.1]')
      end if

      call chk_var(nsigmasf, 'nsigmasf', '[>=0.0]')
      call chk_var(nscalesf, 'nscalesf', '[0.01:6.0]')

    end if ! l_orograin

    ! Check prognostic precipitation fraction options...
    if ( l_mcr_precfrac ) then

      call chk_var(i_update_precfrac, 'i_update_precfrac',                     &
                   [i_homog_areas, i_sg_correl])

      if ( .not. l_mcr_qrain ) then
        ! Error if trying to use prognostic precipitation fraction
        ! without prognostic rain
        ErrorStatus = 100
        comments = 'Cannot use prognostic precipitation fraction' // newline// &
                   '(l_mcr_precfrac)'                             // newline// &
                   'without also using prognostic rain-mass'      // newline// &
                   '(l_mcr_qrain).'
        call ereport(RoutineName, ErrorStatus, comments)
      end if

      if ( l_subgrid_graupel_frac .and.                                        &
           (sediment_loc==rain_sed_end .or. sediment_loc==warm_sed_end) ) then
        ! Error if trying to use combined prognostic fraction for rain and
        ! graupel, but doing sedimentation of rain and graupel at different
        ! points in the timestep
        ErrorStatus = 100
        comments = 'Cannot use combined prognostic fraction for'  // newline// &
                   'rain and graupel'                             // newline// &
                   '(l_mcr_precfrac .and. l_subgrid_graupel_frac)'// newline// &
                   'if doing sedimentation of rain and graupel at'// newline// &
                   'different points in the timestep'             // newline// &
                   '(sediment_loc == rain_sed_end, warm_sed_end).'
        call ereport(RoutineName, ErrorStatus, comments)
      end if

    end if  ! ( l_mcr_precfrac )

    if ( l_subgrid_graupel_frac .and.                                          &
         ( .not. ( graupel_option > no_graupel .and. l_mcr_precfrac ) ) ) then
      ! Error if trying to use sub-grid fraction for graupel
      ! without prognostic graupel and a prognostic precipitation fraction
      ErrorStatus = 100
      comments = 'Cannot use sub-grid fraction for graupel '      // newline// &
                 '(l_subgrid_graupel_frac) '                      // newline// &
                 'without also using prognostic graupel '         // newline// &
                 '(graupel_option > 0) '                          // newline// &
                 'and a prognostic precipitation fraction '       // newline// &
                 '(l_mcr_precfrac).'
      call ereport(RoutineName, ErrorStatus, comments)
    end if

    if ( l_proc_fluxes .and. (.not. sediment_loc==all_sed_start) ) then
      ! If applying microphysics process rates to fluxes as well as prognostics
      ! Error if not doing all the sedimentation calculations at the start
      ErrorStatus = 100
      comments = 'Applying process rates to precip fluxes '       // newline// &
                 '(l_proc_fluxes) requires that all precip '      // newline// &
                 'sedimentation calculations are done before '    // newline// &
                 'process rate calculations.  You must set '      // newline// &
                 'sediment_loc=1 (all_sed_start).'
      call ereport(RoutineName, ErrorStatus, comments)
    end if

  end if ! l_casim

  ! Check tapering values of cloud droplet number
  ! Valid for both CASIM and non-CASIM runs.
  if (l_droplet_tpr) then

    call chk_var(ndrop_surf, 'ndrop_surf', '[1.0E6:375.0E6]')
    call chk_var(z_surf, 'z_surf', '[0.0:150.0]')
    ! Checking that z_surf <= z_peak_nd now hard-wired to 150.0 m

  end if ! l_droplet_tpr

end if ! l_rain

def_src = ''
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine check_run_precip

end module mphys_inputs_mod

