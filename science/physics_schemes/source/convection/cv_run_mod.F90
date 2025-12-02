! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
!  Global data module for switches/options concerned with convection.

module cv_run_mod

  ! Description:
  !   Module containing runtime logicals/options used by the convection code.
  !
  ! Method:
  !   All switches/options which are contained in the &Run_Convection
  !   namelist in the CNTLATM control file are declared in this module.
  !   Default values have been declared where appropriate.
  !
  !   Any routine wishing to use these options may do so with the 'use'
  !   statement.
  !
  ! Code Owner: Please refer to the UM file CodeOwners.txt
  ! This file belongs in section: Convection
  !
  ! Code Description:
  !   Language: FORTRAN 90
  !   This code is written to UMDP3 programming standards.
  !
  ! Declarations:

use missing_data_mod, only: rmdi, imdi
use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim
use errormessagelength_mod, only: errormessagelength

use um_types, only: real_umphys

implicit none
save

! To satisfy the requirements of ROSE, the defaults for all these
! switches/variables have been set to false/MDI.
! The original defaults are noted in comments below.
! These values are also indicated in the ROSE help text.

integer :: i_convection_vn = imdi
                      ! Switch to determine version of convection scheme
                      ! 5 => 5A scheme
                      ! 6 => 6A scheme
                      ! 10 => Lambert-Lewis scheme
                      ! 11 => Simple Betts-Miller scheme
                      ! 12 => CoMorph flexible mass-flux scheme
integer, parameter :: i_convection_vn_5a = 5
integer, parameter :: i_convection_vn_6a = 6
integer, parameter :: i_cv_llcs = 10
integer, parameter :: i_cv_betts = 11
integer, parameter :: i_cv_comorph = 12

integer :: llcs_cloud_precip = imdi ! LLCS NL: cloud vs precip option

! LLCS Cloud vs precipitation options
integer, parameter :: llcs_opt_all_rain = 0
integer, parameter :: llcs_opt_all_cloud = 1
integer, parameter :: llcs_opt_crit_condens = 2
integer, parameter :: llcs_opt_const_frac = 3

!===========================================================================
! Logical switches set from GUI
!===========================================================================

! Comments for true status

logical :: l_fix_udfactor = .false. ! Fix application of UD_FACTOR (lock)

logical :: l_cloud_deep   = .false. ! Use depth criterion for convective anvils
                                    ! Original default true
logical :: l_mom          = .false. ! Use Convective Momentum Transport

logical :: l_eman_dd      = .false. ! Use Emanuel downdraught scheme

logical :: l_rediagnosis  = .false. ! If more than one convection step per model
                                    ! physics timestep then rediagnose cumulus
                                    ! points by recalling conv_diag before each
                                    ! convective sub-step.
                                    ! Original default true

logical :: l_anvil        = .false. ! Apply anvil scheme to 3D convective cloud.
                                    ! No effect unless l_3d_cca = .true.

logical :: l_snow_rain    = .false. ! If .true. allows a mix of snow and rain
                                    ! between the freezing temperature
                                    ! (273.15K) and t_melt_snow

logical :: l_dcpl_cld4pc2 = .false. ! Decouples section 0/5 cloud properties
                                    ! from each other, use with PC2.
                                    ! May change results if used without PC2.

logical :: l_murk_conv    = .false. ! Enables convective mixing of MURK as
                                    ! a tracer

logical :: l_safe_conv    = .false. ! Safer convection, remove negative q
                                    ! before attempting convection, don't
                                    ! add increments for ascents with
                                    ! negative CAPE. (5A/6A only)

logical :: l_ccrad        = .false. ! Main Logical, will include code
                                    ! connected with CCRad.
                                    ! (including bugfixes)

logical :: l_3d_cca       = .false. ! Use 3D convective cloud amount

logical :: l_conv_hist    = .false. ! True if 3 extra prognostics holding
                                    ! convective history information.

logical :: l_param_conv   = .false. ! Run time switch for convection scheme

logical :: l_conv_prog_dtheta  = .false.
                                    ! Use time-smoothed convective theta
                                    ! increments
logical :: l_conv_prog_dq      = .false.
                                    ! Use time-smoothed convective humidity
                                    ! increments
logical :: l_conv_prog_flx     = .false.
                                    ! Use time-smoothed mass flux in deep and
                                    ! mid closures
logical :: l_conv_prog_precip  = .false.
                                    ! Use convection 3D prognostic field of
                                    ! recent convective precipitation.
logical :: adv_conv_prog_dtheta   = .false.
                                    ! Sets whether or not to advect
                                    ! conv_prog_dtheta
logical :: adv_conv_prog_dq       = .false.
                                    ! Sets whether or not to advect
                                    ! conv_prog_dq
logical :: adv_conv_prog_flx      = .false.
                                    ! Sets whether or not to advect
                                    ! conv_prog_flx
logical :: l_prog_pert  = .false.   ! Use convection 3D prognostic field of
                                    ! recent convective precipitation to scale
                                    ! the initial buoyancy perturbation for
                                    ! mid-level convection
logical :: l_jules_flux = .false.
                                    ! Use the real surface fluxes calculated
                                    ! by Jules in the diagnosis rather
                                    ! than the simple diagnosis estimate
logical :: l_fcape        = .false. ! If .true. uses a forced detrainment
                                    ! weighting in the calculation of CAPE
                                    ! for CAPE closure. If not then the
                                    ! calculation is unweighted.
logical :: l_ccp_blv = .false.
                                    ! Logical to enable
                                    ! cold-pool BL velocity perturbation.
logical :: l_ccp_parcel_md = .false.
                                    ! Logical to enable cold-pool perturbation
                                    ! of the mid-level parcel buoyancy.
logical :: l_ccp_parcel_dp = .false.
                                    ! Logical to enable cold-pool perturbation
                                    ! of the deep parcel buoyancy.
logical :: l_ccp_trig = .false.
                                    ! Logical to enable cold-pool
                                    ! modification of the surface-flux based
                                    ! triggering of convection.
logical :: l_ccp_wind = .false.
                                    ! Switch to add ambient wind components
                                    ! to cold-pool propagation vectors.
logical :: l_ccp_seabreeze = .false.
                                    ! Switch to turn on sea-breeze
                                    ! representation in the cold-pool scheme.
logical :: l_reset_neg_delthvu = .false.
                                    ! Switch to suppress diagnosing cumulus
                                    ! where buoyancy integral comes out negative
logical :: l_cvdiag_ctop_qmax = .false.
                                    ! In the convective diagnosis, use revised
                                    ! check for well-mixed q-profile,
                                    ! accounting for local maxima at cloud-top

!===========================================================================
! Logical switches not set from GUI
!===========================================================================

logical :: l_cv_conserve_check = .false.
                                    ! Diagnostic conservation checks. Goes
                                    ! through energy correction code
                                    ! (5A/6A only)

logical :: l_cmt_heating  = .false. ! Include the heating due to the
                                    ! dissipation of kinetic energy
                                    ! from the convective momentum transport

logical :: l_mr_conv = .false.      ! Flag for convection running with
                                    ! mixing ratios rather than specifics.

logical :: l_wvar_for_conv = .false.! Flag to calculate boundary-layer
                                    ! turbulent vertical velocity variance
                                    ! for use by the convection scheme.

!===========================================================================
! Integer options set from GUI
!===========================================================================

! Convection integer options set from gui/namelist Convection Scheme

integer :: n_conv_calls       = imdi ! Number of calls to convection
                                     ! per physics timestep
                                     ! Original default 1

integer :: cld_life_opt       = imdi ! Convective cloud decay time
                                     ! Original default cld_life_constant (0)
integer :: rad_cloud_decay_opt= imdi ! Convective cloud decay
                                     ! Original default rad_decay_off (0)
integer :: anv_opt            = imdi ! Anvil cloud basis
                                     ! Original default anv_model_levels (2)

! Convection Scheme Options (5A)
! NOTE: These options were valid at the time of writing. They are used in
!       development code(5A) and so very likely to change.
!       Users should consult the Convection Group for current available
!       options.

integer :: iconv_shallow  = imdi  ! Shallow (Original default 0)
                           !   0: no scheme,
                           !   1: Gregory-Rowntree scheme
                           !   2: Turbulence scheme (non-precipitating)
                           !   3: Turbulence scheme (precipitating)

integer :: iconv_congestus = imdi ! Congestus (Original default 0)
                           !   0: no scheme,
                           !   1: Gregory-Rowntree scheme untested for possible
                           !      future use.

integer :: iconv_mid       = imdi ! Mid-level (Original default 0)
                           !   0: no scheme,
                           !   1: Gregory-Rowntree scheme
                           !   2: Future use

integer :: iconv_deep      = imdi ! Deep (Original default 0)
                           !   0: no scheme,
                           !   1: Gregory-Rowntree scheme
                           !   2: turbulence scheme

integer :: deep_cmt_opt    = imdi ! Deep CMT tunings (Original default 0)
                           !   0: turbulence based
                           !   1: Operational 70 level (modified turb based)
                           !   2: Gregory-Kershaw scheme
                           !   3: New turbulence based scheme.
                           !   4: Future use
                           !   6: Stabilized Gregory-Kershaw scheme

integer :: mid_cmt_opt     = imdi ! Mid CMT scheme to be used (Orig default 0)
                           !   0: Gregory-Kershaw scheme
                           !   1: Diffusive scheme
                           !   2: Stabilized Gregory-Kershaw scheme

integer :: icvdiag    = imdi ! Diagnosis calculation options (Orig default 1)
                           !   0: No longer allowed
                           !   1: 5A/6A scheme default
                           !   2: 0.55/z entrainment rates
                           !   3: 1/z entrainment
                           !   4: Diurnal cycle over land entrainment rates
                           !   5: Diurnal cycle over land entrainment rates,
                           !      Ocean undilute ascent.
                           ! 6A scheme only
                           !   6: dilute, p/(p*)^2 entrainment rate dependent
                           !      on precipitation based 3d convective
                           !      prognostic at the initiation level, as
                           !      ent_opt_dp=6
                           !   7: dilute, p/(p*)^2 entrainment rate dependent
                           !      on precipitation based 3d convective
                           !      prognostic at the current level, as
                           !      ent_opt_dp=7
                           !   8: fixed p/(p*)^2 entrainment rate, as
                           !      ent_opt_dp=0

integer :: cvdiag_inv = imdi ! Inversion test in convective diagnosis? 5A & 6A
                           ! Original default 1
                           ! When doing an undilute parcel ascent
                           !   0: No inversion test
                           !   1: Original inversion test (Default 5A and 6A)
                           !   2: Future alternative inversion tests

integer :: tv1_sd_opt = imdi ! Standard dev of level virtual temperature options
                           ! Original default 0
                           !   0: Assume BL depth of 300m (Default)
                           !   1: Use calculated BL depth
                           !   2: As (1) plus stability dependence and coastal
                           !      mean

integer :: midtrig_opt = imdi ! Mid-level triggering option. 6A only.
                           !   0: Allow mid-level convection from two levels
                           !      above ntml or two levels above
                           !      shallow or deep termination
                           !   1: Allow mid-level convection from the level
                           !      above ntml or the level above
                           !      shallow or deep termination
                           !   2: Allow mid-level convection from the level
                           !      1 or the level above
                           !      shallow or deep termination

integer :: adapt    = imdi ! Adaptive detrainment/entrainment options
                           ! Original default 0
                           !   0: Original (Default)
                           !   1: Adaptive detrainment: mid and deep convection
                           !   2: Future/Experimental (En/Detrainment)
                           !   3: Adaptive detrainment: deep convection
                           !   4: Adaptive detrainment: shallow, mid and
                           !      deep convection
                           !   5: Smoothed adaptive detrainment:
                           !      mid and deep convection
                           !   6: Smoothed adaptive detrainment:
                           !      shallow, mid and deep convection
                           !   7: Improved smoothed adaptive detrainment:
                           !      mid and deep convection
                           !   8: Improved smoothed adaptive detrainment:
                           !      shallow, mid and deep convection

integer :: fdet_opt = imdi ! Forced detrainment calculation options. 6A only.
                           !   0: calculate forced detrainment rate using
                           !      humidity continuity equation.
                           !   1: calculate forced detrainment rate using
                           !      theta continuity equation.
                           !   2: forced detrainment rate based on theta
                           !      equation with an improved treatment of
                           !      subsaturated conditions.
                           !   3: forced detrainment rate based on theta
                           !      equation with improved treatment of
                           !      subsaturated conditions and improved
                           !      numerical calculation of detrained parcel
                           !      properties.

integer :: ent_opt_dp = imdi  ! Deep entrainment option (Orig default 0)
                           !   0: original Ap/(p*)^2
                           !   1: n/z dependence where n=ent_fac_dp
                           !   2: As Ap/(p*)^2 but multiplied by extra factor
                           !      f=1.+3.*(1.-(100000.-p(k))/(100000.-50000.))
                           !   3: factor * (A/p*)*((p/p*)^ent_dp_power)
                           !   Options 4 & 5 are experimental options that
                           !   depend on the diagnosed depth of convection.
                           !   Available for 6a only.
                           !   4: variable n/z style entrainment
                           !   5: variable p/(p*)^2  style entrainment
                           !   6: variable p/(p*)^2 entrainment dependent on
                           !     precipitation based 3d convective prognostic
                           !     at the initiation level.
                           !   7: variable p/(p*)^2 entrainment dependent on
                           !     precipitation based 3d convective prognostic
                           !      at the current level.
integer :: ent_opt_md = imdi  ! mid entrainment option (Orig default 0)
                           !   0: original Ap/(p*)^2
                           !   1: n/z dependence where n=ent_fac_md
                           !   2: As Ap/(p*)^2 but multiplied by extra factor
                           !      f=1.+3.*(1.-(100000.-p(k))/(100000.-50000.))
                           !   3: factor * (A/p*)*((p/p*)^ent_md_power)
                           !   6: variable p/(p*)^2 entrainment dependent on
                           !     precipitation based 3d convective prognostic
                           !     at the initiation level.
                           !   7: variable p/(p*)^2 entrainment dependent on
                           !     precipitation based 3d convective prognostic
                           !      at the current level.
integer :: mdet_opt_dp = imdi  ! Deep mixing detrainment option - 6a only
                           !   0: orig_amdet_fac*original entrainment/3
                           !   1: amdet_fac*entrainment*(1-rh)
integer :: mdet_opt_md = imdi  ! mid-level mixing detrainment option - 6a only
                           !   0: orig_amdet_fac*original entrainment/3
                           !   1: amdet_fac*entrainment*(1-rh)

! Deep cloud base closure options 5A & 6A
integer :: cldbase_opt_dp = imdi
                           ! This option may have one of the following values:
                           !   0: RH based timescale
                           !   1: RH based timescale (timestep limited)
                           !   2: Fixed timescale
                           !   3: Fixed timescale reduced if w exceeds w max
                           !      CAPE closure
                           !   4: Area scaled CAPE closure
                           !   5: Experimental - not allowed
                           !   6: RH based timescale (timestep limited)
                           !      plus reduction if w exceeds w-max
                           !   7: CAPE timescale dependent on large-scale
                           !      vertical velocity
                           !   8: Closure based on boundary layer fluxes
                           !      and large-scale vertical velocity based CAPE
                           !      timescale
                           !   9: Closure based on boundary layer fluxes
                           !      and large-scale vertical velocity.

integer :: cldbase_opt_md = imdi  ! Mid-level cloud base closure options 5A&6A
                           !   0: RH based timescale
                           !   1: RH based timescale (timestep limited)
                           !   2: Fixed timescale
                           !   3: Fixed timescale reduced if w exceeds w max
                           !      CAPE closure
                           !   4: Area scaled CAPE closure
                           !   5: Experimental - not allowed
                           !   6: RH based timescale (timestep limited)
                           !      plus reduction if w exceeds w-max
                           !   7: CAPE timescale dependent on large-scale
                           !      vertical velocity

integer :: cldbase_opt_sh = imdi  ! Shallow cloud base closure options 5A&6A
                           !   0: Closure based on boundary layer fluxes
                           !   1: Grey zone combination of fixed timescale
                           !      CAPE closure and BL fluxes

integer :: cape_bottom = imdi    ! Start level for w_max in column
                                 ! Original default IMDI

integer :: cape_top    = imdi    ! End   level for w_max in column
                                 ! Original default IMDI

integer :: sh_pert_opt  = imdi   ! Initial perturbation method for
                                 ! shallow cumulus (Orig default 0)
                                 ! 0 = Original code
                                 ! 1 = Revised  code


integer :: md_pert_opt  = imdi   ! Initial perturbation method for
                                 ! mid-level convection Orig default 0)
                                 ! 0 = Original code
                                 ! 1 = specify the ratio of temperature
                                 !     and humidity perturbation via efrac
                                 ! 2 = initial parcel perturbation based on
                                 !     the surface evaporative fraction

integer :: limit_pert_opt = imdi ! Limits convective parcel perturbation
                                 ! to physically sensible values.
                                 ! Orig default 0
                                 ! 0 = original code - no limits
                                 ! 1 = apply limits to main ascent only
                                 ! 2 = apply limits to main ascent and in
                                 !     the convection diagnosis


integer :: dd_opt         = imdi ! Downdraught scheme options
                                 ! Orig default 0
                                 ! 0 = Original code
                                 ! 1 = Revised  code

integer :: termconv       = imdi ! Original default 0
                                 ! 0 for default
                                 ! 1 for modified termination condition

integer :: bl_cnv_mix     = imdi ! Options for mixing convective increments in the BL
                                 ! Original default 0
                                 ! 0: original code
                                 ! 1: only mix the increments from the initial
                                 !    parcel perturbation

integer :: cnv_wat_load_opt=imdi ! Options for including liquid and frozen water loading
                                 ! in the convective updraught buoyancy calculation
                                 ! Original default 0
                                 ! 0: Do not include water loading (default)
                                 ! 1: Include water loading

integer :: cca2d_sh_opt = imdi   ! Method to evaluate cca2d (Shallow)
                                 ! Original default 0
integer :: cca2d_md_opt = imdi   ! Method to evaluate cca2d (Mid-level)
                                 ! Original default 0
integer :: cca2d_dp_opt = imdi   ! Method to evaluate cca2d (Deep)
                                 ! Original default 0

integer :: ccw_for_precip_opt = imdi  ! Option controlling critical cloud water
                                      ! for the formation of precipitation
                                      ! Original default 0
                                      ! 0 - original code (dependent on
                                      !     fac_qsat & qlmin) Land sea
                                      !     dependence with dcrit.
                                      ! 1 - no Dcrit option
                                      ! 2 - Known as Manoj's first function
                                      ! 3 - Known as Manoj's congestus function
                                      ! 4 - no land sea diff dependence on
                                      !     fac_qsat & qlmin.

integer :: plume_water_load = imdi ! Option for water loading in undilute
                                   ! parcel ascent
                                   ! Original default 0
                                   ! 0 - no water removal original undilute
                                   ! 1 - remove any water > 1g/kg
                                   ! 2 - remove any water > profile varying
                                   !     with qsat

integer :: dil_plume_water_load  = imdi ! Option for water loading in dilute
                                        ! parcel ascent
                                        ! Original default 0
                                        ! 0 - no water removal
                                        ! 1 - remove any water > 1g/kg
                                        ! 2 - remove any water > profile varying
                                        !     with qsat

integer :: cnv_cold_pools = imdi ! Convective Cold Pool scheme options
                                 ! 0 - scheme not used
                                 ! 1 - CCP version 1 : propagating
                                 ! 2 - CCP version 2 : single-column

integer :: pr_melt_frz_opt = imdi  ! Options for the treatment of phase
                                   ! changes in convective precipitation
                                   ! 0 - Rain instantly freezes if T<0C
                                   !     and snow instantly melts if T>0C
                                   ! 1 - Allows mixed phase precipitation.
                                   !     No freezing of rain. At T>0C Snow is
                                   !     melted at a rate that depends on the
                                   !     temperature and t_melt_snow.
                                   ! 2 - Allows mixed phase precipitation.
                                   !     At T<0C rain is frozen at a rate
                                   !     that depends on the temperature and
                                   !     t_frez_rain (hardwired).
                                   !     At T>0C Snow is melted at a rate that
                                   !     depends on the temperature and
                                   !     t_melt_snow as for option 1.

!===========================================================================
! Real values set from GUI
!===========================================================================

real(kind=real_umphys) :: cca_sh_knob = rmdi
                               ! Scales Shallow cloud fraction (CCRad)
                               ! Original default 1.0
real(kind=real_umphys) :: cca_md_knob = rmdi
                               ! Scales Mid     cloud fraction (CCRad)
                               ! Original default 1.0
real(kind=real_umphys) :: cca_dp_knob = rmdi
                               ! Scales Deep    cloud fraction (CCRad)
                               ! Original default 1.0
real(kind=real_umphys) :: ccw_sh_knob = rmdi
                               ! Scales Shallow cloud water (CCRad)
                               ! Original default 1.0
real(kind=real_umphys) :: ccw_md_knob = rmdi
                               ! Scales Mid     cloud water (CCRad)
                               ! Original default 1.0
real(kind=real_umphys) :: ccw_dp_knob = rmdi
                               ! Scales Deep    cloud water (CCRad)
                               ! Original default 1.0

real(kind=real_umphys) :: fixed_cld_life = rmdi
                              ! Fixed convective cloud lifetime decay value
                              ! (seconds) (Original default 7200.0)

real(kind=real_umphys) :: cca_min = rmdi
                        ! Threshold value of convective cloud fraction
                        ! below which cca has neglible radiative impact and
                        ! is reset to zero (Original default 0.02)

real(kind=real_umphys) :: r_det   = rmdi
                        ! Parameter controlling adaptive detrainment -
                        ! Orig default 0.75 (operational)
                        ! HadGEM1a recommended 0.5

real(kind=real_umphys) :: cape_min     = rmdi ! Scale dependent min cape
                            ! Original default RMDI
real(kind=real_umphys) :: w_cape_limit = rmdi
                            ! Test w for scale dependent cape timescale
                            ! Original default RMDI
real(kind=real_umphys) :: cape_ts_min  = rmdi
                            ! Minimum CAPE timescale, for w-base CAPE closure
real(kind=real_umphys) :: cape_ts_max  = rmdi
                            ! Maximum CAPE timescale, for w-base CAPE closure
real(kind=real_umphys) :: c_mass_sh       = rmdi
                               ! Cloud-base area scaling applied to shallow
                               ! convection mass-flux (original default 0.03).

real(kind=real_umphys) :: mparwtr      = rmdi
                     ! Maximum value of the function that is used to calculate
                     ! the maximum convective cloud water/ice in a layer (kg/kg)
                     ! Original default 1.0e-3
real(kind=real_umphys) :: qlmin        = rmdi
                     ! Minimum value of the function that is used to calculate
                     ! the maximum convective cloud water/ice in a layer (kg/kg)
                     ! Original default 2.0e-4
real(kind=real_umphys) :: fac_qsat     = rmdi
                            ! Factor used to scale environmental qsat to give
                            ! the maximum convective cloud water/ice in a layer
                            ! Original default RMDI
real(kind=real_umphys) :: mid_cnv_pmin = rmdi
                            ! The minimum pressure (max height) at which mid
                            ! level convection is allowed to initiate (Pa)
                            ! Original default 0.0
real(kind=real_umphys) :: orig_mdet_fac= rmdi
                            ! Factor multiplying the original mixing
                            ! detrainment rate.
                            ! Original default 1.0
real(kind=real_umphys) :: amdet_fac    = rmdi
                            ! Factor multiplying (1-rh) for adaptive mixing
                            ! detrainment rate.
                            ! Original default 1.0
! Scales with cca2d to determine convective cloud amount with anvil
real(kind=real_umphys) :: anvil_factor = rmdi
                            ! x cca2d = max cca, cca is capped at 1.0
                            ! but this does not imply a cap to anvil_factor
                            ! Original default RMDI
real(kind=real_umphys) :: tower_factor = rmdi ! x cca2d = min cca
                            ! Original default RMDI
real(kind=real_umphys) :: ud_factor    = rmdi
                            ! Updraught factor used in calculation of
                            ! convective water path Original default RMDI

real(kind=real_umphys) :: tice         = rmdi
                            ! Phase change temperature in plume
                            ! Original default 273.15
real(kind=real_umphys) :: qstice       = rmdi ! Qsat at phase change temperature
                            ! (freeze/melting temperature)
                            ! Original default 3.5E-3

real(kind=real_umphys) :: t_melt_snow     = rmdi
                               ! Temperature at which to melt all snow in
                               ! the downdraught.

! 5A & 6A code only
real(kind=real_umphys) :: ent_fac_sh      = rmdi
                               ! Factor multiplying entrainment rate - shallow
                               ! Original default 1.0
real(kind=real_umphys) :: ent_fac_dp      = rmdi
                               ! Factor multiplying entrainment rate - deep
                               ! Original default 1.0
real(kind=real_umphys) :: ent_fac_md      = rmdi
                               ! Factor multiplying entrainment rate - mid-level
                               ! Original default 1.0
real(kind=real_umphys) :: ent_dp_power    = rmdi
                               ! Power n for (p/p*)^n for entrainment option 3
                               ! Original default 2.0
real(kind=real_umphys) :: ent_md_power    = rmdi
                               ! Power n for (p/p*)^n for entrainment option 3
                               ! Original default 2.0
real(kind=real_umphys) :: cape_timescale  = rmdi ! Timescale for CAPE closure.
                               ! Original default RMDI
real(kind=real_umphys) :: cvdiag_sh_wtest = rmdi
                               ! w for air above shallow convection must be
                               ! less than this. (Default value is 0.0)
                               ! Original default 0.0
real(kind=real_umphys) :: cpress_term     = rmdi
                               ! Coefficient for pressure term relative to
                               ! shear term as found in Kershaw & Gregory 1997

! 6A code only
real(kind=real_umphys) :: eff_dcfl        = rmdi
                               ! Factor that defines the efficiency by which
                               ! detrained liquid condensate creates liquid
                               ! cloud fraction
real(kind=real_umphys) :: eff_dcff        = rmdi
                               ! Factor that defines the efficiency by which
                               ! detrained frozen condensate creates frozen
                               ! cloud fraction
real(kind=real_umphys) :: efrac           = rmdi
                               ! The ratio of latent heat flux to sensible
                               ! plus latent heat flux owing to the initial
                               ! parcel perturbation for mid-level convection.
real(kind=real_umphys) :: max_pert_scale  = rmdi
                               ! The maximum prognostic dependent mid-level
                               ! initial perturbation scaling if l_prog_pert=T
real(kind=real_umphys) :: thpixs_mid      = rmdi
                               ! The initial excess potential temperature (k)

! Convective prognostic options / variables
real(kind=real_umphys) :: tau_conv_prog_precip = rmdi
                               ! decay timescale for prognostic field of
                               ! recent convective precipitation
real(kind=real_umphys) :: tau_conv_prog_dtheta = rmdi
                               ! decay timescale for time-smoothed
                               ! convective theta increment
real(kind=real_umphys) :: tau_conv_prog_dq     = rmdi
                               ! decay timescale for time-smoothed
                               ! convective humidity increment
real(kind=real_umphys) :: tau_conv_prog_flx    = rmdi
                               ! decay timescale for time-smoothed
                               ! convective mass flux

! Factors used to calculate entrainment scaling from 3d prognostic field
! based on surface precipitation. Used under ent_opt_md,ent_opt_dp=6 or 7
real(kind=real_umphys) :: prog_ent_grad   = rmdi ! gradient term
real(kind=real_umphys) :: prog_ent_int    = rmdi ! intercept term
real(kind=real_umphys) :: prog_ent_max    = rmdi ! maximum scaling
real(kind=real_umphys) :: prog_ent_min    = rmdi ! minimum scaling

! Convective cold pool options

real(kind=real_umphys) :: kappa_g    = rmdi
! g'_DD tuning parameter
real(kind=real_umphys) :: kappa_h    = rmdi
! h_DD tuning parameter
real(kind=real_umphys) :: remain_max = rmdi
! Limit on size of remain counters
real(kind=real_umphys) :: phi_ccp    = rmdi
! Scaling constant for convective cold-pool vertical kinetic energy scale
real(kind=real_umphys) :: ccp_nonelastic = rmdi
! non-elasticity of cold-pool collisions
real(kind=real_umphys) :: entrain_max = rmdi
! maximum value of cold-pool dependent entrainment rate
real(kind=real_umphys) :: entrain_min = rmdi
! minimum value of cold-pool dependent entrainment rate
real(kind=real_umphys) :: ccp_buoyancy = rmdi
! factor to control cold-pool contribution to parcel buoyancy

!===========================================================================
! Real fields (not set from GUI)
!===========================================================================

! Convective cold pool geometry-based field, set during initialisation
!
real(kind=real_umphys) :: del_max = rmdi
! gridbox max dimension for calculation purposes (m)

!------------------------------------------------------------------------------
! Switches added for Warm TCS scheme
!------------------------------------------------------------------------------
! The following have been removed from the RUN_convection name list
!  - defaults are hard wired

integer :: iwtcs_diag1 = 0  ! Options for WTCS diagnosis
integer :: iwtcs_diag2 = 0  ! Options for WTCS diagnosis

!------------------------------------------------------------------------------
! Further changes

real(kind=real_umphys), parameter :: qmin_conv = 1.0e-8
                                          ! Minimum allowed value of q after
   !  convection, also used for negative q check, and used in UKCA to ensure
   !  that stratospheric water vapour does not trigger messages from convection

!Betts namelist variables
real(kind=real_umphys) :: betts_rh  = rmdi    ! critical RH
real(kind=real_umphys) :: betts_tau = rmdi    ! relaxation time scale
real(kind=real_umphys) :: betts_cld_evap_const = rmdi
                                    ! re-evaporation scalar constant

! LLCS namelist
real(kind=real_umphys) :: llcs_rhcrit = rmdi ! critical RH
real(kind=real_umphys) :: llcs_timescale = rmdi ! relaxation time scale
real(kind=real_umphys) :: llcs_detrain_coef = rmdi ! detrainment coef
real(kind=real_umphys) :: llcs_rain_frac = rmdi ! Rain fraction


!------------------------------------------------------------------------------
! Define namelist &Run_Convection read in from CNTLATM control file.
! Changes made to this list will affect both the Full UM and the SCM
!------------------------------------------------------------------------------

namelist/Run_Convection/                                                       &

i_convection_vn,                                                               &

! Logical switches
l_mom, l_fix_udfactor,                                                         &
l_snow_rain,                                                                   &
l_eman_dd, l_cloud_deep,                                                       &
          l_rediagnosis,  l_dcpl_cld4pc2,                                      &
l_anvil,  l_murk_conv, l_safe_conv, l_cv_conserve_check,                       &
l_ccrad,              l_3d_cca,               l_conv_hist,                     &
l_param_conv,                                                                  &
l_cmt_heating,                                                                 &
l_conv_prog_dtheta,   l_conv_prog_dq,         l_conv_prog_flx,                 &
l_conv_prog_precip,                                                            &
adv_conv_prog_dtheta, adv_conv_prog_dq,       adv_conv_prog_flx,               &
l_prog_pert,          l_jules_flux,         l_fcape,                           &
l_ccp_blv, l_ccp_parcel_md, l_ccp_parcel_dp, l_ccp_trig,                       &
l_ccp_wind, l_ccp_seabreeze,                                                   &
l_reset_neg_delthvu, l_cvdiag_ctop_qmax,                                       &

! General scheme options/variables
n_conv_calls,                                                                  &
sh_pert_opt,          md_pert_opt,                                             &
dd_opt,               deep_cmt_opt,           mid_cmt_opt,                     &
termconv,             adapt,                  fdet_opt,                        &
r_det,                                                                         &
tice,                 qstice,                                                  &
ent_fac_sh,           ent_fac_dp,             ent_fac_md,                      &
ent_opt_dp,           ent_opt_md,                                              &
ent_dp_power,         ent_md_power,           efrac,                           &
max_pert_scale,       thpixs_mid,mdet_opt_dp,            mdet_opt_md,          &
bl_cnv_mix,                                                                    &
mid_cnv_pmin,         amdet_fac,              orig_mdet_fac,                   &
ccw_for_precip_opt,                                                            &
cnv_wat_load_opt,     tv1_sd_opt,             limit_pert_opt,                  &
cnv_cold_pools,       midtrig_opt,                                             &

! Convective diagnosis options
icvdiag,              plume_water_load,       dil_plume_water_load,            &
cvdiag_inv,           cvdiag_sh_wtest,                                         &

! Closure & Cape related options/variables
cldbase_opt_dp,       cldbase_opt_md,         cldbase_opt_sh,                  &
cape_bottom,          cape_top,               cape_timescale,                  &
w_cape_limit,         cape_min,               cape_ts_min,                     &
cape_ts_max,          c_mass_sh,                                               &

! Convective cloud options/variables
cld_life_opt,         rad_cloud_decay_opt,    cca_min,                         &
fixed_cld_life,       ud_factor,              mparwtr,                         &
qlmin,                fac_qsat,               eff_dcfl,                        &
eff_dcff,                                                                      &

! CMT variables
cpress_term,                                                                   &

! CCRad options options/variables
cca2d_sh_opt,         cca_sh_knob,            ccw_sh_knob,                     &
cca2d_md_opt,         cca_md_knob,            ccw_md_knob,                     &
cca2d_dp_opt,         cca_dp_knob,            ccw_dp_knob,                     &

! Anvil scheme options
anvil_factor,         tower_factor,           anv_opt,                         &

! Downdraught and precipitation phase change options and variables
pr_melt_frz_opt,      t_melt_snow,                                             &

! Convection type options
iconv_shallow,        iconv_mid,              iconv_deep,                      &
iconv_congestus,                                                               &

! Convection prognostic options / variables
tau_conv_prog_precip,                                                          &
tau_conv_prog_dtheta,  tau_conv_prog_dq, tau_conv_prog_flx,                    &

! Factors used to calculate entrainment scaling from 3d prognostic field
prog_ent_grad,        prog_ent_int,           prog_ent_max,                    &
prog_ent_min,                                                                  &

! Convective cold pool options
kappa_g,              kappa_h,                remain_max,                      &
phi_ccp,              ccp_buoyancy,           ccp_nonelastic,                  &
entrain_max,          entrain_min,                                             &

!free parameters for Bett-Miller scheme
betts_tau, betts_rh, betts_cld_evap_const,                                     &

! Parameters for Lambert-Lewis Convection Scheme
llcs_cloud_precip, llcs_detrain_coef, llcs_rhcrit,                             &
llcs_timescale, llcs_rain_frac

!------------------------------------------------------------------------------

!DrHook-related parameters
integer(kind=jpim), parameter, private :: zhook_in  = 0
integer(kind=jpim), parameter, private :: zhook_out = 1

character(len=*), parameter, private :: ModuleName='CV_RUN_MOD'

contains

subroutine check_run_convection()

use chk_opts_mod, only: chk_var, def_src
use cv_param_mod, only: ccp_off, ccp_prop, ccp_sc
use ereport_mod,  only: ereport
use timestep_mod, only: timestep

implicit none

integer :: icode       ! error code
integer, parameter :: len_cond = 10  ! length of condition string

character(len=*), parameter :: RoutineName='CHECK_RUN_CONVECTION'
character (len=errormessagelength) :: cmessage  ! used for ereport
character (len=errormessagelength) :: mymessage ! Addtional error message
character (len=len_cond) :: chkcond             ! Condition for check

real(kind=jprb)                    :: zhook_handle

!---------------------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!---------------------------------------------------------------------------
! required by chk_var routine
def_src = ModuleName//':'//RoutineName

! Which convection scheme

call  chk_var(i_convection_vn,'i_convection_vn',                               &
              [i_convection_vn_5a, i_convection_vn_6a,                         &
               i_cv_llcs, i_cv_betts, i_cv_comorph])

!---------------------------------------------------------------------------
! Check convective diagnosis options even when running without convection
!---------------------------------------------------------------------------

call chk_var(icvdiag,'icvdiag','[1:8]')
if (i_convection_vn /= i_convection_vn_6a) then
  ! Check that not trying to use a 6A option
  call chk_var(icvdiag,'icvdiag','[1:5]')
end if
! Don't allow using dilute diagnosis options designed to couple to certain
! 5A/6A convection scheme deep closure options, if not using the 5A/6A scheme!
if ( icvdiag > 3 .and. ( .not. ( l_param_conv .and.                            &
                                 ( i_convection_vn == i_convection_vn_6a .or.  &
                                   i_convection_vn == i_convection_vn_5a )     &
                               ) ) ) then
  icode = 10
  write(cmessage,"(A,I2,A)") "Convective diagnosis option icvdiag = ", icvdiag,&
                             " can only be used in conjunction with" //        &
                             " parameterised convection using the" //          &
                             " 5A or 6A convection scheme."
  call ereport(RoutineName,icode,cmessage)
end if

call chk_var(cvdiag_inv,'cvdiag_inv','[0,1]')

call chk_var(cvdiag_sh_wtest,'cvdiag_sh_wtest','[-1.0:1.0]')

call chk_var(limit_pert_opt,'limit_pert_opt','[0,1,2]')
call chk_var(tv1_sd_opt,'tv1_sd_opt','[0,1,2]')

! These are for plume water loading in the diagnosis parcel:
call chk_var(plume_water_load,'plume_water_load','[0,1,2]')
call chk_var(dil_plume_water_load,'dil_plume_water_load','[0,1,2]')

!---------------------------------------------------------------------------
! Only check the rest of the namelist if running with convection on
!---------------------------------------------------------------------------
if (l_param_conv) then

  ! Options which only apply to the 5A and 6A schemes...
  if (i_convection_vn == i_convection_vn_5a .or.                               &
      i_convection_vn == i_convection_vn_6a) then
    ! Sub-stepping - don't want to encourage a high value so now limiting
    ! to 4
    mymessage='Calling convection more than 4 times per timestep is not allowed'
    call chk_var(n_conv_calls,'n_conv_calls','[1:4]', cmessage=mymessage)

    ! Values only applying to 6A scheme
    if (i_convection_vn == i_convection_vn_6a) then

      ! Shallow convection cloud-base closure scaling and entrainment factor
      call chk_var(c_mass_sh,'c_mass_sh','[0.01:0.09]')
      call chk_var(ent_fac_sh,'ent_fac_sh','[0.33:3.0]')

      ! Scaling for deep and mid-level mixing detrainment (original scheme)
      call chk_var(orig_mdet_fac,'orig_mdet_fac','[0.0:5.0]')

      ! l_cmt_heating
      call chk_var(mdet_opt_dp,'mdet_opt_dp','[0:1]')
      call chk_var(mdet_opt_md,'mdet_opt_md','[0:1]')
      call chk_var(fdet_opt,'fdet_opt','[0:3]')
      call chk_var(eff_dcff,'eff_dcff','[0.0:10.0]')
      call chk_var(eff_dcfl,'eff_dcfl','[0.0:10.0]')

      ! Convective memory
      if (l_conv_prog_precip) then
        ! Note meta data had no upper bound but checking routines require one
        ! so put in a day in seconds as not expecting to use values as large
        ! as this.
        call chk_var(tau_conv_prog_precip,'tau_conv_prog_precip',              &
                     '[1.0:86400.0]')
        call chk_var(prog_ent_grad,'prog_ent_grad','[-3.0:0.0]')
        call chk_var(prog_ent_int,'prog_ent_int','[-10.0:5.0]')
        call chk_var(prog_ent_max,'prog_ent_max','[0.0:5.0]')
        call chk_var(prog_ent_min,'prog_ent_min','[0.0:5.0]')
        call chk_var(w_cape_limit,'w_cape_limit','[0.0:10000.0]')
      end if

      !Mid-level convection only options
      if (iconv_mid == 1) then
        call chk_var(midtrig_opt,'midtrig_opt','[0:3]')
        call chk_var(md_pert_opt,'md_pert_opt','[0,1,2]')
        call chk_var(thpixs_mid, 'thpixs_mid', '[0.0:2.0]')
        if (md_pert_opt >= 1) then
          call chk_var(efrac,'efrac','[-10.0:10.0]')
        end if
        if (l_prog_pert) then
          call chk_var(max_pert_scale,'max_pert_scale','[1.0:20.0]')
        end if
      end if

      if (timestep > 0.0) then
        !Only check timescales if timestep is defined
        if (l_conv_prog_dtheta) then
          !Time smoothed theta increments
          !timescale ge model timestep
          write(chkcond, '(A3,F6.1,A1)')'[>=',timestep,']'
          call chk_var(tau_conv_prog_dtheta,'tau_conv_prog_dtheta',chkcond)
        end if
        if (l_conv_prog_dq) then
          !Time smoothed humidity increments
          !timescale ge model timestep
          write(chkcond, '(A3,F6.1,A1)')'[>=',timestep,']'
          call chk_var(tau_conv_prog_dq,'tau_conv_prog_dq',chkcond)
        end if
        if (l_conv_prog_flx) then
          !Time smoothed convective mass flux
          !timescale ge model timestep
          write(chkcond, '(A3,F6.1,A1)')'[>=',timestep,']'
          call chk_var(tau_conv_prog_flx,'tau_conv_prog_flx',chkcond)
        end if
      end if

      ! CMT checking
      if (l_mom) then
        if (mid_cmt_opt  == 1 .or. mid_cmt_opt == 3 .or.                       &
            deep_cmt_opt == 2 .or. deep_cmt_opt ==6) then
          !cpress_term only needed with Gregory-Kershaw CMT
          call chk_var(cpress_term,'cpress_term','[0.0:1.0]')
        end if
      end if

      ! Convective cold-pool scheme
      if ((cnv_cold_pools == ccp_prop) .or. (cnv_cold_pools == ccp_sc)) then
        call chk_var(ccp_buoyancy,'ccp_buoyancy','[0.0:50.0]')
        call chk_var(entrain_max,'entrain_max','[0.0:50.0]')
        call chk_var(entrain_min,'entrain_min','[0.0:50.0]')
        call chk_var(kappa_g,'kappa_g','[0.0:1.0e7]')
        call chk_var(kappa_h,'kappa_h','[0.0:1.0e7]')
        call chk_var(phi_ccp,'phi_ccp','[0.0:1.0]')
      end if
      if (cnv_cold_pools == ccp_prop) then
        call chk_var(ccp_nonelastic,'ccp_nonelastic','[0.0:1.0]')
        call chk_var(remain_max,'remain_max','[1.0:100.0]')
      end if

    end if  ! (i_convection_vn == i_convection_vn_6a)

    ! General convection panel
    ! Allowing setting of 1 to be used for testing mass flux congestus code
    call chk_var(iconv_congestus,'iconv_congestus','[0,1]')
    call chk_var(iconv_deep,'iconv_deep','[0:2]')
    call chk_var(iconv_mid,'iconv_mid','[0,1]')
    call chk_var(iconv_shallow,'iconv_shallow','[0:3]')

    call chk_var(bl_cnv_mix,'bl_cnv_mix','[0,1]')
    call chk_var(cnv_wat_load_opt,'cnv_wat_load_opt','[0,1]')
    call chk_var(mid_cnv_pmin,'mid_cnv_pmin','[0.0:1.0e7]')

    call chk_var(qlmin,'qlmin','[0.0:1.0e-3]')
    call chk_var(qstice,'qstice','[0.0:0.5]')
    call chk_var(tice,'tice','[200.0:300.0]')

    ! Updraughts panel

    call chk_var(sh_pert_opt,'sh_pert_opt','[0,1]')
    if (i_convection_vn == i_convection_vn_6a) then
      call chk_var(termconv,'termconv','[0:2]')
    else
      call chk_var(termconv,'termconv','[0,1]')
    end if

    call chk_var(ccw_for_precip_opt,'ccw_for_precip_opt','[0:4]')
    call chk_var(mparwtr,'mparwtr','[1.0e-5:1.0e-2]')
    if (ccw_for_precip_opt == 4) then
      call chk_var(fac_qsat,'fac_qsat','[0.0:5.0]')
    end if
    call chk_var(ent_opt_dp,'ent_opt_dp','[0:9]')
    call chk_var(ent_fac_dp,'ent_fac_dp','[0.0:2.0]')
    if (ent_opt_dp == 3) then
      call chk_var(ent_dp_power,'ent_dp_power','[0.0:10.0]')
    end if

    call chk_var(ent_opt_md,'ent_opt_md','[0,1,2,3,6,7,8]')
    call chk_var(ent_fac_md,'ent_fac_md','[0.0:2.0]')
    if (ent_opt_md == 3) then
      call chk_var(ent_md_power,'ent_md_power','[0.0:10.0]')
    end if

    call chk_var(adapt,'adapt','[0:8]')
    call chk_var(r_det,'r_det','[0.0:1.0]')
    call chk_var(amdet_fac,'amdet_fac','[0.0:20.0]')

    ! Convective closure checking

    call chk_var(cldbase_opt_dp,'cldbase_opt_dp','[0,1,2,3,4,6,7,8,9]')
    call chk_var(cldbase_opt_md,'cldbase_opt_md','[0,1,2,3,4,6,7]')
    call chk_var(cldbase_opt_sh,'cldbase_opt_sh','[0,1]')
    call chk_var(cape_timescale,'cape_timescale','[1.0:9999999.0]')

    ! Variables only need checking for certain settings
    if (cldbase_opt_dp == 3 .or. cldbase_opt_dp == 4 .or.                      &
        cldbase_opt_dp == 6 .or. cldbase_opt_md == 3 .or.                      &
        cldbase_opt_md == 4 .or. cldbase_opt_md == 6) then
      call chk_var(cape_bottom,'cape_bottom','[1:10000]')
      call chk_var(cape_top,'cape_top','[2:10000]')
      ! Also check cape_top > cape_bottom or force error
      if (cape_bottom >= cape_top) then
        cmessage='cape_top must be > cape_bottom'
        icode = 1
        call ereport(RoutineName,icode,cmessage)
      end if
      call chk_var(w_cape_limit,'w_cape_limit','[0.0:10000.0]')
    end if
    if (cldbase_opt_dp == 4 .or. cldbase_opt_md == 4 ) then
      call chk_var(cape_min,'cape_min','[0.0:2000.0]')
    end if
    if (cldbase_opt_dp == 7 .or. cldbase_opt_dp == 8 .or.                      &
        cldbase_opt_md == 7) then
      call chk_var(cape_ts_min,'cape_ts_min','[0.0:14400.0]')
      call chk_var(cape_ts_max,'cape_ts_max','[0.0:14400.0]')
      if (cape_ts_min > cape_ts_max) then
        cmessage='cape_ts_min must be <= cape_ts_max'
        icode = 1
        call ereport(RoutineName,icode,cmessage)
      end if
    end if

    ! Downdraught and precipitation phase change options and variables
    call chk_var(dd_opt,'dd_opt','[0,1]')

    call chk_var(pr_melt_frz_opt,'pr_melt_frz_opt','[0,1,2]')

    if (pr_melt_frz_opt == 1 .or.  pr_melt_frz_opt == 2 ) then
      call chk_var(t_melt_snow,'t_melt_snow','[273.15:300.0]')
    end if

    ! CMT checking
    if (l_mom) then
      call chk_var(deep_cmt_opt,'deep_cmt_opt','[0:6]')
      call chk_var(mid_cmt_opt,'mid_cmt_opt','[0:2]')
    end if

    ! Cloud panel checks
    if (l_3d_cca) then

      if (l_anvil) then
        call chk_var(anv_opt,'anv_opt','[0,1,2]')
        call chk_var(anvil_factor,'anvil_factor','[0.0:3.0]')
        call chk_var(tower_factor,'tower_factor','[0.0:3.0]')
      end if

      if (l_ccrad) then
        call chk_var(cca2d_dp_opt,'cca2d_dp_opt','[0,1,2]')
        call chk_var(cca2d_md_opt,'cca2d_md_opt','[0,1,2]')
        call chk_var(cca2d_sh_opt,'cca2d_sh_opt','[0,1,2]')
        call chk_var(cca_dp_knob,'cca_dp_knob','[0.0:10.0]')
        call chk_var(cca_md_knob,'cca_md_knob','[0.0:10.0]')
        call chk_var(cca_sh_knob,'cca_sh_knob','[0.0:10.0]')
        call chk_var(ccw_dp_knob,'cca_dp_knob','[0.0:10.0]')
        call chk_var(ccw_md_knob,'cca_md_knob','[0.0:10.0]')
        call chk_var(ccw_sh_knob,'cca_sh_knob','[0.0:10.0]')
      end if
    end if

    call chk_var(rad_cloud_decay_opt,'rad_cloud_decay_opt','[0,1,2]')

    if (rad_cloud_decay_opt==1 .or. rad_cloud_decay_opt==2) then
      call chk_var(cld_life_opt,'cld_life_opt','[0,1]')
      call chk_var(cca_min,'cca_min','[0.001:0.4]')
      call chk_var(fixed_cld_life,'fixed_cld_life','[0.0:99999.0]')
    end if

    ! Updraught factor only used for these switches
    if (l_fix_udfactor .or. (ccw_for_precip_opt == 0) .or.                     &
       (ccw_for_precip_opt == 1) ) then
      call chk_var(ud_factor,'ud_factor','[0.0:1.0]')
    end if

  end if ! 5a or 6a scheme checking

  if ( i_convection_vn == i_cv_betts) then
    ! check the betts critical RH is not a percentage
    call chk_var(betts_rh,'betts_rh','[0.0:1.0]')
    !check the rate of precip evap as cloud is not a percentage (for betts)
    call chk_var(betts_cld_evap_const,'betts_cld_evap_const','[0.0:1.0]')
  end if

  ! LLCS checks
  if (i_convection_vn == i_cv_llcs) then
    ! Check that the cloud-precip switch is within options available
    call chk_var(llcs_cloud_precip, 'llcs_cloud_precip', '[0,1,2,3]')
    ! Check the critical RH for LLCS is not a percentage
    call chk_var(llcs_rhcrit, 'llcs_rhcrit', '[0.0:1.0]')
    ! Check the timescale
    call chk_var(llcs_timescale, 'llcs_timescale', '[1.0:99999.0]')
    if (llcs_cloud_precip /= llcs_opt_all_rain) then
      ! Check detrainment coefficient
      call chk_var(llcs_detrain_coef, 'llcs_detrain_coef', '[0.0:1.0]')
    end if
    if (llcs_cloud_precip == llcs_opt_crit_condens) then
      ! Check condensate profile limits
      call chk_var(qlmin, 'qlmin', '[0.0:1.0e-3]')
      call chk_var(mparwtr, 'mparwtr', '[1.0e-5:1.0e-2]')
      ! Check the factor controlling the threshold condensate profile
      call chk_var(fac_qsat, 'fac_qsat', '[0.0:5.0]')
    end if
    if (llcs_cloud_precip == llcs_opt_const_frac) then
      ! Check the rain fraction for LLCS is not a percentage
      call chk_var(llcs_rain_frac, 'llcs_rain_frac', '[0.0:1.0]')
    end if
  end if

  ! Variable checking for the CoMorph convection scheme
  if ( i_convection_vn == i_cv_comorph ) then

    ! CoMorph currently not set up to allow turning off wind increments
    if ( .not. l_mom ) then
      icode = 2
      write(cmessage,"(a)") "l_mom must be set to .true. if using "         // &
        "the CoMorph convection scheme, as it needs the wind "              // &
        "profiles to be interpolated onto the p-grid for input."
      call ereport(RoutineName,icode,cmessage)
    end if

    ! CoMorph needs the prognostic convective cloud fields to be
    ! enabled as 3-D fields, so throw an error if not l_ccrad and l_3d_cca
    if ( .not. ( l_ccrad .and. l_3d_cca ) ) then
      icode = 3
      write(cmessage,"(a)") "l_ccrad and l_3d_cca must be set to "          // &
        ".true. if using the CoMorph convection scheme, as it "             // &
        "always uses the 3-D prognostic convective cloud fraction and "     // &
        "convective cloud water content."
      call ereport(RoutineName,icode,cmessage)
    end if

  end if  ! ( i_convection_vn == i_cv_comorph )

else  ! Not planning to call a convection scheme so values should be
      ! trigger ignored which means logical must be false.
  ! Check no logicals set to true - at present prints warnings but should really
  ! be stopping the run.
  if (l_mom) then
    icode = -1
    write(cmessage,'(a29,a52)')' l_mom  SHOULD BE set .false.',                &
         ' as not calling convection - run will treat as false'
    call ereport(RoutineName,icode,cmessage)
  end if
  if (l_3d_cca) then
    icode = -2
    write(cmessage,'(a32,a43)')' l_3d_cca  SHOULD BE set .false.',             &
         ' as not calling convection so wasting space'
    call ereport(RoutineName,icode,cmessage)
  end if
  if (l_ccrad) then
    icode = -3
    write(cmessage,'(a31,a43)')' l_ccrad  SHOULD BE set .false.',              &
         ' as not calling convection so wasting space'
    call ereport(RoutineName,icode,cmessage)
  end if

end if

! Reset
def_src = ''

! Further checking
if (L_ccrad) then

  if (.not. l_3d_cca) then
    icode = 100
    CMessage    = '**error**: CCRad is not yet available without'//            &
                            ' the anvil scheme (L_3D_CCA = .true.)'

    call ereport(RoutineName, icode, CMessage)
  end if

  if (l_fix_udfactor) then
    icode = 101
    CMessage    = '**error**: L_CCRad and l_fix_udfactor'//                    &
                            ' should not be both set to true.'

    call ereport(RoutineName, icode, CMessage)
  end if

end if
!

! convective cold pools
!
if (cnv_cold_pools > ccp_off) then

  if (.not. l_param_conv) then
    icode = 103
    CMessage    = '**error**: Cnv_cold_pools requires'//                       &
                            ' l_param_conv true.'
    call ereport(RoutineName, icode, CMessage)
  end if
  !
  if (i_convection_vn /= i_convection_vn_6a) then
    icode = 104
    CMessage    = '**error**: Cnv_cold_pools requires'//                       &
                            ' the 6A convection scheme.'
    call ereport(RoutineName, icode, CMessage)
  end if
  !
  if (.not. l_ccrad) then
    icode = 106
    CMessage    = '**error**: Cnv_cold_pools requires'//                       &
                            ' l_ccrad true.'
    call ereport(RoutineName, icode, CMessage)
  end if

end if

!---------------------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!---------------------------------------------------------------------------
return
end subroutine check_run_convection


subroutine print_nlist_run_convection()
use umPrintMgr, only: umPrint
implicit none
character(len=50000) :: lineBuffer
real(kind=jprb) :: zhook_handle

character(len=*), parameter :: RoutineName='PRINT_NLIST_RUN_CONVECTION'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

call umPrint('Contents of namelist run_convection',                            &
    src='cv_run_mod')

write(lineBuffer,*)' i_convection_vn = ',i_convection_vn
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' l_mom = ',l_mom
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' l_fix_udfactor = ',l_fix_udfactor
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' l_snow_rain = ',l_snow_rain
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' l_eman_dd = ',l_eman_dd
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' l_cloud_deep = ',l_cloud_deep
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' l_rediagnosis = ',l_rediagnosis
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' l_dcpl_cld4pc2 = ',l_dcpl_cld4pc2
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' l_anvil = ',l_anvil
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' l_murk_conv = ',l_murk_conv
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' l_safe_conv = ',l_safe_conv
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' l_ccrad = ',l_ccrad
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' l_3d_cca = ',l_3d_cca
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' l_conv_hist = ',l_conv_hist
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' l_param_conv = ',l_param_conv
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' l_cv_conserve_check = ',l_cv_conserve_check
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' l_cmt_heating = ',l_cmt_heating
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A22,L1)')' l_conv_prog_dtheta = ',l_conv_prog_dtheta
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A22,L1)')' l_conv_prog_dq     = ',l_conv_prog_dq
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A22,L1)')' l_conv_prog_flx    = ',l_conv_prog_flx
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A22,L1)')' l_conv_prog_precip = ',l_conv_prog_precip
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A24,L1)')' adv_conv_prog_dtheta = ',adv_conv_prog_dtheta
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A24,L1)')' adv_conv_prog_dq     = ',adv_conv_prog_dq
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A24,L1)')' adv_conv_prog_flx    = ',adv_conv_prog_flx
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A,L1)')' l_prog_pert  = ',l_prog_pert
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A,L1)')' l_jules_flux = ',l_jules_flux
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A,L1)')' l_ccp_blv = ',l_ccp_blv
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A,L1)')' l_ccp_parcel_md = ',l_ccp_parcel_md
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A,L1)')' l_ccp_parcel_dp = ',l_ccp_parcel_dp
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A,L1)')' l_ccp_trig = ',l_ccp_trig
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A,L1)')' l_ccp_wind = ',l_ccp_wind
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A,L1)')' l_ccp_seabreeze = ',l_ccp_seabreeze
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A,L1)')' l_reset_neg_delthvu = ', l_reset_neg_delthvu
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A,L1)')' l_cvdiag_ctop_qmax = ', l_cvdiag_ctop_qmax
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A11,L1)')' l_fcape = ',l_fcape
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' n_conv_calls = ',n_conv_calls
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' sh_pert_opt = ',sh_pert_opt
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A,I0)')' md_pert_opt = ',md_pert_opt
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' dd_opt = ',dd_opt
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' deep_cmt_opt = ',deep_cmt_opt
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' mid_cmt_opt = ',mid_cmt_opt
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' termconv = ',termconv
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' adapt = ',adapt
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' fdet_opt = ',fdet_opt
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' r_det = ',r_det
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' tice = ',tice
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' qstice = ',qstice
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A,F0.4)')' ent_fac_sh = ',ent_fac_sh
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' ent_fac_dp = ',ent_fac_dp
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' ent_fac_md = ',ent_fac_md
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' ent_opt_dp = ',ent_opt_dp
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' ent_opt_md = ',ent_opt_md
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' ent_dp_power = ',ent_dp_power
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' ent_md_power = ',ent_md_power
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A,F0.3)')' efrac = ',efrac
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A,F0.3)')' max_pert_scale = ',max_pert_scale
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A,F0.3)')' thpixs_mid = ',thpixs_mid
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' bl_cnv_mix = ',bl_cnv_mix
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' mid_cnv_pmin = ',mid_cnv_pmin
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A,F0.4)')' orig_mdet_fac = ',orig_mdet_fac
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' amdet_fac = ',amdet_fac
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' ccw_for_precip_opt = ',ccw_for_precip_opt
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' cnv_wat_load_opt = ',cnv_wat_load_opt
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' tv1_sd_opt = ',tv1_sd_opt
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A15,I0)')' midtrig_opt = ',midtrig_opt
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' limit_pert_opt = ',limit_pert_opt
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' icvdiag = ',icvdiag
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' plume_water_load = ',plume_water_load
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' dil_plume_water_load = ',dil_plume_water_load
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' cvdiag_inv = ',cvdiag_inv
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' cvdiag_sh_wtest = ',cvdiag_sh_wtest
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' cldbase_opt_dp = ',cldbase_opt_dp
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' cldbase_opt_md = ',cldbase_opt_md
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' cldbase_opt_sh = ',cldbase_opt_sh
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' cape_bottom = ',cape_bottom
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' cape_top = ',cape_top
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' cape_timescale = ',cape_timescale
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' w_cape_limit = ',w_cape_limit
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' cape_min = ',cape_min
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A,F0.3)')' cape_ts_min = ',cape_ts_min
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A,F0.3)')' cape_ts_max = ',cape_ts_max
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A,F0.4)')' c_mass_sh = ',c_mass_sh
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' cld_life_opt = ',cld_life_opt
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' rad_cloud_decay_opt = ',rad_cloud_decay_opt
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' cca_min = ',cca_min
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' fixed_cld_life = ',fixed_cld_life
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' ud_factor = ',ud_factor
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' mparwtr = ',mparwtr
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' qlmin = ',qlmin
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' fac_qsat = ',fac_qsat
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' eff_dcfl = ',eff_dcfl
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' eff_dcff = ',eff_dcff
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(a15,e13.6)')' cpress_term = ',cpress_term
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' cca2d_sh_opt = ',cca2d_sh_opt
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' cca_sh_knob = ',cca_sh_knob
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' ccw_sh_knob = ',ccw_sh_knob
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' cca2d_md_opt = ',cca2d_md_opt
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' cca_md_knob = ',cca_md_knob
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' ccw_md_knob = ',ccw_md_knob
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' cca2d_dp_opt = ',cca2d_dp_opt
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' cca_dp_knob = ',cca_dp_knob
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' ccw_dp_knob = ',ccw_dp_knob
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' anvil_factor = ',anvil_factor
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' tower_factor = ',tower_factor
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' anv_opt = ',anv_opt
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' t_melt_snow = ',t_melt_snow
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' iconv_shallow = ',iconv_shallow
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' iconv_mid = ',iconv_mid
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' iconv_deep = ',iconv_deep
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' iconv_congestus = ',iconv_congestus
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,*)' tau_conv_prog_precip = ',tau_conv_prog_precip
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A24,F15.3)')' tau_conv_prog_dtheta = ',tau_conv_prog_dtheta
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A24,F15.3)')' tau_conv_prog_dq     = ',tau_conv_prog_dq
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A24,F15.3)')' tau_conv_prog_flx    = ',tau_conv_prog_flx
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(a17,f15.3)')' prog_ent_grad = ',prog_ent_grad
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(a17,f15.3)')' prog_ent_int =  ',prog_ent_int
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(a17,f15.3)')' prog_ent_max =  ',prog_ent_max
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(a17,f15.3)')' prog_ent_min =  ',prog_ent_min
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A19,I0)')' pr_melt_frz_opt = ',pr_melt_frz_opt
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A,I0)')' cnv_cold_pools = ',cnv_cold_pools
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A,E10.3)')' kappa_g = ',kappa_g
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A,E10.3)')' kappa_h = ',kappa_h
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A,E10.3)')' remain_max = ',remain_max
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A,E10.3)')' phi_ccp = ',phi_ccp
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A,E10.3)')' ccp_buoyancy = ',ccp_buoyancy
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A,E10.3)')' ccp_nonelastic = ',ccp_nonelastic
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A,E10.3)')' entrain_max = ',entrain_max
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A,E10.3)')' entrain_min = ',entrain_min
call umPrint(lineBuffer,src='cv_run_mod')

write(lineBuffer,'(A,F13.1)')' betts_tau =  ',betts_tau
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A,F17.5)')' betts_rh =  ',betts_rh
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A,F17.5)')' betts_cld_evap_const =  ',betts_cld_evap_const
call umPrint(lineBuffer,src='cv_run_mod')

write(lineBuffer,'(A,I0)')' llcs_cloud_precip = ',llcs_cloud_precip
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A,F17.5)')' llcs_detrain_coef = ',llcs_detrain_coef
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A,F17.5)')' llcs_rhcrit = ',llcs_rhcrit
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A,F13.1)')' llcs_timescale = ',llcs_timescale
call umPrint(lineBuffer,src='cv_run_mod')
write(lineBuffer,'(A,F13.1)')' llcs_rain_frac = ',llcs_rain_frac
call umPrint(lineBuffer,src='cv_run_mod')

call umPrint('- - - - - - end of namelist - - - - - -',                        &
    src='cv_run_mod')

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine print_nlist_run_convection


end module cv_run_mod
