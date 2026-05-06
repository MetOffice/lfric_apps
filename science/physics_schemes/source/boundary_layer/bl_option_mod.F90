! *****************************COPYRIGHT********************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT********************************

!  Data module for switches/options concerned with the BL scheme.
! Description:
!   Module containing runtime options/data used by the boundary
!   layer scheme

! Method:
!   Switches and associated data values used by the boundary layer
!   scheme are defined here and assigned default values. These may
!   be overridden by namelist input.

!   Any routine wishing to use these options may do so with the 'use'
!   statement.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer

! Code Description:
!   Language: FORTRAN 95

module bl_option_mod

use missing_data_mod, only: rmdi, imdi
use visbty_constants_mod, only: calc_prob_of_vis
use mym_option_mod, only: bdy_tke, tke_levels, l_local_above_tkelvs,           &
      l_my_initialize, l_my_ini_zero, my_ini_dbdz_min,                         &
      l_adv_turb_field, l_my_condense, l_shcu_buoy, shcu_levels,               &
      wb_ng_max, my_lowest_pd_surf, l_my_prod_adj,                             &
      my_z_limit_elb, l_print_max_tke, tke_cm_mx, tke_cm_fa, tke_dlen,         &
      l_3dtke
use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim
use errormessagelength_mod, only: errormessagelength
use chk_opts_mod, only: chk_var, def_src

use um_types, only: real_umphys, r_bl

implicit none

!=======================================================================
! BL namelist options - ordered as in meta-data sort-key, with namelist
! option and possible values
!=======================================================================

! 01 Switch to determine version of boundary layer scheme
integer :: i_bl_vn = imdi
! 0 => no boundary layer scheme
integer, parameter :: i_bl_vn_0  = 0
! 1 => 1A prognostic TKE based turbulent closure
integer, parameter :: i_bl_vn_1a = 1
! 3 => 9C revised entr fluxes plus new scalar flux-grad
integer, parameter :: i_bl_vn_9c = 3

! 02 Stable boundary layer stability function option
integer :: sbl_op= imdi
! Options for stable boundary layers
integer, parameter ::  long_tails           = 0
integer, parameter ::  sharpest             = 1
integer, parameter ::  sharp_sea_long_land  = 2
integer, parameter ::  mes_tails            = 3
integer, parameter ::  louis_tails          = 4
integer, parameter ::  depth_based          = 5 ! not available in LFRic
integer, parameter ::  sharp_sea_mes_land   = 6
integer, parameter ::  lem_stability        = 7
integer, parameter ::  sharp_sea_louis_land = 8
integer, parameter ::  equilibrium_sbl      = 9

! 03 Sets transition Ri for Sharpest function
real(kind=r_bl) :: ritrans = rmdi  ! 1/g0 so 0.1 for original Sharpest
                        ! (closer to 0.2 makes it sharper)

! 04 Switch for what interpolations to do to calculate the local Kh
integer :: i_interp_local = imdi
! Allowed values:
! 1) Interpolate the vertical gradients of sl,qw and calculate
!    stability dbdz and Kh on theta-levels
integer, parameter :: i_interp_local_gradients = 1
! 2) Interpolate cloud-fraction and calculate stability dbdz on
!    rho-levels, but still calculate Kh on theta-levels
integer, parameter :: i_interp_local_cf_dbdz = 2

! 05 Options for convective BL stability functions
integer :: cbl_op = imdi
! UM_std   (=0) => Standard UM
integer, parameter ::  um_std     = 0
! neut_cbl(=1) => Keep neutral stability for Ri<0
integer, parameter ::  neut_cbl   = 1
! LEM_conv (=2) => Conventional LEM
integer, parameter ::  lem_conven = 2
! LEM_std (=3) => Standard LEM
integer, parameter ::  lem_std    = 3
! LEM_adjust (=4) => Adjustable LEM (adjustable between 'Conventional' and
! 'Standard' models)
integer, parameter ::  lem_adjust  = 4
! 05a Parameter to control mixing in convective boundary layers between
! 'Conventional' and 'Standard' models
real(kind=r_bl) :: cbl_mix_fac_nml = rmdi

! 06 Minimum value of the mixing length (m)
real(kind=r_bl) :: lambda_min_nml = rmdi

! 07 Switch to mix qcf in the vertical
logical :: l_bl_mix_qcf = .false.

! 08 Switch to reset q and qcf to zero if negative after imp_solver in
!    ni_imp_ctl (non-conservative)
logical :: l_reset_neg_q=.false.

! 09 Switch for free atmospheric mixing options
integer :: local_fa= imdi
! to_sharp_across_1km (=1) => smoothly switch to sharpest across 1km AGL
integer, parameter :: to_sharp_across_1km = 1
! ntml_level_corrn (=2) => correct level for ntml_local (downwards)
integer, parameter :: ntml_level_corrn    = 2
! free_trop_layers (=3) => as "ntml_level_corrn" but also diagnose
! FA turbulent layer depths
integer, parameter :: free_trop_layers    = 3

! 10 Switch to keep local mixing in free atmosphere
integer :: Keep_Ri_FA = imdi
! Only reduce local K to zero at discontinuous inversions
integer, parameter :: except_disc_inv = 2

! 11 SBL mixing dependent on sg orography
integer :: sg_orog_mixing = imdi
!   1 - extending length of SHARPEST tail following McCabe and Brown (2007)
integer, parameter :: extended_tail = 1
!   2 - include subgrid drainage shear in Ri and diffusion coefficients,
integer, parameter :: sg_shear = 2
!   3 - as 2 + orographically enhanced lambda_m (note lambda_h enhancement
integer, parameter :: sg_shear_enh_lambda = 3
! not included in bdy_expl2 in error but now operational in UKV) and smooth
! decrease to lambda_min above

! 12 Switch to apply heating source from turbulence dissipation
integer :: fric_heating = imdi
! OFF (=0) => not used
! ON  (=1) => used

! 13 Options for surface fluxes solver (explicit vs implicit)
integer :: flux_bc_opt = imdi
! Options for surface fluxes solver (explicit vs implicit)
integer, parameter :: interactive_fluxes  = 0
! use land surface model to calculate fluxes interactively with PBL
integer, parameter :: specified_fluxes_only  = 1
! only specify sensible and latent surface fluxes, flux_h and flux_e
integer, parameter :: specified_fluxes_tstar = 2
! as (1) plus specify surface T, tstar (aka t_surf in atmos_phys2)
integer, parameter :: specified_fluxes_cd    = 3
! as (2) plus use specified ustar to give CD (still implicit solver)

! 14 parameter in visibility diagnostic (from visbty_constants_mod)
! calc_prob_of_vis

! 15 Tunable parameter used in the calculation of the wind gust
real(kind=real_umphys) :: c_gust=rmdi  ! used to be 4.0
                                       ! EC have 7.2/2.29=3.14

!=======================================================================
! Options specific to 9C (non-local) schemes
!=======================================================================

! 01 Height of number of levels for non_local scheme
real(kind=r_bl) :: z_nl_bl_levels= rmdi

! 02 Switch to use zi/L in the diagnosis of shear-driven boundary layers:
integer :: idyndiag = imdi
! OFF (=0) => not used
! DynDiag_ZL (=1) => over sea uses zi/L but gives problems when
! BL_LEVELS is high
integer, parameter :: DynDiag_ZL = 1
! DynDiag_ZL_corrn (=2) => as 1 but copes with high BL_LEVELS
integer, parameter :: DynDiag_ZL_corrn = 2
! DynDiag_ZL_CuOnly (=3) => as 2 but only applied to points diagnosed
! with Cumulus and strictly for sea points (fland<0.01, cf 0.5)
integer, parameter :: DynDiag_ZL_CuOnly = 3
! DynDiag_Ribased (=4) => as 3 but also overrides Cumulus diagnosis if
!     ZH(Ri) > ZLCL+zhloc_depth_fac*(ZHPAR-ZLCL)
! Note that here Ri accounts for gradient adjustment by the non-local scheme.
integer, parameter :: DynDiag_Ribased = 4

! 03 For idyndiag options DynDiag_Ribased: the fractional height into
! the cloud layer reached by the local Ri-based BL depth calculation
real(kind=r_bl) :: zhloc_depth_fac = rmdi ! suggested 0.5

! 03a For all non-zero idyndiag options: the value of z/l below which
!     to diagnose shear driven, as opposed to cumulus, layers
real(kind=r_bl) :: near_neut_z_on_l = rmdi ! Suggest 1.6

! 04
logical :: l_use_surf_in_ri = .false.
                       ! dbdz in Ri calculation on theta level 1 (ri(k=2))
                       ! false => use dbdz between levels 1 and 2
                       ! true  => use Tstar to give centred difference

! 05 Switch for revised flux-gradient relationships (9C Only)
integer :: flux_grad= imdi
! Flux gradients as in Lock et al. (2000)
integer, parameter :: Locketal2000   = 0
! Flux gradients as in Lock et al (2000) but using coefficients from
! Holtslag and Boville (1993)
integer, parameter :: HoltBov1993 = 1
! Flux gradients as in Lock and Whelan (2006)
integer, parameter :: LockWhelan2006 = 2

! 06 Method of doing entrainment at dsc top (9C only)
integer :: entr_smooth_dec = imdi
! 0 - old method, include surface terms only in weakly coupled,
! no cumulus situations
! 1 - taper method - no hard limit but reduce surface terms according
! to svl difference
! 2 - as 1 but also taper height-scale used in w* between zh and zhsc,
! instead of jumping straight to zhsc when coupled.
integer, parameter :: entr_taper_zh = 2

! 19 Switch to use miscellaneous minor fixes and improvements in the
! treatment of the surface mixed-layer and decoupled stratocumulus and their
! coupling / decoupling.
logical :: l_use_sml_dsc_fixes = .false.

! 07 Switch to reset the decoupling threshold to dec_thres_cloud for
! the Ktop iteration (in case it has been set to cu or clear values)
logical :: l_reset_dec_thres = .false.

! 07a Decoupling threshold for cloudy boundary layers
!    (larger makes decoupling less likely) (9C Only)
real(kind=r_bl) :: dec_thres_cloud = rmdi ! Suggest 0.1

! 08 Options for diagnosing stratocumulus layers
integer :: sc_diag_opt = imdi
! 0) Original code: diagnose Sc based on svl gradient at non-cumulus points
!    but set to the conv_diag parcel top at shallow cumulus points.
integer, parameter :: sc_diag_orig = 0
! 1) As (0) but also set the Sc-top to the conv_diag parcel top at
!    "deep" cumulus points.
integer, parameter :: sc_diag_cu_relax = 1
! 2) At all cumulus points, diagnose Sc-top based on max total-water RH
!    instead of using the conv_diag parcel top.
integer, parameter :: sc_diag_cu_rh_max = 2
! 3) Diagnose Sc-top based on max total-water RH at all points
!    instead of using the svl gradient.
integer, parameter :: sc_diag_all_rh_max = 3

! 09 Options for non-local mixing across the LCL in Cu (9C only)
integer :: kprof_cu  = imdi
! klcl_entr (=1) => set K-profile at the LCL using the standard CBL
! entrainment parametrization
integer, parameter :: klcl_entr = 1
! buoy_integ (=2) => use buoyancy flux integration, as used in decoupling
integer, parameter :: buoy_integ = 2
! buoy_integ_low (=3) => use buoyancy flux integration and allow top of mixing
!                        to be below the lcl
integer, parameter :: buoy_integ_low = 3

! 10 Switch for parametrizing fluxes across a resolved inversion (9C only)
integer :: bl_res_inv= imdi      ! OFF (=0) => not used
integer, parameter :: cosine_inv_flux = 1
! cosine_inv_flux (=1) => specify the inversion flux profile shape as a cosine
integer, parameter :: target_inv_profile = 2
! target_inv_profile (=2) => specify the inversion flux profile to bring the
!                            svl profile to a target piecewise linear shape

! 11 Options for blending with 3D Smagorinsky
integer :: blending_option = imdi
! blend_allpoints (=1) => blending everywhere
integer, parameter :: blend_allpoints = 1
! blend_except_cu (=2) => blending everywhere except in cumulus layers
integer, parameter :: blend_except_cu = 2
! blend_gridindep_fa (=3) => as blending everywhere but relax to 1D BL
! above the BL and in a way that is independent of grid size
integer, parameter :: blend_gridindep_fa =3
! blend_cth_shcu_only (=4) => as blend_gridindep_fa but inc Smag mixing length
! as minimum, and use cloud-top as blending length-scale only for shallow cu
integer, parameter :: blend_cth_shcu_only =4

! 11a max permitted cloud top height in metres for shallow cu, used for
!     blending_option=blend_cth_shcu_only
real(kind=r_bl) :: shallow_cu_maxtop = rmdi

! 12 empirical parameter scaling ustar entrainment term (9C only)
real(kind=r_bl) :: a_ent_shr_nml = rmdi ! suggested 1.6

! 12a empirical parameter scaling evaporative entrainment term
real(kind=r_bl) :: a_ent_2 = rmdi ! suggested 0.056

! 13 Improved method to calculate cloud-top radiative flux jump
logical :: l_new_kcloudtop = .false.

! 15 Multiplicative tuning factor in TKE diagnosis, usually 1.0
real(kind=r_bl) :: tke_diag_fac = rmdi

! 16 Include convective effects in TKE diagnostic
logical :: l_conv_tke = .false.

! 17 Use separate rhokm_ent arrays to hold momentum entrainment coefficients
!    (preferred)
logical :: l_use_var_fixes = .false.
! Use fixes (separate rhokm_ent variables instead of rhokm array,
! don't create DSC layers when no top-driven turbulence and
! set zh to max(zh, zsml_top), not just ntml to keep consistent

! 18 Switch to ignore cloud ice (qcf) in the BL scheme
logical :: l_noice_in_turb = .false.

! 19 Options for blending height level
! Implicitly couple at model level one
integer, parameter :: blend_level1  = 0
! Implicitly couple to a specified (fixed) blending height
integer, parameter :: blend_levelk  = 1
! Implicitly couple at a prognostically determined blending height
integer, parameter :: blend_levelvar = 2
! Switch to determine blending height for implicit coupling
integer :: blend_height_opt = imdi

! 19a Blending height (m) used to calculate the level for surface coupling
real(kind=r_bl) :: h_blend_fix = rmdi

! Options for non-gradient stress scheme following Brown and Grant (1997)
! Not a parameter as needs to be switched off for 3D Smag
! (divide by zhnl=0 in ex_flux_uv)
integer :: ng_stress = imdi
! Brown and Grant (1997), modified version (applied between 0.1*ZH and ZH,
! where ZH is the ABL depth) with no limit placed on its size
integer, parameter :: BrownGrant97 = 1
! Brown and Grant (1997), modified version, including a limit on its size
integer, parameter :: BrownGrant97_limited = 2
! Brown and Grant (1997) as originally defined, applied between the
! ground and ZH, including a limit on its size
integer, parameter :: BrownGrant97_original = 3

! Switch for vertical discretization of the Sc-top radiatively-cooled layer,
! in the event that the default depth dzrad is too shallow compared to
! the model vertical resolution.
integer :: dzrad_disc_opt = imdi
! Set base of layer at the theta-level below ntml (1.5-2.5 levels below zh)
integer, parameter :: dzrad_ntm1 = 1
! Set base of layer 1.5 levels below zh (smoothly varies without jumping)
integer, parameter :: dzrad_1p5dz = 2

! Number of iterations used to converge buoyancy-flux criteria
integer :: num_sweeps_bflux = imdi

! Improve accuracy of mixed-layer depths found via buoyancy-flux integration,
! by accounting for sensitivity of gradient adjustment to the depth.
logical :: l_converge_ga = .false.

!=======================================================================
! Implicit solver options
!=======================================================================
! 01 Options for where in the timestep to do the boundary-layer implicit
!     solve and add on its increments
integer :: i_impsolve_loc = imdi
integer, parameter :: i_impsolve_befconv=1  ! Call ni_imp_ctl before convection
integer, parameter :: i_impsolve_aftconv=2  ! Call ni_imp_ctl after convection

! 02 Time weights for boundary layer levels
real(kind=real_umphys) :: alpha_cd_in (2) = rmdi
real(kind=real_umphys), allocatable :: alpha_cd(:)

! Parameters for uncond stable numerical solver
! 03 Puns : used in an unstable BL column
real(kind=real_umphys) :: puns = rmdi ! suggested 0.5
! 04 Pstb : used in an stable BL column
real(kind=real_umphys) :: pstb = rmdi ! suggested 2.0

!=======================================================================
! a few idealised switches available with the 0A scheme
!=======================================================================
! Choice of Rayleigh friction
integer :: fric_number = imdi
! Options for Rayleigh friction
integer, parameter :: fric_none = 0
integer, parameter :: fric_HeldSuarez_1 = 1
integer, parameter :: fric_HeldSuarez_2 = 2

!=======================================================================
! Options not in namelist, Integers
!=======================================================================

! Generic switch options
integer, parameter :: off = 0  ! Switch disabled
integer, parameter :: on  = 1  ! Switch enabled

! Number of levels for non_local sheme
integer :: nl_bl_levels= imdi

! Options for Prandtl number (in local Ri scheme)
integer, parameter ::  LockMailhot2004 = 1
! Switch for Prandtl number options
integer, parameter :: Prandtl= LockMailhot2004

! Switch for shear-dominated b.l.
! Not a parameter as needs switching off with 1A scheme
integer :: ishear_bl = on

! Switch on the non-local scheme
! Not a parameter as needs switching off for 3D Smag
integer :: non_local_bl = on

! Switch to use implicit weights of 1 for tracers (if on)
integer, parameter :: trweights1 = on

! Switch to enhance entrainment in decoupled stratocu over cu (9C Only)
! OFF (=0) => not used
! Buoyrev_feedback (=1) => buoyancy reversal feedback enhanced by Cu
! for 0<D<0.1
integer, parameter :: Buoyrev_feedback = 1
integer, parameter :: entr_enhance_by_cu = Buoyrev_feedback

! 03 Switch to allow different critical Richardson numbers in the diagnosis
! of BL depth
! 0 => RiC=1 everywhere
! 1 => RiC=0.25 for SHARPEST, 1 otherwise
integer, parameter :: Variable_RiC = on

!=======================================================================
! Options not in namelist, Reals
!=======================================================================
! Weighting of the Louis tail towards long tails:
real(kind=r_bl), parameter :: WeightLouisToLong = 0.0_r_bl
! 0.0 = Louis tails
! 1.0 = Long tails

! Maximum value of the mixing length (m)
real(kind=r_bl), parameter :: lambda_max_nml = 500.0_r_bl

! fraction of BL depth for mixing length
real(kind=r_bl), parameter :: lambda_fac = 0.15_r_bl

! Maximum implied stress gradient across the boundary layer, used to limit
! the explicit stress applied in non-local scalings (m/s2)
real(kind=r_bl), parameter :: max_stress_grad = 0.05_r_bl

! Max height above LCL for K profile in cumulus
real(kind=r_bl), parameter :: max_cu_depth = 250.0_r_bl

! timescale for drainage flow development (s)
real(kind=r_bl), parameter :: t_drain = 1800.0_r_bl

! Horizontal scale for drainage flows (m)
! Currently hard-wired for UKV run over UK
real(kind=r_bl), parameter :: h_scale = 1500.0_r_bl

! Maximum allowed value of the Prandtl number with the LockMailhot2004 option
real(kind=r_bl), parameter :: pr_max  = 5.0_r_bl

! Parameters governing the speed of transition from 1D BL to Smagorinsky.
! beta = 0.15 matches Honnert et al well in the BL but a faster transition
! seems appropriate in the free atmosphere above.
real(kind=r_bl), parameter :: beta_bl = 0.15_r_bl
real(kind=r_bl), parameter :: beta_fa = 1.0_r_bl

! Constants in linear part of the blending weight function
! linear0 and linear1 are the values of delta/z_scale where linear factor
! equals 0 and 1 respectively
real(kind=r_bl), parameter :: linear1=0.25_r_bl
real(kind=r_bl), parameter :: linear0=4.0_r_bl
real(kind=r_bl), parameter :: rlinfac=1.0_r_bl/(linear0-linear1)

! Critical Ri when SHARPEST stability functions used
real(kind=r_bl), parameter :: RiCrit_sharp = 0.25_r_bl

! Parameter used in gradient adjustment
real(kind=r_bl), parameter :: a_grad_adj = 3.26_r_bl

! Parameter used in gradient adjustment
real(kind=r_bl), parameter :: max_t_grad = 1.0e-3_r_bl

! Powers to use in prescription of equilibrium profiles of stress and
! buoyancy flux in Equilib. SBL model
real(kind=r_bl), parameter :: Muw_SBL= 1.5_r_bl
real(kind=r_bl), parameter :: Mwt_SBL= 1.0_r_bl

! Pre-calculated powers for speeding up calculations
real(kind=r_bl), parameter :: one_third = 1.0_r_bl/3.0_r_bl
real(kind=r_bl), parameter :: two_thirds = 2.0_r_bl * one_third
real(kind=r_bl), parameter :: zero = 0.0_r_bl
real(kind=r_bl), parameter :: one = 1.0_r_bl
real(kind=r_bl), parameter :: one_half = 1.0_r_bl/2.0_r_bl

! Maximum value of TKE to me used in TKE diagnostic
real(kind=r_bl), parameter :: max_tke = 5.0_r_bl

! Decoupling threshold for cumulus sub-cloud layers
! (larger makes decoupling less likely) (9C Only)
real(kind=r_bl) :: dec_thres_cu = rmdi

! cloud fraction threshold required for a cloud layer to be diagnosed
real(kind=r_bl), parameter :: sc_cftol=0.1_r_bl
!=======================================================================
! Options not in namelist, logical
!=======================================================================
! Switch for coupled gradient method in Equilibrium SBL model
logical, parameter :: L_SBLco = .true.

! logical for whether to skip calculations based on start of timestep
! quantities when using semi-lagrangian cycling with Endgame
! This is set to true in dynamics_input_mod as Endgame always uses it,
! except for the 1a BL scheme where it is false
! N.B. results should bit-compare whether this logical is true or false
! so changing it there is a good test of whether new code has been
! added correctly
logical :: l_quick_ap2 = .false.

! Logical for whether to compute the explicit momentum fluxes on the
! p-grid on all model-levels.  These might be wanted for a new convection
! scheme, or just for diagnostic purposes.
logical :: l_calc_tau_at_p = .false.

!=======================================================================
!run_bl namelist
!=======================================================================

namelist/run_bl/ i_bl_vn, sbl_op, cbl_op, cbl_mix_fac_nml,                     &
    l_use_surf_in_ri, lambda_min_nml, ritrans, c_gust,                         &
    dzrad_disc_opt, num_sweeps_bflux, l_converge_ga,                           &
    local_fa, Keep_Ri_FA, l_bl_mix_qcf, l_conv_tke, l_use_var_fixes,           &
    l_reset_neg_q, tke_diag_fac, i_interp_local,                               &
    sg_orog_mixing, fric_heating, calc_prob_of_vis, z_nl_bl_levels,            &
    idyndiag, zhloc_depth_fac, flux_grad, entr_smooth_dec,                     &
    sc_diag_opt, kprof_cu, bl_res_inv, blending_option,                        &
    a_ent_shr_nml, a_ent_2, alpha_cd_in, puns, pstb,                           &
    l_use_sml_dsc_fixes, l_reset_dec_thres,                                    &
    bdy_tke, tke_levels, l_local_above_tkelvs, l_my_initialize,                &
    l_my_ini_zero, my_ini_dbdz_min, l_adv_turb_field, l_my_condense,           &
    l_shcu_buoy, shcu_levels, wb_ng_max, my_lowest_pd_surf,                    &
    l_my_prod_adj, my_z_limit_elb, l_print_max_tke, tke_cm_mx,                 &
    tke_cm_fa, tke_dlen, l_new_kcloudtop, fric_number, flux_bc_opt,            &
    i_impsolve_loc,                                                            &
    l_3dtke, l_noice_in_turb, shallow_cu_maxtop, dec_thres_cloud,              &
    near_neut_z_on_l, blend_height_opt, h_blend_fix, ng_stress

!DrHook-related parameters
integer(kind=jpim), parameter, private :: zhook_in  = 0
integer(kind=jpim), parameter, private :: zhook_out = 1

character(len=*), parameter, private :: ModuleName='BL_OPTION_MOD'

contains
subroutine print_nlist_run_bl()
use umPrintMgr, only: umPrint
implicit none
character(len=50000) :: lineBuffer
real(kind=jprb)      :: zhook_handle

character(len=*), parameter :: RoutineName='PRINT_NLIST_RUN_BL'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

call umPrint('Contents of namelist run_bl',                                    &
    src='bl_option_mod')

write(lineBuffer,'(A,I0)') 'i_bl_vn = ',i_bl_vn
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,I0)') 'sbl_op = ',sbl_op
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,I4)') 'cbl_op = ',cbl_op
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,ES12.4)') 'cbl_mix_fac_nml = ',cbl_mix_fac_nml
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,L1)') 'l_use_surf_in_ri = ',l_use_surf_in_ri
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,F0.5)') 'ritrans = ',ritrans
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,F0.5)') 'lambda_min_nml = ',lambda_min_nml
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,L1)') 'l_bl_mix_qcf = ',l_bl_mix_qcf
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,ES12.4)') 'c_gust = ',c_gust
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,L1)') 'l_reset_neg_q = ',l_reset_neg_q
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,I0)') 'local_fa = ',local_fa
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,I0)') 'Keep_Ri_FA = ',Keep_Ri_FA
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,I0)') 'sg_orog_mixing = ',sg_orog_mixing
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,I0)') 'fric_heating = ',fric_heating
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,F0.5)') 'calc_prob_of_vis = ',calc_prob_of_vis
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,ES12.4)') 'z_nl_bl_levels = ',z_nl_bl_levels
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,I0)') 'idyndiag = ',idyndiag
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,F0.5)') 'zhloc_depth_fac = ',zhloc_depth_fac
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,I0)') 'flux_grad = ',flux_grad
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,I0)') 'entr_smooth_dec = ',entr_smooth_dec
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,L1)') 'l_use_sml_dsc_fixes = ',l_use_sml_dsc_fixes
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,L1)') 'l_reset_dec_thres = ',l_reset_dec_thres
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,L1)') 'l_new_kcloudtop = ',l_new_kcloudtop
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,I0)') 'sc_diag_opt = ',sc_diag_opt
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,I0)') 'kprof_cu = ',kprof_cu
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,I0)') 'bl_res_inv = ',bl_res_inv
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,F0.5)') 'a_ent_shr_nml = ',a_ent_shr_nml
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,F0.5)') 'a_ent_2 = ',a_ent_2
call umPrint(lineBuffer,src='bl_option_mod')
write(linebuffer,'(A,I0)') 'blending_option = ',blending_option
call umprint(linebuffer,src='bl_option_mod')
write(linebuffer,'(A,L1)') 'l_conv_tke = ',l_conv_tke
call umprint(linebuffer,src='bl_option_mod')
write(linebuffer,'(A,I1)') 'dzrad_disc_opt = ',dzrad_disc_opt
call umprint(linebuffer,src='bl_option_mod')
write(linebuffer,'(A,I4)') 'num_sweeps_bflux = ',num_sweeps_bflux
call umprint(linebuffer,src='bl_option_mod')
write(linebuffer,'(A,L1)') 'l_converge_ga = ',l_converge_ga
call umprint(linebuffer,src='bl_option_mod')
write(linebuffer,'(A,L1)') 'l_use_var_fixes = ',l_use_var_fixes
call umprint(linebuffer,src='bl_option_mod')
write(linebuffer,'(A,ES12.4)') 'tke_diag_fac = ',tke_diag_fac
call umprint(linebuffer,src='bl_option_mod')
write(linebuffer,'(A,I4)') 'i_interp_local = ',i_interp_local
call umprint(linebuffer,src='bl_option_mod')
write(linebuffer,'(A,L1)') 'l_noice_in_turb = ',l_noice_in_turb
call umprint(linebuffer,src='bl_option_mod')
write(lineBuffer,'(A,ES12.4)') 'shallow_cu_maxtop = ',shallow_cu_maxtop
call umprint(linebuffer,src='bl_option_mod')
write(lineBuffer,'(2(A,F0.5))') 'alpha_cd_in = ',                              &
                                    alpha_cd_in(1),' ',alpha_cd_in(2)
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,F0.5)') 'puns = ',puns
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,F0.5)') 'pstb = ',pstb
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,I4)') 'flux_bc_opt = ',flux_bc_opt
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,I4)') 'i_impsolve_loc = ',i_impsolve_loc
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,I0)') 'bdy_tke = ',bdy_tke
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,I0)') 'tke_levels = ',tke_levels
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,L1)') 'l_local_above_tkelvs = ',l_local_above_tkelvs
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,L1)') 'l_my_initialize = ',l_my_initialize
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,L1)') 'l_my_ini_zero = ',l_my_ini_zero
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,F0.5)') 'my_ini_dbdz_min = ',my_ini_dbdz_min
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,L1)') 'l_adv_turb_field = ',l_adv_turb_field
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,L1)') 'l_my_condense = ',l_my_condense
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,L1)') 'l_shcu_buoy = ',l_shcu_buoy
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,I0)') 'shcu_levels = ',shcu_levels
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,F0.5)') 'wb_ng_max = ',wb_ng_max
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,I0)') 'my_lowest_pd_surf = ',my_lowest_pd_surf
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,L1)') 'l_my_prod_adj = ',l_my_prod_adj
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,F0.5)') 'my_z_limit_elb = ',my_z_limit_elb
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,L1)') 'l_print_max_tke = ',l_print_max_tke
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,F0.5)') 'tke_cm_mx = ',tke_cm_mx
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,F0.5)') 'tke_cm_fa = ',tke_cm_fa
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,I0)') 'tke_dlen = ',tke_dlen
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,I0)')'fric_number = ',fric_number
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,L1)')'l_3dtke = ',l_3dtke
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,F0.5)')'dec_thres_cloud = ',dec_thres_cloud
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,F0.5)')'near_neut_z_on_l = ',near_neut_z_on_l
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,I3)')'blend_height_opt = ',blend_height_opt
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,ES12.4)') 'h_blend_fix = ',h_blend_fix
call umPrint(lineBuffer,src='bl_option_mod')
write(lineBuffer,'(A,I0)') 'ng_stress = ',ng_stress
call umPrint(lineBuffer,src='bl_option_mod')

call umPrint('- - - - - - end of namelist - - - - - -',                        &
    src='bl_option_mod')

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine print_nlist_run_bl


end module bl_option_mod
