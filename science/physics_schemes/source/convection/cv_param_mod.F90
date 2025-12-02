! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
!  Global data module for switches/options concerned with convection.

module cv_param_mod

  ! Description:
  !   Module containing parameters used by the convection code.
  !
  ! Method:
  !   Default values have been declared where appropriate.
  !
  !   Any routine wishing to use these options may do so with the 'Use'
  !   statement.
  !
  ! Code Owner: Please refer to the UM file CodeOwners.txt
  ! This file belongs in section: Convection
  !
  ! Code Description:
  !   Language: FORTRAN 90
  !   This code is written to UMDP3 v7.4 programming standards.
  !
  ! Declarations:

use conversions_mod, only: rsec_per_day  ! Get seconds per day, to easily
                                         ! express constants in units per day
                                         ! but store in units per second.

use um_types, only: real_umphys, r_bl

implicit none
save

!===========================================================================
! Integer parameters to remove use of magic numbers
!===========================================================================

! Basis for 2d convective cloud amount
!   Total Condensed Water(TCW) original 4A
!   LES mb/wsc scalings (Grant and Lock, 2004) (shallow)
!   Surface rain from respective convective cloud type (Dp and Md-level)
!   Average mass flux based (Dp and Md-level)
integer, parameter :: total_condensed_water = 0
integer, parameter :: grant_lock_over       = 1 ! Sh cu with overlap (CCRad)
integer, parameter :: grant_lock            = 2 ! Sh cu no overlap (CCRad)
integer, parameter :: srf_precip            = 1 ! mid/deep cnv (CCRad)
integer, parameter :: avg_mass_flux         = 2 ! avg mass flux: md/dp cnv

! Convective cloud decay
integer, parameter :: rad_decay_off           = 0
integer, parameter :: rad_decay_full_timestep = 1
integer, parameter :: rad_decay_conv_substep  = 2 ! CCRad only


! Convective cloud decay timescale
integer, parameter :: cld_life_constant = 0
integer, parameter :: cld_life_func_hgt = 1 ! CCRad only

! For all the options below except 3), the anvil base is at the freezing
! Level
integer, parameter :: anv_pressure      = 0
integer, parameter :: anv_height        = 1
integer, parameter :: anv_model_levels  = 2
integer, parameter :: anv_limited_pressure_depth = 3

! Energy correction for 6a convection
integer, parameter :: method_en_rho=1     ! Correct energy but not water on
                                          ! the rho grid
integer, parameter :: method_en_mx_rho=2  ! Correct water and then energy on
                                          ! the rho grid
integer, parameter :: method_en_qx_p=3    ! Correct water and then energy on
                                          ! the pressure level grid. Not
                                          ! recommended for normal use a
                                          ! useful check when developing
                                          ! scheme

! Termination condition
integer, parameter :: term_orig  = 0      ! Original termination condition
integer, parameter :: term_simp  = 1      ! Simplified termination condition
integer, parameter :: term_undil = 2      ! Simplified plus test on undilute
                                          ! parcel buoyancy

! Mid-level trigger options
integer, parameter :: mtrig_ntmlplus2 =0
integer, parameter :: mtrig_ntmlplus1 =1
integer, parameter :: mtrig_ntml      =2
integer, parameter :: mtrig_surface   =3

! Options for the treatment of phase changes in convective precipitation
integer, parameter :: pr_melt_frz_inst= 0 ! Instantaneously freeze rain
                                          ! and melt snow. Mixed phase precip
                                          ! not allowed.
integer, parameter :: pr_melt_tdep    = 1 ! Snow is melted at a rate dependent
                                          ! on the temperature.
integer, parameter :: pr_melt_frz_tdep= 2 ! Snow is melted and rain is frozen
                                          ! at rates dependent on the
                                          ! temperature.

! Mid-level perturbation options
integer, parameter :: md_pert_orig  = 0   ! Original perturbation
integer, parameter :: md_pert_efrac = 1   ! Perturbation partition between
                                          ! temperature and humidity based
                                          ! on specified efrac
integer, parameter :: md_pert_bowen = 2   ! Perturbation partition between
                                          ! temperature and humidity based
                                          ! on surface Bowen ratio

! Shallow scheme cloud base closure options
integer, parameter ::  sh_wstar_closure = 0 ! standard closure on wstar
integer, parameter ::  sh_grey_closure  = 1 ! grey zone closure

! Parameters used to relate cca_2d of cld to its precipitation rate
real(kind=real_umphys), parameter    :: a_land = 0.3
real(kind=real_umphys), parameter    :: a_sea  = 0.3
real(kind=real_umphys), parameter    :: b_land = 0.025
real(kind=real_umphys), parameter    :: b_sea  = 0.025

! Parameter for mass flux based cca_2d
real(kind=real_umphys), parameter    :: w_cca  = 1.0
                                    ! assumed vertical velocity in cca_2d
                                    ! calculation (m/s)

! Application of convective cloud anvils
real(kind=real_umphys), parameter :: deep_dp = 50000.0
                                      ! Min. depth for anvil criteria(Pa)

! Critical depth of cloud for the formation of
! convective precipitation over sea (m)
real(kind=real_umphys), parameter :: critdsea = 1.5e3

! critical depth of cloud for the formation of convective
! precipitation over land (m)
real(kind=real_umphys), parameter :: critdlnd = 4.0e3

! critical depth of a glaciated cloud for the formation of
! convective precipitation (m)
real(kind=real_umphys), parameter :: critdice = 1.0e3

! Parcel ascent in diagnosis, cloud water condensate
real(kind=r_bl), parameter :: qlcrit = 1.0e-3_r_bl   ! critical cloud water

! Timestep frequency for calling convection - hardwired to call every timestep
integer,parameter :: a_conv_step = 1

! Threshold magnitude of convective temperature tendency to determine whether
! convection is "active" on a given model-level (used for the conv_prog_precip
! experimental prognostic field).
real(kind=real_umphys)                                                         &
    ,parameter :: dthetadt_conv_active_threshold = 0.001 / rsec_per_day
                               ! Set to 0.001 K day-1

! Minimum bound applied to the convective precipitation rate when updating the
! conv_prog_precip experimental prognostic field / kg m-2 s-1
real(kind=real_umphys),parameter :: conv_prog_precip_min_threshold = 1.0e-5

!============================================================================
! Parcel ascent
!============================================================================

  ! coefficients used in calculation of entrainment rate
real(kind=real_umphys), parameter :: ae1     = 1.0        ! Not used
real(kind=real_umphys), parameter :: ae2     = 1.5

! Reference depth for variable entrainment and mixing detrainment for deep
! convection. It is the depth at which the variable entrainment
! and m.det will be the same as the original rates.
real(kind=real_umphys), parameter :: refdepth_dp = 8000.0

! Reference saturated humidity in kg/kg for variable entrainment.
real(kind=real_umphys), parameter :: refqsat = 0.02

! minimum excess buoyancy to continue parcel ascent (K)

real(kind=real_umphys),parameter:: xsbmin = 0.2       ! Used in mid scheme

! initial excess potential temperature (k) and mixing ratio
! (kg/kg) for deep convection
real(kind=real_umphys), parameter :: thpixs_deep= 0.2
real(kind=real_umphys), parameter :: qpixs_deep =0.0

! initial excess potential temperature (k) and mixing ratio
! (kg/kg) for shallow convection
real(kind=real_umphys), parameter :: thpixs_shallow = 0.2
real(kind=real_umphys), parameter :: qpixs_shallow  = 0.0

! initial excess mixing ratio
! (kg/kg) for mid-level convection
real(kind=real_umphys), parameter :: qpixs_mid =0.0

! Minimum parcel buoyancy/layer thickness (K/Pa)
real(kind=real_umphys), parameter :: mparb = 1.0            ! Used in mid scheme

! Constants to determine initial convective mass flux from parcel buoyancy
! Deep convection
real(kind=real_umphys), parameter :: c_deep = 5.17e-4       ! No longer used
real(kind=real_umphys), parameter :: d_deep = 0.0           ! No longer used

! Shallow convection
real(kind=real_umphys), parameter :: c_shallow = 5.17e-4    ! No longer used
real(kind=real_umphys), parameter :: d_shallow = 0.0        ! No longer used

! Mid convection
real(kind=real_umphys), parameter :: c_mid = 5.17e-4
real(kind=real_umphys), parameter :: d_mid = 0.0

! limits on the initial convective parcel perturbations
real(kind=r_bl), parameter :: max_diag_thpert  =  2.0_r_bl
real(kind=real_umphys), parameter :: max_dp_thpert    =  2.0
real(kind=real_umphys), parameter :: min_dp_thpert    = -2.0
real(kind=real_umphys), parameter :: max_sh_thpert    =  2.0
real(kind=real_umphys), parameter :: min_sh_thpert    = -2.0
real(kind=real_umphys), parameter :: max_dp_qpert_fac =  0.2
real(kind=real_umphys), parameter :: max_sh_qpert_fac =  0.2

! mparfl = 1E-3 * minimum parcel buoyancy * mass flux parameter c
real(kind=real_umphys), parameter :: mparfl = 1.0e-3 * 1.0 * 3.33e-4

! Difference in potential temperature between levels above which the
! atmosphere is assumed to be too stable to convect (K)
real(kind=real_umphys), parameter :: delthst = 0.5

! Entrainment coefficient

real(kind=real_umphys), parameter :: entcoef =3.0

!============================================================================
! Convective closure
!============================================================================

  ! Coefficient relating sub-cloud convective velocity scale to cumulus
  ! mass flux for shallow convection
real(kind=real_umphys), parameter :: c_mass=0.03

! Tuneable factor used in denominator of W_CAPE timescale calculation
real(kind=real_umphys), parameter :: wcape_fac = 3.0

! Parameter governing the speed of transition from parametrized to
! resolved convection, for use with cldbase_opt_sh = sh_grey_closure
real(kind=real_umphys) :: beta_cu = 0.15

! CAPE closure based on large-scale w
real(kind=real_umphys), parameter :: a_cape = 3600.0*0.08
                                         ! value from CASCADE fit
real(kind=real_umphys), parameter :: b_cape = -0.7
                                        ! value from CASCADE fit

! Deep cloud closure based on wstar and large-scale w
! wup_cb   = a_cb * wstar                           (units m/s)
! sigma_cb = b_cb + c_cb * wls                      (fractional area)
! mf_cb    = rho_cb * sigma_cb * wup_cb             (Kg/m2/s)
! massflux_cb = g * rho_cb * sigma_cb * wup_cb         (Pa/s)
real(kind=real_umphys), parameter :: a_cb = 2.6
                                       ! from fit to CRM simulations
real(kind=real_umphys), parameter :: b_cb = 0.02         !
real(kind=real_umphys), parameter :: c_cb = 0.2          !

!============================================================================
! Downdraught and evaporation below cloud base calculations
!============================================================================

  ! Coefficients used in calculation of downdraught entrainment
  ! rates
real(kind=real_umphys), parameter :: ddcoef1 = 1.8e6
real(kind=real_umphys), parameter :: ddcoef2 = 3.0

! Thickness level used in calculation of mixing detrainment for
! downdraught  (pa)
real(kind=real_umphys), parameter :: det_lyr = 10000.0

! exponents used in calculation of evaporation of liquid
real(kind=real_umphys), parameter :: p_lq1 = 0.52
real(kind=real_umphys), parameter :: p_lq2 = 0.67

! exponents used in calculation of evaporation of ice
real(kind=real_umphys), parameter :: p_ice1 = 0.55
real(kind=real_umphys), parameter :: p_ice2 = 0.76

! exponents and constants associated with density term in
! evaporation of liquid
real(kind=real_umphys), parameter :: rho_lqp1 = 0.26
real(kind=real_umphys), parameter :: rho_lqp2 = 0.59
real(kind=real_umphys), parameter :: rho_lqa  = 108.80
real(kind=real_umphys), parameter :: rho_lqb  = 830.73

! exponents and constants associated with density term in
! evaporation of ice
real(kind=real_umphys), parameter :: rho_icp1 = 0.28
real(kind=real_umphys), parameter :: rho_icp2 = 0.63
real(kind=real_umphys), parameter :: rho_icea = 1569.52
real(kind=real_umphys), parameter :: rho_iceb = 32069.02

! constants used in quadratic formula for evaporation of liquid
real(kind=real_umphys), parameter :: lq_a = 2.008e-9
real(kind=real_umphys), parameter :: lq_b = -1.385e-6
real(kind=real_umphys), parameter :: lq_c = 2.424e-4

! constants used in quadratic formula for evaporation of ice
real(kind=real_umphys), parameter :: ice_a = -5.2e-9
real(kind=real_umphys), parameter :: ice_b = 2.5332e-6
real(kind=real_umphys), parameter :: ice_c = -2.911e-4
real(kind=real_umphys), parameter :: ice_d = 1.7405e-5
                                      ! value of A(T,p) at T=243.58
                                      ! See Gregory 1995 in UM DOC 27

! Scaling factor used to calculate area occupied by precip below cloud-base
! when no downdraught, used in evaporation and melting calculations
real(kind=real_umphys), parameter :: precip_area_fac = 1.0
                                         ! precip area = CCA * precip_area_fac

! Scaling factor used to calculate area occupied by downdraught,
! used in evaporation calc
real(kind=real_umphys), parameter :: dd_area_fac = 0.5
                                         ! dd area = CCA * dd_area_fac

! downdraught precipitation transfer efficiency factor
real(kind=real_umphys), parameter :: ddptef = 2.0

! Characteristic temperature at which rain should freeze
real(kind=real_umphys), parameter :: t_frez_rain = 270.15

!============================================================================
! Convective momentum transport (CMT) calculations
!============================================================================

! Deep turbulent CMT scheme
real(kind=real_umphys), parameter :: dp_cmt_gamma = 1.63
real(kind=real_umphys), parameter :: dp_cmt_delta = 2.0
real(kind=real_umphys), parameter :: top_press    = 15000.0




! Parameters specifying methods in calculation of w_eqn
integer, parameter :: SimpsonWiggert = 1
integer, parameter :: BuoyScale = 2

!------------------------------------------------------------------------------
! Current parameter settings for w-eqn options
! Currently hardwired for development
!------------------------------------------------------------------------------

real(kind=real_umphys),    parameter :: w2pi           = 1.0
                                              ! Initial w^2 at cloud base
real(kind=real_umphys),    parameter :: gamma_in_w_eqn = 0.5
                                              ! Virtual mass coefficient
real(kind=real_umphys),    parameter :: cumulus_r      = 100.0
                                              ! Radius of cumulus tower ??
real(kind=real_umphys),    parameter :: drag_coeff     = 1.25
                                              ! Aerodynamic drag coefficient
                                              ! for solid spheres
real(kind=real_umphys),    parameter :: k2_const       = 0.71
                                              ! Entrainment coefficient
real(kind=real_umphys),    parameter :: fix_alpha      = 0.0001
                                              ! Value for alpha if constant
real(kind=real_umphys),    parameter :: gamma_b        = 130.0
                                              ! Scaling value for buoyancy
real(kind=real_umphys),    parameter :: NegBuoyMinW    = 0.1
                                              ! Min. threshold of buoyancy

integer, parameter :: watload_opt    = 2      ! How to apply water loading
                                              ! for SimpsonWiggert W calc
                                              ! 1) SimpsonWiggert
                                              ! 2) As in UM
integer, parameter :: NegBuoyOpt     = 0      ! Option to decide what to do
                                              ! when encountering -ve buoyancy
integer, parameter :: alpha_opt      = 2      ! How to set alpha
                                              ! 1) As with SimpsonWiggert
                                              ! 2) Us UM entrainment values
                                              ! 3) Set as constant value
integer, parameter :: wCalcMethod    = 1      ! Sets method to calculate w
                                              ! 1) SimpsonWiggert
                                              ! 2) Buoyancy Scaling


!------------------------------------------------------------------------------
! end w-eqn options
!------------------------------------------------------------------------------

!============================================================================
! Diagnostics passed to the perturbation forecast (PF) model in VAR
! This parameter will not affect model evolution on its own but will
! alter the evolution through the analysis cycle if l_fix_conv_diags_var=True
! because it affects the mass flux diagnostic that is passed to VAR.
! The chosen value works well but may not be optimal. It is possible that the
! optimal value of this parameter will change slightly with vertical
! resolution.
!============================================================================

real(kind=real_umphys), parameter :: max_mf_fall = 0.3
                                              ! The maximum fractional fall
                                              ! in mass flux between k and k+1
                                              ! before the mass flux is masked
                                              ! out

!============================================================================
! Prognostic based initial perturbation
!============================================================================
real(kind=real_umphys), parameter :: min_pert_scale = 1.0
                                              ! The minimum perturbation
                                              ! scaling
real(kind=real_umphys), parameter :: logp_min = -5.0
                                              ! log10 of the precipitation rate
                                              ! where the minimum perturbation
                                              ! scaling is applied
real(kind=real_umphys), parameter :: logp_max = -3.0
                                              ! log10 of the precipitation rate
                                              ! where the maximum perturbation
                                              ! scaling is applied

!============================================================================
! Convective cold pools
!============================================================================

! convective cold pool switch options
integer, parameter :: ccp_off = 0
! "propagating" version
integer, parameter :: ccp_prop = 1
! "single column" version
integer, parameter :: ccp_sc = 2

real(kind=real_umphys), parameter :: g_max = 1.0
! max. allowable initial reduced gravity (m s^-2)
real(kind=real_umphys), parameter :: h_max = 750.0
! max. allowable initial height (m)

! parameter for cold-pool modification of entrainment
! Linear dependence on ccp_strength with an exponent of 1.0, etc.
real(kind=real_umphys), parameter :: ccp_exp = 0.25  ! exponent of ccp_strength

! coefficient in height decay of cold-pool parcel buoyancy perturbations
real(kind=real_umphys), parameter :: ccp_h_coef = 0.1  ! (dimensionless)

! tuning coefficient used in cold-pool strength calculation
real(kind=real_umphys), parameter :: strength_coeff = 3.0e2  ! (m^2 s^-2)

! tuning coefficients used in fb_surf perturbation
real(kind=real_umphys), parameter :: ccp_fb_coeff  =   0.01 ! (units as fb_surf)
real(kind=real_umphys), parameter :: ccp_fb_offset =   0.2  ! (dimensionless)

! value to initialise array when the ccp strength in unset
real(kind=real_umphys), parameter :: ccp_strength_unset_value = 0.0

! height over which to average background winds (fixed at present)
real(kind=real_umphys), parameter :: z_ccp_wind = 300.0 ! (m)


end module cv_param_mod
