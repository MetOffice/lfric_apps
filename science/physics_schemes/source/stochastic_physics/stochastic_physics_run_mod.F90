! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

module stochastic_physics_run_mod

!  Global data module for switches/options concerned with stochastic physics
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: stochastic_physics

  ! Description:
  !   Module containing runtime logicals/options used by the SKEB2
  !   and RP2 schemes.
  !   and the stochastic_physics_run_setup subroutine to control logic of
  !   selected options.

  ! Method:
  !   All switches/options which are contained in the &RUN_Stochastic
  !   sub-namelist in the CNTLATM control file are declared in this module.
  !   Default values have been declared where appropriate, but mostly
  !   rmdi values are set to force users to specify values through Rose.

  !   Any routine wishing to use these options may do so with the 'Use'
  !   statement.

  !   Note, flags for other stochastic physics routines, which were
  !   embedded within other namelists, have now been moved here.

use missing_data_mod, only: rmdi, imdi
use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

use errormessagelength_mod, only: errormessagelength

! npft_max is defined in Jules code which the Utils do not depend on.
! fieldcalc uses this module, but not the arrays being declared as size
! npft_max. So rather than include a dependency on Jules in the Utils, the
! requirement to include npft_max is excluded from Utils builds.
#if !defined(UTILIO)
use max_dimensions, only: npft_max
use jules_vegetation_mod, only: l_triffid, l_phenol
#endif

use um_types, only: real_umphys

implicit none
!======================================================================
! Logical switches set from RUN_Stochastic namelist via UI
!======================================================================

logical :: l_skeb2  = .false.                                                  &
           ! switch controls use of SKEB2 scheme
,           l_skeb2_psicdisp = .false.                                         &
           ! true = incl streamfunction modulation by convection
,           l_skeb2_psisdisp = .false.                                         &
           ! true = incl streamfunction modulation by SL-advection
,           l_skeb2_velpot = .false.                                           &
           ! true = Calc divergent-wind incr from VelPot forcing
,           l_rp2 = .false.                                                    &
           ! switch controls use of RP2 scheme
,           l_skebsmooth_adv = .false.                                         &
           ! true = Perform advanced smoothing of energy diss fields
,           l_skebprint = .false.                                              &
           ! true = Print global KE backscattered at each timestep
,           l_skeb2_conv_disp_mod=.false.                                      &
           ! true= A horizontal resolution factor modulates the
           ! convective dissipation rate.


,           l_x_eq_sin_x= .false.                                              &
           ! true= Activates the approximation sin(x)=x for the
           ! computation of Lengrende Polynomials
,           l_rp2_cycle_out = .false.                                          &
           ! true = Write RP parameters to a file after 12 hours
,           l_rp2_cycle_in = .false.                                           &
           ! true = Initialise RP parameters from a specified file

!=========================================================================
! Logical switches for Stochastic Peturbation of Tendencies scheme (SPT)
!    (true: Included, false: Not included)
!=========================================================================

 ,           l_spt = .false.                                                   &
            ! Switch controls use of SPT scheme
 ,           l_spt_cfl=.false.                                                 &
             ! Switch to activate the CFL filtering to those
             ! points that breach CFL criteria
 ,           l_spt_rain = .false.                                              &
             ! Switch to add LS rain (microphysics) tendency to SPT
 ,           l_spt_rad= .false.                                                &
             ! Switch to add Radiation tendency to SPT
 ,           l_spt_gwd= .false.                                                &
             ! Switch to add Gravity wave drag tendency to SPT
 ,           l_spt_conv= .false.                                               &
             ! Switch to add convection tendency to SPT
 ,           l_spt_conv_mom= .false.                                           &
             ! Switch to add tendency from convective momentum
             ! transport to SPT
 ,           l_spt_qcons = .false.                                             &
             ! Conserve total water in the column
 ,           l_spt_mse = .false.
             ! Conserve moist static energy

! Define numerical and convective scheme type, only one can be active
! at any time.
integer :: type_cape = 4 ! Calc Conv Dissipation using CAPE
integer :: type_mflx = 5 ! Calc Conv Dissipation using Mass-Flux
integer :: skeb2_cdisp = imdi ! type of streamfn modulation by conv
integer :: type_smag = 6 ! Calc num dissipation using Smagorinsky
integer :: type_bihm = 7 ! Calc num dissipation using bi-harmonic
integer :: skeb2_sdisp = imdi ! type of numerical dissipation scheme

!======================================================================
! Integer parameters to identify type of model run
!======================================================================
integer, parameter :: firsttimestep_true  = 1
         ! This is the first timestep in a new run (NRUN)
integer, parameter :: firsttimestep_false = 2
         ! This is after the first timestep
integer, parameter :: firsttimestep_crun  = -1
         ! This is the first timestep in a continuation run (CRUN)

!======================================================================
! Integers to indicate which stochastic physics information is
! output to the dump header
!======================================================================

integer, parameter :: stph_seed_present          = 1
integer, parameter :: stph_spt_data_present      = 2
integer, parameter :: stph_skeb2_data_present    = 4
integer, parameter :: stph_rp2_data_present      = 8

integer :: stph_header_flag = stph_seed_present
           ! Flag to pass to ih_stochastic_flag indicating which
           ! stochastic physics information is output to the dump header.
           ! For each stochastic physics scheme in use, the parameter
           ! stph_<scheme>_data_present is added to stph_header_flag,
           ! resulting in a single integer value that can be used to determine
           ! which schemes have data stored in the dump header.
           ! It is initialised to stph_seed_present since the stochastic
           ! seed is stored if any stochastic physics scheme is in use.
           ! Data for the stochastic BL perturbation scheme is stored as
           ! a prognostic so is not included here.
integer :: stph_spt_data_check = 0
integer :: stph_skeb2_data_check = 0
integer :: stph_rp2_data_check = 0

integer, parameter :: rp_max = 25  ! Max number of RP fields
                                   ! Used to set the size of the array
                                   ! for storing RP values in the dump

integer, parameter :: n_rp = 3  ! Length of RP entry
integer, parameter :: rp_min_idx = 1  ! Index of RP value for min_rp
integer, parameter :: rp_idx     = 2  ! Index of RP value for rp
integer, parameter :: rp_max_idx = 3  ! Index of RP value for max_rp

#if defined(UTILIO)
integer, parameter :: npft_max = 1
logical :: l_triffid = .false.
logical :: l_phenol = .false.
#endif

!======================================================================
! Integer options set from RUN_Stochastic namelist via UI
!======================================================================
!
!   stph_n1 and stph_n2 define the range of spherical harmonic orders over
!     which backscatter is applied (remember: N(N+1)/R*R is the
!     effective wavenumber squared). For ENDGame the suggested range
!     is [20;60] - horizontal wavelengths in the range 500 km to 2000 km.
!
!   skeb2_toplev defines the top level of the streamfunction modulation
!     field. This will change with vertical resolution. Typically it is
!     set to a level near 18km to avoid making changes above the stratosphere
!   skeb2_botlev defines the bottom level of the streamfunction
!     modulating field (we may want to avoid changes in the boundary
!     layer, certainly in level one)
!   nsmooth is the number of 1-2-1 spatial smoothing iterations (5)
!   rhcrit_ref_level is usually level 3

integer :: stphseed = imdi       ! Control variable for random seed options:
           ! 0 => Use ensemble member and date/time of dump
           ! 1 => Use date/time from computer clock
           ! 2 => Use seed/coefficient values stored in dump file

integer :: stph_n1=imdi          ! minimum wavenumber for backscatter
integer :: stph_n2=imdi          ! maximum wavenumber for backscatter
integer :: rp2_callfreq = imdi   ! RP2 calling frequency (in seconds)
integer :: rp2_cycle_tm = imdi   ! time at which to write out RP2 parameters
                                 ! (in seconds) when l_rp2_cycle_out = true

integer :: i_rp_scheme = imdi    ! Switch to specify RP scheme option:
           ! 0 => RP2 scheme
integer, parameter :: i_rp2 = 0
           ! 1 => RP2b scheme (updated algorithm and additional parameters)
integer, parameter :: i_rp2b = 1

integer :: skeb2_toplev=imdi     ! Top level of SKEB2 calculations
integer :: skeb2_botlev=imdi     ! Bottom level of SKEB2 calculations
integer :: nsmooth = imdi        ! Iteration count for spatial smoothing
integer :: rhcrit_ref_level=imdi ! RHCrit reference level for RP2
integer :: offx_stph = imdi      ! Halo size used in spatial smoothing
integer :: offy_stph = imdi      ! Halo size used in spatial smoothing
                                 ! Both set to nsmooth in sthp_setup
integer :: ran_max = imdi        ! number of independent AR1 processes in RP2
integer :: ran_count = imdi      ! counter for random number array in RP
integer :: stph_nens = imdi      ! ensemble member number
character(len=8) :: ens_member   ! ensemble member (8-character string)

! ---- SPT integers ---
integer :: spt_top_tap_lev= imdi ! Top level of spt calculations.
                                 ! Including tapering zone
integer :: spt_toplev= imdi      ! Top level of spt calculations.
                                 ! Excluding tapering zone
integer :: spt_bot_tap_lev= imdi ! Bottom level of spt calculations.
                                 ! Including tapering zone
integer :: spt_botlev= imdi      ! Bottom level of spt calculations (min 2)
                                 ! Excluding tapering zone
integer :: nsmooth_spt = imdi    ! Iteration count for spatial smoothing
                                 ! of SPT increments
integer :: offx_spt, offy_spt    ! Halo size used in spatial smoothing
                                 ! for SPT

 !======================================================================
 ! Real values set from  RUN_Stochastic namelist via GUI
 !======================================================================
real(kind=real_umphys) :: tau_skeb = rmdi
                                 ! decorrelation time for SKEB AR1 (in secs)
real(kind=real_umphys) :: tau_spt  = rmdi
                                 ! decorrelation time for SPT AR1 (in secs)
real(kind=real_umphys) :: tot_backscat = rmdi
                                 ! global-mean rate of energy backscatter
                                 ! in m**2 s**(-3) (usually 1.0e-4)
real(kind=real_umphys) :: br = rmdi
                                 ! backscatter ratio (frac of assumed
                                 ! dissipated energy - usually 0.2)
real(kind=real_umphys) :: sdispfac = rmdi
                                 ! Multiplication factor for numerical diss
                                 ! (determined empirically as 2.0)
real(kind=real_umphys) :: cdispfac = rmdi
                                 ! Multiplication factor for convection diss
                                 ! (determined empirically as 1.0)
real(kind=real_umphys) :: alphac = rmdi
                                 ! Updraught proportion of gridbox
                                 ! (0.2% typical)

real(kind=real_umphys) :: rp2_decorr_ts = rmdi
                                  ! RP2 de-correlation timescale (in seconds)

! RHCrit minimum and maximum default values (varies with resolution)
real(kind=real_umphys) :: rhcrit_max          = rmdi
real(kind=real_umphys) :: rhcrit_min          = rmdi

!------------------------------------

! Stochastic physics ice fallspeed multiplier m_ci,
! Min / default / max values:
real(kind=real_umphys) :: m_ci_rp(n_rp) = rmdi ! suggested values 0.6, 1.0, 1.4

! Rain particle size distribution, x1r
! Min / default / max values:
real(kind=real_umphys) :: x1r_rp(n_rp) = rmdi
                            ! suggested values 0.07, 0.22, 0.52

! Surface droplet number, ndrop_surf
! Min / default / max values:
real(kind=real_umphys) :: ndrop_surf_rp(n_rp) = rmdi ! suggested values
                                   ! 2.0e+07, 7.5e+07, 10.0e+07

! Snow fallspeed parameter snow_fspd_rp,
! Min / default / max values:
real(kind=real_umphys) :: snow_fspd_rp(n_rp) = rmdi
! suggested values 7.2, 12.0, 16.8

! Ice fallspeed parameter ice_fspd_rp,
! Min / default / max values:
real(kind=real_umphys) :: ice_fspd_rp(n_rp) = rmdi
! suggested values 3600000.0, 6000000.0, 8400000.0

! Mixed-phase overlap factor mpof, used in CASIM
! Min / default / max values:
real(kind=real_umphys) :: mpof_rp(n_rp) = rmdi
                          ! suggested values 0.1, 0.5, 0.6

! Fixed cloud number, used in CASIM
! Min / default / max values:
real(kind=real_umphys) :: fxd_cld_num_rp(n_rp) = rmdi
                          ! suggested values 70.0e+06, 150.0e+06, 300.0e+06

! Lagrangian structure-function constant mp_czero,
! Min / default / max values:
real(kind=real_umphys) :: mp_czero_rp(n_rp) = rmdi
                            ! suggested values 4.0, 10.0, 10.0

!------------------------------------

! ROSE - BOUNDARY LAYER stochastic items moved here from bl_option_mod
! Maximum and minimum values for the STPH_RP scheme
! Boundary Layer

! Neutral mixing length, par_mezcla
! Min / default / max values:
real(kind=real_umphys) :: par_mezcla_rp(n_rp) = rmdi
                                   ! suggested values 0.03, 0.15, 0.45

! Flux profile parameter, g0
! Min / default / max values:
real(kind=real_umphys) :: g0_rp(n_rp) = rmdi ! suggested values 5.0, 10.0, 20.0

! Charnock parameter for roughness length of sea points
! Max and min values:
real(kind=real_umphys) :: charnock_max= rmdi       ! suggested value = 0.026
real(kind=real_umphys) :: charnock_min= rmdi       ! suggested value = 0.01

! Minimum mixing length, lambda_min
! Min / default / max values:
real(kind=real_umphys) :: lambda_min_rp(n_rp) = rmdi
                                   ! suggested values 8.0, 40.0, 120.0

! Critical Richardson number, Ricrit
! Min / default / max values:
real(kind=real_umphys) :: ricrit_rp(n_rp) = rmdi
                               ! suggested values 0.25 / 1.0 / 1.0

! Entrainment parameter A1
! Min / default / max values:
real(kind=real_umphys) :: a_ent_1_rp(n_rp) = rmdi
                                ! suggested values 0.1 / 0.23 / 0.4

! Shear entrainment parameter, A1_shr
! Max and min values (default not used - linked to a_ent_1)
real(kind=real_umphys) :: a_ent_shr_rp_max = rmdi  ! suggested value = 5.0
real(kind=real_umphys) :: a_ent_shr_rp  = rmdi     ! suggested value = 1.6

! Velocity scale parameter, g1
real(kind=real_umphys) :: g1_rp(n_rp) = rmdi ! suggested values 0.5 / 0.85 / 1.5
! Turbulent mixing in the convective BL, cbl_mix_fac
! Min / default / max values:
real(kind=real_umphys) :: cbl_mix_fac_rp(n_rp) = rmdi
                                    ! suggested values 0.0, 0.5, 1.0

! Smagorinsky coefficient, cs_rp
! Min / default / max values:
real(kind=real_umphys) :: cs_rp(n_rp) = rmdi ! suggested values 0.01, 0.2, 1.0

! Fractional cloud height reached by local BL depth calulcation,
! zhloc_depth_fac_rp
! Min / default / max values:
real(kind=real_umphys) :: zhloc_depth_fac_rp(n_rp) = rmdi
                                        ! suggested values 0.1, 0.5, 1.0

! Land-surface parameters for use with i_rp_scheme == i_rp2b

! Drag coefficient for orographic form drag
! Default, max and min values for orog_drag_param_rp
! Min / default / max values:
! suggested values 0.001, 0.15, 0.2
real(kind=real_umphys) :: orog_drag_param_rp(n_rp) = rmdi
!
! These are arrays corresponding to plant function types.
!
! Default, max and min values for JULES parameter z0v
real(kind=real_umphys) :: z0v_rp(npft_max)
real(kind=real_umphys) :: z0v_rp_max(npft_max)
real(kind=real_umphys) :: z0v_rp_min(npft_max)
data z0v_rp / npft_max * rmdi /
data z0v_rp_max / npft_max * rmdi /
data z0v_rp_min / npft_max * rmdi /
! JULES roughness length for bare soil z0_nvg(2), denoted
! here as z0_soil
! Varies with z0v
! Min / default / max values:
real(kind=real_umphys) :: z0_soil_rp(n_rp) = rmdi
! JULES roughness length for urban canyon and roof
! tiles are perturbed by multiplying by a random factor
! Varies with z0v
! Default, max and min values for z0_urban_mult
real(kind=real_umphys) :: z0_urban_mult_rp(n_rp) = rmdi
! Default, max and min values for JULES parameter z0hm_pft
! Varies with z0v
real(kind=real_umphys) :: z0hm_pft_rp(npft_max)
real(kind=real_umphys) :: z0hm_pft_rp_max(npft_max)
real(kind=real_umphys) :: z0hm_pft_rp_min(npft_max)
data z0hm_pft_rp / npft_max * rmdi /
data z0hm_pft_rp_max / npft_max * rmdi /
data z0hm_pft_rp_min / npft_max * rmdi /
! JULES parameter for bare soil z0hm_nvg(2), denoted
! here as z0hm_soil
! Varies with z0hm_pft
! Min / default / max values:
real(kind=real_umphys) :: z0hm_soil_rp(n_rp) = rmdi
! LAI multiplier - used to scale lai
! Default, max and min values for lai_mult_rp
! Note lai_mult_rp default value should be 1.0
real(kind=real_umphys) :: lai_mult_rp(npft_max)
real(kind=real_umphys) :: lai_mult_rp_max(npft_max)
real(kind=real_umphys) :: lai_mult_rp_min(npft_max)
data lai_mult_rp / npft_max * rmdi /
data lai_mult_rp_max / npft_max * rmdi /
data lai_mult_rp_min / npft_max * rmdi /
! Surface vegetation tile albedo parameters (4 separate ones)
! Default, max and min values for JULES parameter alnir
real(kind=real_umphys) :: alnir_rp(npft_max)
real(kind=real_umphys) :: alnir_rp_max(npft_max)
real(kind=real_umphys) :: alnir_rp_min(npft_max)
data alnir_rp / npft_max * rmdi /
data alnir_rp_max / npft_max * rmdi /
data alnir_rp_min / npft_max * rmdi /
! Default, max and min values for JULES parameter alpar
real(kind=real_umphys) :: alpar_rp(npft_max)
real(kind=real_umphys) :: alpar_rp_max(npft_max)
real(kind=real_umphys) :: alpar_rp_min(npft_max)
data alpar_rp / npft_max * rmdi /
data alpar_rp_max / npft_max * rmdi /
data alpar_rp_min / npft_max * rmdi /
! Default, max and min values for JULES parameter omega
real(kind=real_umphys) :: omega_rp(npft_max)
real(kind=real_umphys) :: omega_rp_max(npft_max)
real(kind=real_umphys) :: omega_rp_min(npft_max)
data omega_rp / npft_max * rmdi /
data omega_rp_max / npft_max * rmdi /
data omega_rp_min / npft_max * rmdi /
! Default, max and min values for JULES parameter omnir
real(kind=real_umphys) :: omnir_rp(npft_max)
real(kind=real_umphys) :: omnir_rp_max(npft_max)
real(kind=real_umphys) :: omnir_rp_min(npft_max)
data omnir_rp / npft_max * rmdi /
data omnir_rp_max / npft_max * rmdi /
data omnir_rp_min / npft_max * rmdi /

!==============================================
! Stochastic forcing of theta in the BL options
!==============================================
! Switch to perturb theta: 0 = Off
!                          1 = with mag_pert_theta
!                          2 = using heat flux
!                          3 = theta and moisture
integer :: i_pert_theta = imdi
integer, parameter :: pert_theta_mag = 1
integer, parameter :: pert_theta_star = 2
integer, parameter :: pert_theta_and_moist = 3

! Switch to perturb all points
! If false, perturb only cumulus points
logical :: l_pert_all_points = .false.

! Switch to add vertical shape to perturbation profile if true
logical :: l_pert_shape = .false.

! Switch to specify time variation of perturbations
! 0 = random sequence
! 1 = time correlated
integer :: i_pert_theta_type = imdi
integer, parameter :: pert_theta_random_seq = 0
integer, parameter :: pert_theta_correl_seq = 1

! Maximum perturbation to theta (K)
real(kind=real_umphys) :: mag_pert_theta = rmdi

! Decorrelation timescale for theta perturbation (seconds)
real(kind=real_umphys) :: decorr_ts_pert_theta = rmdi

! Heights between which to pertub theta, converted to minlev_pert_theta
! and maxlev_pert_theta in setcona
real(kind=real_umphys) :: zmin_pert_theta = rmdi
real(kind=real_umphys) :: zmax_pert_theta = rmdi

! Model levels between which to perturb theta, calculated in setcona
integer :: minlev_pert_theta = imdi
integer :: maxlev_pert_theta = imdi

! Number of points to apply the same perturbation over
integer :: npts_pert_theta = imdi

! Switch to enable and the number of points from the domain edge to apply
! perturbations
logical :: l_only_pert_near_edge = .false.
integer :: npts_pert_from_edge = imdi

  !======================================================================
  ! Arrays from Convection required by SKEB2
  !======================================================================
real(kind=real_umphys), allocatable ::                                         &
      skeb2_up_flux(:, :, :)                                                   &
           !  updraught mass flux
,     skeb2_dwn_flux(:, :, :)                                                  &
           ! downdraught mass flux
,     skeb2_cape(:, :)                                                         &
           ! CAPE
,     Ymn(:,:,:)
           ! Matrix for the value of Legendre Polynomials


!========================================================================
! Standard deviation of the forcing pattern for the different
! parametrizations included in SPT.
!========================================================================

real(kind=real_umphys) ::                                                      &

   rain_std= rmdi                                                              &
           ! Sigma for rain ECMWF 0.72
,  rad_std= rmdi                                                               &
           ! Sigma for radiation ECMWF 0.33
,  gwd_std= rmdi                                                               &
           ! Sigma for Gravity wave drag ECMWF 0.52
,  conv_std= rmdi                                                              &
           ! Sigma for convection ECMWF 0.62
,  sd_orog_thres = rmdi                                                        &
           ! Threshold of orographic sd to apply the capping
,  psif_orog_thres = rmdi
           ! Threshold of psif above the one SPT pert are set to 0 over
           ! regions of high sd orog (specified by psif_orog_thres)

  !======================================================================
  ! Define Namelist RUN_Stochastic
  !======================================================================
namelist/run_stochastic/                                                       &
      l_skeb2, l_rp2, rp2_callfreq, i_rp_scheme, rp2_decorr_ts                 &
,     l_rp2_cycle_out, l_rp2_cycle_in, rp2_cycle_tm                            &
,     stph_n1, stph_n2, br, tot_backscat, tau_skeb, alphac                     &
,     l_skeb2_psicdisp, l_skeb2_psisdisp                                       &
,     sdispfac, cdispfac, skeb2_sdisp, skeb2_cdisp, nsmooth                    &
,     skeb2_toplev, skeb2_botlev, l_skeb2_velpot, ran_max                      &
,     rhcrit_ref_level, l_skebsmooth_adv, l_skebprint, stphseed                &
,     l_x_eq_sin_x, rhcrit_max, rhcrit_min                                     &
,     m_ci_rp, x1r_rp, ndrop_surf_rp                                           &
,     snow_fspd_rp, ice_fspd_rp, mpof_rp, fxd_cld_num_rp, mp_czero_rp          &
,     par_mezcla_rp, lambda_min_rp, cbl_mix_fac_rp, cs_rp                      &
,     zhloc_depth_fac_rp                                                       &
,     ricrit_rp, a_ent_1_rp, g1_rp, g0_rp                                      &
,     charnock_max, charnock_min, a_ent_shr_rp, a_ent_shr_rp_max               &
,     orog_drag_param_rp, lai_mult_rp, lai_mult_rp_max, lai_mult_rp_min        &
,     z0v_rp_max, z0v_rp, z0v_rp_min, z0_soil_rp, z0_urban_mult_rp             &
,     alnir_rp_max, alnir_rp, alnir_rp_min                                     &
,     alpar_rp_max, alpar_rp, alpar_rp_min                                     &
,     omega_rp_max, omega_rp, omega_rp_min                                     &
,     omnir_rp_max, omnir_rp, omnir_rp_min                                     &
,     z0hm_pft_rp_max, z0hm_pft_rp, z0hm_pft_rp_min, z0hm_soil_rp              &
,     l_skeb2_conv_disp_mod, l_pert_all_points, l_pert_shape                   &
,     decorr_ts_pert_theta, i_pert_theta_type                                  &
,     i_pert_theta, mag_pert_theta, zmin_pert_theta, zmax_pert_theta           &
,     npts_pert_theta, l_only_pert_near_edge, npts_pert_from_edge              &
! SPT namelist
,     l_spt, tau_spt, l_spt_rain, l_spt_rad, l_spt_gwd, l_spt_conv             &
,     l_spt_conv_mom, l_spt_cfl, spt_toplev, spt_botlev                        &
,     spt_bot_tap_lev, spt_top_tap_lev, nsmooth_spt                            &
,     rain_std, rad_std, gwd_std, conv_std,l_spt_qcons, l_spt_mse              &
,     sd_orog_thres, psif_orog_thres

  !======================================================================
  ! Array for spatial 1-2-1 smoothing
  !======================================================================
real(kind=real_umphys), allocatable, save ::                                   &
      mask_pdamp(:, :)                                                         &
           !  Array of pattern damping coefficients for SKEB2
,     mask_smooth(:, :)                                                        &
           !  Array of smoothing coefficients for SKEB2
,     mask_smooth_spt(:, :)
           !  Array of smoothing coefficients for SPT
!======================================================================
! logical switches not set in namelist
!======================================================================


integer  :: stphseed_unit = imdi
             ! Random seed file unit

!DrHook-related parameters
integer(kind=jpim), parameter, private :: zhook_in  = 0
integer(kind=jpim), parameter, private :: zhook_out = 1

character(len=*), parameter, private :: ModuleName='STOCHASTIC_PHYSICS_RUN_MOD'

contains

subroutine check_run_stochastic()

! Description:
!   Subroutine to apply logic controls and set control variables based on the
!   options selected in the run_stochastic namelist.

! Dr Hook Modules
use ereport_mod, only: ereport

use umPrintMgr, only: umPrint, umMessage, newline
implicit none

integer :: icode    ! error code for ereport
! local temporary arrays
character(len=errormessagelength)       :: cmessage      ! out error message
character(len=*), parameter  :: RoutineName='CHECK_RUN_STOCHASTIC'

real(kind=jprb)               :: zhook_handle

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! --------- Check that lai_mult_rp_min is not equal to zero in RP scheme
if (l_rp2 .and. i_rp_scheme == i_rp2b) then
  if ( l_triffid .or. l_phenol ) then
    icode   = 1
    write (cmessage,'(A)')                                                     &
      '  **********************************  '//             newline//         &
      'The RP2b scheme is on with either l_triffid'//        newline//         &
      'or l_phenol set to true.  The RP2b scheme should'//   newline//         &
      'not be used with either of these schemes.This is'//   newline//         &
      'because LAI is a prognostic in the triffid and'//     newline//         &
      'phenol schemes, and the RP2b scheme perturbs LAI'//   newline//         &
      'with a method that assumes it is not a prognostic.'
    call ereport(routinename, icode, cmessage)
  end if

  if ( any (abs(lai_mult_rp_min(:)) < tiny(0.0)) ) then
    icode   = 1
    write (cmessage,'(A)')                                                     &
      '  **********************************  '//             newline//         &
      'RP2b scheme is called with lai_mult_rp_min set'//     newline//         &
      'to zero.  This may result in the model failing.'
    call ereport(routinename, icode, cmessage)
  end if
end if
! --------- Check that appropriate convective dissipation scheme is active
if (l_skeb2_psicdisp) then
  if (skeb2_cdisp /= type_mflx .and. skeb2_cdisp /= type_cape) then
    write(umMessage,'(A)')'  **********************************  '
    call umPrint(umMessage,src='stochastic_physics_run_mod')
    write(umMessage,'(A)')'SKEB Conv dissipation is on, but neither'
    call umPrint(umMessage,src='stochastic_physics_run_mod')
    write(umMessage,'(A)')'CAPE or mass-flux scheme is selected.'
    call umPrint(umMessage,src='stochastic_physics_run_mod')
    icode   = 1
    write (cmessage,'(A)')'SKEB2 CAPE and mass-flux convective '               &
          // 'dissipation are both off, but l_skeb2_psicdisp is true.'
    call ereport(routinename, icode, cmessage)
  end if
end if
! --------- Check that appropriate numerical dissipation scheme is active
if (l_skeb2_psisdisp) then
  if (skeb2_sdisp /= type_smag .and. skeb2_sdisp /= type_bihm) then
    write(umMessage,'(A)')'  **********************************  '
    call umPrint(umMessage,src='stochastic_physics_run_mod')
    write(umMessage,'(A)')'SKEB Numerical dissipation is on, but neither'
    call umPrint(umMessage,src='stochastic_physics_run_mod')
    write(umMessage,'(A)')'Smagorinsky or Biharmonic scheme is selected.'
    call umPrint(umMessage,src='stochastic_physics_run_mod')
    icode   = 1
    write (cmessage,'(A)')'SKEB2 Biharmonic and Smagorinsky numerical '        &
          // 'dissipation are both off, but l_skeb2_psisdisp is true.'
    call ereport(routinename, icode, cmessage)
  end if
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine check_run_stochastic

subroutine print_nlist_run_stochastic()
use umPrintMgr, only: umPrint, ummessage
implicit none
real(kind=jprb) :: zhook_handle

character(len=*), parameter :: RoutineName='PRINT_NLIST_RUN_STOCHASTIC'
character(len=30) :: fmt_lsfc ! Format character for writing out
                              ! land-surface RP namelist values
character(len=30) :: fmt_rp ! Format character for writing out
                            ! atmos RP namelist values

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set write statement for land-surface random parameters based
! on the value of npft_max
write(fmt_lsfc,'("(A,",I0, "ES12.4)")') npft_max
write(fmt_rp,'("(A,",I0, "ES12.4)")') n_rp

call umPrint('Contents of namelist run_stochastic',                            &
    src='stochastic_physics_run_mod')

write(umMessage,'(A,L1)') 'l_skeb2 = ',l_skeb2
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,L1)') 'l_rp2 = ',l_rp2
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,I0)') 'rp2_callfreq = ',rp2_callfreq
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,I10)') 'i_rp_scheme = ',i_rp_scheme
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,ES12.4)') 'rp2_decorr_ts = ',rp2_decorr_ts
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,L1)') 'l_rp2_cycle_out = ',l_rp2_cycle_out
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,L1)') 'l_rp2_cycle_in = ',l_rp2_cycle_in
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,I10)') 'rp2_cycle_tm = ',rp2_cycle_tm
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,I0)') 'stph_n1 = ',stph_n1
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,I0)') 'stph_n2 = ',stph_n2
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,ES12.4)') 'br = ',br
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,ES12.4)') 'tot_backscat = ',tot_backscat
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,ES12.4)') 'tau_skeb = ',tau_skeb
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,ES12.4)') 'tau_spt = ',tau_spt
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,ES12.4)') 'alphac = ',alphac
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,L1)') 'l_skeb2_psicdisp = ',l_skeb2_psicdisp
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,L1)') 'l_skeb2_psisdisp = ',l_skeb2_psisdisp
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,ES12.4)') 'sdispfac = ',sdispfac
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,ES12.4)') 'cdispfac = ',cdispfac
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,I0)') 'skeb2_sdisp = ',skeb2_sdisp
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,I0)') 'skeb2_cdisp = ',skeb2_cdisp
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,I0)') 'nsmooth = ',nsmooth
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,I0)') 'skeb2_toplev = ',skeb2_toplev
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,I0)') 'skeb2_botlev = ',skeb2_botlev
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,L1)') 'l_skeb2_velpot = ',l_skeb2_velpot
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,I0)') 'ran_max = ',ran_max
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,I0)') 'rhcrit_ref_level = ',rhcrit_ref_level
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,L1)') 'l_skebsmooth_adv = ',l_skebsmooth_adv
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,L1)') 'l_skebprint = ',l_skebprint
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,L1)') 'l_x_eq_sin_x = ',l_x_eq_sin_x
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,I0)') 'stphseed = ',stphseed
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_rp) 'm_ci_rp = ',m_ci_rp
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_rp) 'x1r_rp = ',x1r_rp
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_rp) 'ndrop_surf_rp = ',ndrop_surf_rp
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_rp) 'snow_fspd_rp = ',snow_fspd_rp
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_rp) 'ice_fspd_rp = ',ice_fspd_rp
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_rp) 'mpof_rp = ',mpof_rp
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_rp) 'fxd_cld_num_rp = ',fxd_cld_num_rp
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_rp) 'mp_czero_rp = ',mp_czero_rp
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,ES12.4)') 'rhcrit_max = ',rhcrit_max
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,ES12.4)') 'rhcrit_min = ',rhcrit_min
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_rp) 'par_mezcla_rp = ',par_mezcla_rp
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_rp) 'g0_rp = ',g0_rp
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,ES12.4)') 'charnock_max = ',charnock_max
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,ES12.4)') 'charnock_min = ',charnock_min
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_rp) 'lambda_min_rp = ',lambda_min_rp
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_rp) 'cbl_mix_fac_rp = ',cbl_mix_fac_rp
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_rp) 'cs_rp = ',cs_rp
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_rp) 'zhloc_depth_fac_rp = ',zhloc_depth_fac_rp
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_rp) 'ricrit_rp = ',ricrit_rp
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_rp) 'a_ent_1_rp = ',a_ent_1_rp
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_rp) 'a_ent_shr_rp = ',a_ent_shr_rp
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_rp) 'g1_rp = ',g1_rp
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_rp) 'orog_drag_param_rp = ',orog_drag_param_rp
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_lsfc) 'lai_mult_rp = ',lai_mult_rp
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_lsfc) 'lai_mult_rp_max = ',lai_mult_rp_max
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_lsfc) 'lai_mult_rp_min = ',lai_mult_rp_min
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_lsfc) 'z0v_rp_max = ',z0v_rp_max
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_lsfc) 'z0v_rp = ',z0v_rp
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_lsfc) 'z0v_rp_min = ',z0v_rp_min
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_rp) 'z0_soil_rp = ',z0_soil_rp
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_rp) 'z0_urban_mult_rp = ',z0_urban_mult_rp
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_lsfc) 'z0hm_pft_rp_max = ',z0hm_pft_rp_max
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_lsfc) 'z0hm_pft_rp = ',z0hm_pft_rp
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_lsfc) 'z0hm_pft_rp_min = ',z0hm_pft_rp_min
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_lsfc) 'alnir_rp_max = ',alnir_rp_max
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_lsfc) 'alnir_rp = ',alnir_rp
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_lsfc) 'alnir_rp_min = ',alnir_rp_min
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_lsfc) 'alpar_rp_max = ',alpar_rp_max
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_lsfc) 'alpar_rp = ',alpar_rp
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_lsfc) 'alpar_rp_min = ',alpar_rp_min
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_lsfc) 'omega_rp_max = ',omega_rp_max
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_lsfc) 'omega_rp = ',omega_rp
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_lsfc) 'omega_rp_min = ',omega_rp_min
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_lsfc) 'omnir_rp_max = ',omnir_rp_max
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_lsfc) 'omnir_rp = ',omnir_rp
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_lsfc) 'omnir_rp_min = ',omnir_rp_min
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,fmt_rp) 'z0hm_soil_rp = ',z0hm_soil_rp
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,L1)') 'l_skeb2_conv_disp_mod = ',l_skeb2_conv_disp_mod
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,I0)') 'i_pert_theta = ',i_pert_theta
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,I0)') 'i_pert_theta_type = ',i_pert_theta_type
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,ES12.4)') 'decorr_ts_pert_theta = ',decorr_ts_pert_theta
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,ES12.4)') 'mag_pert_theta = ',mag_pert_theta
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,ES12.4)') 'zmin_pert_theta = ',zmin_pert_theta
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,ES12.4)') 'zmax_pert_theta = ',zmax_pert_theta
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,I0)') 'npts_pert_theta = ',npts_pert_theta
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,L1)') 'l_only_pert_near_edge = ',l_only_pert_near_edge
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,I0)') 'npts_pert_from_edge = ',npts_pert_from_edge
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,L1)') 'l_pert_all_points = ',l_pert_all_points
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,L1)') 'l_pert_shape = ',l_pert_shape
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,L1)') 'l_spt = ',l_spt
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,L1)') 'l_spt_rain = ',l_spt_rain
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,L1)') 'l_spt_rad = ',l_spt_rad
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,L1)') 'l_spt_gwd = ',l_spt_gwd
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,L1)') 'l_spt_conv = ',l_spt_conv
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,L1)') 'l_spt_conv_mom = ',l_spt_conv_mom
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,L1)') 'l_spt_cfl = ',l_spt_cfl
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,L1)') 'l_spt_qcons = ',l_spt_qcons
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,L1)') 'l_spt_mse = ',l_spt_mse
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,I0)') 'spt_toplev = ',spt_toplev
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,I0)') 'spt_botlev = ',spt_botlev
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,I0)') 'spt_bot_tap_lev = ',spt_bot_tap_lev
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,I0)') 'spt_top_tap_lev = ',spt_top_tap_lev
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,I0)') 'nsmooth_spt  = ',nsmooth_spt
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,ES12.4)') 'rain_std  = ',rain_std
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,ES12.4)') 'rad_std  = ',rad_std
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,ES12.4)') 'gwd_std  = ',gwd_std
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,ES12.4)') 'conv_std  = ',conv_std
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,ES12.4)') 'sd_orog_thres = ',sd_orog_thres
call umPrint(umMessage,src='stochastic_physics_run_mod')
write(umMessage,'(A,ES12.4)') 'psif_orog_thres = ',psif_orog_thres
call umPrint(umMessage,src='stochastic_physics_run_mod')

call umPrint('- - - - - - end of namelist - - - - - -',                        &
    src='stochastic_physics_run_mod')

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine print_nlist_run_stochastic

end module stochastic_physics_run_mod
