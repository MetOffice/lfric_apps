! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: Module to hold all UKCA variables in RUN_UKCA
!          namelist
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA_UM
!
! Code description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------
module ukca_option_mod

use ukca_config_specification_mod, only:                                       &
  ukca_chem_offline => i_ukca_chem_offline,                                    &
  ukca_chem_offline_be => i_ukca_chem_offline_be,                              &
  ukca_chem_tropisop => i_ukca_chem_tropisop,                                  &
  ukca_chem_strattrop => i_ukca_chem_strattrop,                                &
  ukca_chem_strat => i_ukca_chem_strat,                                        &
  ukca_chem_cristrat => i_ukca_chem_cristrat,                                  &
  ukca_activation_arg => i_ukca_activation_arg,                                &
  ukca_activation_jones => i_ukca_activation_jones,                            &
  ukca_chem_off => i_ukca_chem_off,                                            &
  ukca_lightning_ext => i_light_param_ext
  !!!! Should use 'ukca_api_mod' instead once its indirect dependencies on
  !!!! 'ukca_option_mod' have been removed

use missing_data_mod,      only: rmdi, imdi
use ukca_tracer_stash,     only: a_max_ukcavars

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim
use filenamelength_mod, only: filenamelength
use errormessagelength_mod, only: errormessagelength

implicit none


private :: a_max_ukcavars

! Declarations for UKCA sub-model
! -----------------------------------------------------------------------------

! Namelist items

logical :: l_ukca           =.false. ! True when UKCA is switched on
logical :: l_ukca_aie1      =.false. ! True when 1st aerosol ind effect required
logical :: l_ukca_aie2      =.false. ! True when 2nd aerosol ind effect required

! Main chemistry namelist inputs:
integer :: i_ukca_chem         = 0      ! chemistry scheme to use
logical :: l_ukca_chem_aero    =.false. ! add aerosol precursors to chemistry
logical :: l_ukca_trophet      =.false. ! T for tropospheric heterogeneous chem
logical :: l_ukca_mode         =.false. ! True for UKCA-MODE aerosol scheme
logical :: l_ukca_dust         =.false. ! True for UKCA-mode dust aerosol
logical :: l_ukca_qch4inter    =.false. ! True for interact wetland CH4 ems
logical :: l_ukca_emsdrvn_ch4  =.false. ! T when running UKCA in
                                        ! CH4 emissions-driven mode
logical :: l_ukca_het_psc      =.false. ! True for Het/PSC chemistry
logical :: l_ukca_limit_nat    =.false. ! True for limiting NatPSC formation
                                        ! below specified height
logical :: l_ukca_sa_clim      =.false. ! True to use SPARC surface area density
logical :: l_ukca_h2o_feedback =.false. ! True for H2O feedback from chem
logical :: l_ukca_rado3        =.false. ! T when using UKCA O3 in radiation
logical :: l_ukca_radch4       =.false. ! T when using UKCA CH4 in radiation
logical :: l_ukca_radn2o       =.false. ! T when using UKCA N2O in radiation
logical :: l_ukca_radf11       =.false. ! T when using UKCA CFC-11 in radn
logical :: l_ukca_radf12       =.false. ! T when using UKCA CFC-12 in radn
logical :: l_ukca_radf113      =.false. ! T when using UKCA CFC-113 in radn
logical :: l_ukca_radf22       =.false. ! T when using UKCA HCFC-22 in radn
logical :: l_ukca_radaer       =.false. ! Radiative effects of UKCA aerosols

logical :: l_ukca_radaer_sustrat =.false. ! Use H2SO4 for stratospheric sulphate
logical :: l_ukca_intdd        =.false. ! T when using interact dry deposition
logical :: l_ukca_ddepo3_ocean =.false. ! when T use Luhar et al. (2018)
                                        ! oceanic O3 dry-deposition scheme
logical :: l_ukca_prescribech4 =.false. ! T when prescribing surface ch4
logical :: l_ukca_set_trace_gases =.false. ! T to use UM values for fCO2 etc
logical :: l_ukca_use_background_aerosol =.false. ! use bg aerosol climatology
logical :: l_ukca_dry_dep_so2wet =.false. ! Accounting for wet surfaces in SO2
                                          ! dry deposition.

! Option codes for 'i_ukca_radaer_prescribe_ssa'
! Prescribe SSA in RADAER
!   0: No prescribed single scattering albedo
!   1: Prescribed on 1 wb
!   2: Prescribed on radiation wavebands
! By default set this to 0 to turn off the scheme in RADAER
INTEGER, PARAMETER :: do_not_prescribe = 0
INTEGER, PARAMETER :: prescribe_one_wb = 1
INTEGER, PARAMETER :: prescribe_all_wb = 2
INTEGER :: i_ukca_radaer_prescribe_ssa = 0

! Tuning options for BC absorption
integer :: i_ukca_tune_bc      = imdi
! 0 = No tuning. BC density at default value, standard volume-mixing method
!     used for incorporating BC in the refractive index calculation.
! 1 = BC density at tuned value, standard volume-mixing method still used
!     for incorporating BC in the refractive index calculation.
! 2 = BC density at a different tuned value, Maxwell-Garnet method used for
!     incorporating BC in the refractive index calculation.

! Configuration of heterogeneous chemistry scheme
integer :: i_ukca_hetconfig = 0  ! 0 = default, 1 = JPL-15 recommended coeff.
                                 ! 2 = JPL-15 + br reactions
integer :: i_ukca_topboundary = 1
! 0: Do nothing
! 1: Overwrite top 2 levels with 3rd (except H2O)
! 2: Overwrite top level with 2nd
! 3: Overwrite top level of CO, NO, O3 with ACE-FTS climatology
! 4: Overwrite top level of CO, NO, O3, H2O with ACE-FTS climatology

! T to pass columns to ASAD rather than theta_field
logical :: l_ukca_asad_columns =.false.

! T to pass full domain to ASAD rather than columns or theta_field
logical :: l_ukca_asad_full = .false.

! T to use log(p) to distribute lightning NOx in the vertical
logical :: l_ukca_linox_scaling =.false.
! Select lightning flash frequency parameterisation
integer :: i_ukca_light_param = 1
! 0: No lightning
! 1: Use original Price & Rind parameterisation
! 2: Use updated parameterisation of Luhar et al from CSIRO
! 3: Use external lightning scheme

! T to enable additional print statements to debug asad chemistry solver
logical :: l_ukca_debug_asad = .false.

! T to stop transport of peroxy radicals (StratTrop/CRI only)
logical :: l_ukca_ro2_ntp      = .false.
! T to turn on RO2-permutation chemistry (StratTrop/CRI only)
logical :: l_ukca_ro2_perm     = .false.

! T to use interactive cloud pH routine. F for global pH of 5
logical :: l_ukca_intph = .false.
real :: ph_fit_coeff_a = rmdi ! cloud pH fit parameter a
real :: ph_fit_coeff_b = rmdi ! cloud pH fit parameter b
real :: ph_fit_intercept = rmdi ! cloud pH fit intercept

integer :: chem_timestep = imdi         ! Chemical timestep in seconds for N-R
                                        ! and Offline oxidant schemes
integer :: i_chem_timestep_halvings = imdi ! Integer number of times to half the
                                        ! ASAD chemistry timestep
integer :: dts0 = 300                   ! Default Backward Euler timestep
integer :: nit  = 8                     ! Number of iterations of BE Solver

integer :: nrsteps = imdi

integer :: i_ukca_photol = 0            ! Photolysis scheme to use

integer :: nerupt = imdi                ! Number of explosive eruptions
                                        ! to consider

! Directory pathname for 2d photolysis rates
character (len=filenamelength) :: phot2d_dir  = 'phot2d dir is unset'

integer :: fastjx_numwl = imdi        ! No. of wavelengths to use (8, 12, 18)
integer :: fastjx_mode  = imdi        ! 1 = use just 2D above prescutoff,
                                      ! 2 = merge, 3 = just fastjx)
real :: fastjx_prescutoff = rmdi      ! Press for 2D stratospheric photolysis
character(len=filenamelength) :: jvspec_file ='jvspec file is unset'
                                      ! FastJX spectral file
character(len=filenamelength) :: jvscat_file ='jvscat file is unset'
                                      ! FastJX scatter file
character(len=filenamelength) :: jvsolar_file ='jvsolar file is unset'
                                      ! FastJX solar file
character(len=filenamelength) :: jvspec_dir  ='jvspec dir is unset'
                                      ! Dir for jvspec file

! Use of external photolysis rates in UKCA
logical :: l_environ_jo2 = .false.    ! True if using external O2 -> O(3P) rate
logical :: l_environ_jo2b = .false.   ! True if using external O2 -> O(1D) rate

! Dir for stratospheric aerosol file
character(len=filenamelength) :: dir_strat_aer  = 'dir_strat_aer is unset'
! File for stratospheric aerosol file
character(len=filenamelength) :: file_strat_aer = 'file_strat_aer is unset'

! Filepath for explosive volcanic SO2 emissions
character(len=filenamelength) :: file_volc_so2 = 'file_volc_so2 is unset'

! Switch to choose scheme to use for interactive sea-air exchange of DMS
integer :: i_ukca_dms_flux = imdi
                                 ! 1=LissMerl; 2=Wannin, 3=Nightingale

! Switch to choose scheme to use for interactive sea-salt emissions
integer :: i_primss_method = imdi
                               ! 1=Smith; 2=Gong-Monahan, 3=Combined, 4=Jaegle

! UKCA_MODE control features:
logical :: l_ukca_primsu    =.false. ! T for primary sulphate aerosol emissions
logical :: l_ukca_primss    =.false. ! T for primary sea-salt aerosol emissions
logical :: l_ukca_primbcoc  =.false. ! T for primary BC/OC aerosol emissions
logical :: l_ukca_prim_moc  =.false. ! T for primary marine OC aerosol emissions
logical :: l_ukca_primdu    =.false. ! T for primary dust aerosol emissions
logical :: l_bcoc_ff        =.false. ! T for primary fossil fuel BC/OC emiss.
logical :: l_bcoc_bf        =.false. ! T for primary biofuel BC/OC emissions
logical :: l_bcoc_bm        =.false. ! T for primary biomass BC/OC emissions
logical :: l_ukca_scale_biom_aer_ems = .false. ! Apply scaling factor to
                                     ! biomass burning BC/OC aerosol emissions
logical :: l_ukca_scale_sea_salt_ems = .false. ! Apply scaling factor to
                                               ! sea salt emissions
logical :: l_ukca_scale_marine_pom_ems = .false. ! Apply scaling factor to
                                                 ! Marine POM emissions
logical :: l_ukca_scale_seadms_ems = .false. ! Apply scaling to marine DMS
                                     ! emissions.
logical :: l_mode_bhn_on    =.false.  ! T for binary sulphate nucleation
logical :: l_mode_bln_on    =.false.  ! T for BL sulphate nucleation
integer :: i_ukca_activation_scheme =imdi ! 0 - OFF
!                                         ! 1 - Use AR&G aerosol activation
!                                         ! 2 - Use Jones CDNC
logical :: l_ukca_sfix      =.false. ! T for diagnosing UKCA CCN at
                                     ! fixed supersaturation
! These are switches for the UKCA-iBVOC coupling (in ukca_emission_ctl)
logical :: l_ukca_ibvoc     =.false. ! True for interactive bVOC emissions

! These are switches for the UKCA-INFERNO coupling
logical :: l_ukca_inferno      =.false. ! True for INFERNO fire emissions
logical :: l_ukca_inferno_ch4   =.false. ! True for INFERNO CH4 fire emissions
integer :: i_inferno_emi        = imdi  ! maximum INFERNO emission level

logical :: l_ukca_scale_soa_yield_mt = .false. ! Apply scaling factor to SOA
                                            ! production from monoterpene
logical :: l_ukca_scale_soa_yield_isop = .false. ! Apply scaling factor to SOA
                                            ! production from isoprene

integer :: i_mode_setup     = imdi     ! Defines MODE aerosol scheme
integer :: i_mode_nzts      = imdi     ! No. of substeps for nucleation/
                                       ! sedimentation
integer :: i_mode_bln_param_method = 1 ! 1=activ; 2=kinetc; 3=PNAS/Metzer
                                       ! maps to IBLN in GLOMAP

integer :: i_ukca_nwbins = imdi         ! Controls value of nwbins in Activate
                                        !  See Rosalind West paper for details
                                        !  doi:10.5194/acp-14-6369-2014

! Nitrate emissions scheme control
logical :: l_ukca_fine_no3_prod = .false.
logical :: l_ukca_coarse_no3_prod = .false.
logical :: l_no3_prod_in_aero_step = .false.
real    :: hno3_uptake_coeff = rmdi

! Flag to turn on the Slinn impaction scavenging scheme for dust
logical :: l_dust_mp_slinn_impc_scav = .false.

! Flag to turn on dust ageing (coag, nucl) and activation
logical :: l_dust_mp_ageing = .false.

! Microplastic emissions scheme control
logical :: l_ukca_mp_fragment = .false.
logical :: l_ukca_mp_fibre = .false.

! Not included in namelist at present:
integer, parameter :: i_mode_nucscav = 3 ! Choice of nucl. scavenging co-effs:
                                         ! 1=original, 2=ECHAM5-HAM
                                         ! 3=as(1) but no scav of modes 6&7
real, parameter :: max_z_for_offline_chem = 20000.0
                                         ! Maximum height at which to integrate
                                         ! chemistry with the explicit B-E
                                         ! Offline Oxidants chemistry scheme
logical :: l_ukca_plume_scav = .false.
                                         ! use plume scavenging for aerosol
                                         ! tracers
logical, parameter :: l_ukca_conserve_h = .false.
                                         ! Include hydrogen conservation when
                                         ! l_ukca_h2o_feedback is true
logical :: l_ukca_persist_off = .false.
                                         ! Turn off persistence of spatial
                                         ! arrays

integer :: i_ukca_chem_version = imdi   ! Lowest possible value = 107
! This is the chemical mechanism version identifier, used in ukca_chem_master.
! Set here to the UM version when rates were added to the scheme.
! Any entries can be considered whose version identifier is <= this global
! version number. If there are multiple of those, take the one with the
! largest version number.

real :: mode_parfrac         = rmdi ! Fraction of SO2 emissions as aerosol(%)
real :: mode_aitsol_cvscav   = rmdi ! Plume scavenging fraction for AITSOL
real :: mode_activation_dryr = rmdi ! Activation dry radius in nm
real :: mode_incld_so2_rfrac = rmdi
! fraction of in-cloud oxidised SO2 removed by precipitation
real :: biom_aer_ems_scaling = rmdi ! Biomass-burning aerosol emissions scaling
real :: soa_yield_scaling_mt = rmdi ! Monoterpene SOA yield scaling factor
real :: soa_yield_scaling_isop = rmdi ! Isoprene SOA yield scaling factor

real :: ukca_MeBrMMR         = rmdi ! UKCA trace gas mixing value
real :: ukca_MeClMMR         = rmdi ! UKCA trace gas mixing value
real :: ukca_CH2Br2MMR       = rmdi ! UKCA trace gas mixing value
real :: ukca_H2MMR           = rmdi ! UKCA trace gas mixing value
real :: ukca_N2MMR           = rmdi ! UKCA trace gas mixing value
real :: ukca_CFC115MMR       = rmdi ! UKCA trace gas mixing value
real :: ukca_CCl4MMR         = rmdi ! UKCA trace gas mixing value
real :: ukca_MeCCl3MMR       = rmdi ! UKCA trace gas mixing value
real :: ukca_HCFC141bMMR     = rmdi ! UKCA trace gas mixing value
real :: ukca_HCFC142bMMR     = rmdi ! UKCA trace gas mixing value
real :: ukca_H1211MMR        = rmdi ! UKCA trace gas mixing value
real :: ukca_H1202MMR        = rmdi ! UKCA trace gas mixing value
real :: ukca_H1301MMR        = rmdi ! UKCA trace gas mixing value
real :: ukca_H2402MMR        = rmdi ! UKCA trace gas mixing value
real :: ukca_COSMMR          = rmdi ! UKCA trace gas mixing value

! Variables for new UKCA emission system (based on NetCDF input files)
integer, parameter :: nr_cdf_files      = 100     ! Max nr of NetCDF files
                                                  ! allowed in name list
! path of emiss files
character(len=filenamelength) :: ukca_em_dir = 'ukca_em_dir is unset'
character(len=filenamelength)  :: ukca_em_files(nr_cdf_files) = [              &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset',  &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset',  &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset',  &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset',  &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset',  &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset',  &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset',  &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset',  &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset',  &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset',  &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset',  &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset',  &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset',  &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset',  &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset',  &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset',  &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset',  &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset',  &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset',  &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset',  &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset',  &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset',  &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset',  &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset',  &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset',  &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset',  &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset',  &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset',  &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset',  &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset',  &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset',  &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset',  &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset',  &
'ukca_em_files is unset'                                                       &
]
! Names of emission files

! Information on netCDF files for offline oxidants
integer, parameter :: max_offline_files = 5            ! Max no of NetCDF files
character(len=filenamelength) :: ukca_offline_dir = 'ukca_em_dir is unset'
                                                       ! directory
character(len=filenamelength)  :: ukca_offline_files(max_offline_files) = [    &
 'ukca_offline_files is unset', 'ukca_offline_files is unset',                 &
 'ukca_offline_files is unset', 'ukca_offline_files is unset',                 &
 'ukca_offline_files is unset']
! Names of offline oxidants files

! Control of lower boundary condition scenario
integer ::            i_ukca_scenario = imdi           ! 0=UM; 1=WMOA1; 2=RCP
character(len=filenamelength) :: ukca_RCPdir     = 'ukca_RCPdir is unset'
                                                       ! file location
character(len=filenamelength) :: ukca_RCPfile    = 'ukca_RCPfile is unset'
                                                       ! file name

! Option codes for 'i_ukca_scenario'
integer, parameter, public :: i_ukca_scenario_um = 0
                                ! Take LBC values from UM or namelist
integer, parameter, public :: i_ukca_scenario_wmoa1 = 1
                                ! Use internal UKCA values from WMO A1 for LBCs
integer, parameter, public :: i_ukca_scenario_rcp   = 2
                                ! Take LBC values from RCP file

! UKCA LBC inputs =1 if tr in lbc file
integer ::  tc_lbc_ukca(a_max_ukcavars) = 0

! Options to use with ENDGAME
integer :: i_ukca_conserve_method = 0
                   ! Use separate conservation method for UKCA ?
                   ! 0 = Default tracer conservation
                   ! 1 = Priestley (old) kept for continuity
                   ! 2 = optimised Priestley (recommended default)
                   ! 3 = Conservation Off - only applied to UKCA
                   !      -- was 'conserve_ukca_tracers?'

integer :: i_ukca_hiorder_scheme = imdi
                   ! Use different scheme for High order interpolation?
                   ! Use the same codes as Moisture/tracers
                   ! only active for conserve_meth 1 & 2

logical :: L_ukca_src_in_conservation = .true.
                   ! physics2 sources in conservation ?
                   ! only active for conserve_meth 1 & 2

! Option values for i_ukca_conserve_method
integer, parameter :: ukca_conserve_um  = 0   ! Original/ ADAS scheme
integer, parameter :: priestley_old     = 1
integer, parameter :: priestley_optimal = 2   ! recommended
integer, parameter :: ukca_no_conserve  = 3   ! Do not conserve tracers

logical :: l_ukca_ageair =   .false.   ! Allows user to include the
                                       ! Age-of-air tracer on its own or
                                       ! with any chemistry scheme

! Allows the use of heterogeneous chemistry on aerosol surfaces
! from CLASSIC within UKCA.
logical :: l_ukca_classic_hetchem = .false.

! change resistance based dry deposition scheme to apply deposition
! losses only in level 1
logical :: l_ukca_ddep_lev1 = .false.

! RADAER lookup tables and optical properties namelists.
character(len=filenamelength) :: ukcaaclw = 'ukcaaclw is unset'
                               !  Aitken + Insol acc mode (LW)
character(len=filenamelength) :: ukcaacsw = 'ukcaacsw is unset'
                               !  Aitken + Insol acc mode (SW)
character(len=filenamelength) :: ukcaanlw = 'ukcaanlw is unset'
                               !  Soluble accum mode (LW)
character(len=filenamelength) :: ukcaansw = 'ukcaansw is unset'
                               !  Soluble accum mode (SW)
character(len=filenamelength) :: ukcacrlw = 'ukcacrlw is unset'
                               !  Coarse mode (LW)
character(len=filenamelength) :: ukcacrsw = 'ukcacrsw is unset'
                               !  Coarse mode (SW)
character(len=filenamelength) :: ukcacnlw = 'ukcacnlw is unset'
                               !  Coarse narrow mode (LW)
character(len=filenamelength) :: ukcacnsw = 'ukcacnsw is unset'
                               !  Coarse narrow mode (SW)
character(len=filenamelength) :: ukcasulw = 'ukcasulw is unset'
                               !  Super-coarse mode (LW)
character(len=filenamelength) :: ukcasusw = 'ukcasusw is unset'
                               !  Super-coarse mode (SW)
character(len=filenamelength) :: ukcaprec = 'ukcaprec is unset'
                               !  Precomputed values

character(len=filenamelength) :: ukcasto3 = 'ukcasto3 is unset'
                               ! UKCA Standard temp and O3 file
character(len=filenamelength) :: ukcastrd = 'ukcastrd is unset'
                               ! UKCA Photolysis table

real  :: lightnox_scale_fac = rmdi   ! Lightning NOX ems scale factor

logical :: l_ukca_so2ems_expvolc = .false.  ! If True, SO2 emissions from
                ! specific explosive volcanic eruptions are included for
                ! StratTrop + GLOMAP configuration
logical :: l_ukca_so2ems_plumeria = .false.    ! If True, call plumeria to
                ! calculate plume height of the explosive eruptions
                ! using instantaneous atmospheric conditions
integer :: i_ukca_solcyc = 0  ! Use solar cycle in photolysis
integer :: i_ukca_solcyc_start_year = imdi ! First year of solar cycle data
real :: seadms_ems_scaling = rmdi     ! Marine DMS emission scaling factor
real :: sea_salt_ems_scaling = rmdi   ! Sea salt emission scaling factor
real :: marine_pom_ems_scaling = rmdi ! Marine POM emission scaling factor

! UKCA RADAER prescriptions
! Number of distributions in each spectrum: extinction, absorption
integer, parameter :: n_ukca_radaer = 2

! Information on netCDF files for UKCA RADAER prescriptions
integer, parameter :: n_ukca_radaer_files = n_ukca_radaer*2 ! = 4 netCDF files
character (len=filenamelength) :: ukca_radaer_dir = 'unset'
character (len=filenamelength) :: ukca_radaer_swext_file = 'unset'
character (len=filenamelength) :: ukca_radaer_swabs_file = 'unset'
character (len=filenamelength) :: ukca_radaer_lwext_file = 'unset'
character (len=filenamelength) :: ukca_radaer_lwabs_file = 'unset'

! options for quasi-Newton (Broyden) Method to reduce number of iterations
! in asad_spimpmjp
logical :: l_ukca_quasinewton       = .false.
         ! F=do not perform, T=perform
integer :: i_ukca_quasinewton_start = imdi
         ! iter to start quasi-Newton step (>=2,<=50 2 recommended)
integer :: i_ukca_quasinewton_end   = imdi
         ! iter to stop quasi-Newton step (>=2,<=50 3 recommended)

integer :: i_ukca_sad_months = imdi  ! Used for length of SAD ancil file
integer :: i_ukca_sad_start_year = imdi ! Used to set start year of SAD file

! Options to control how the near-surface values in Age-of-air tracer
! are reset to zero (by level no. or height)
integer :: i_ageair_reset_method = imdi
integer :: max_ageair_reset_level = imdi  ! Max level to which to reset
real    :: max_ageair_reset_height = rmdi ! Max height (m) to which to reset

! Scaling parameters for perturbed parameter ensembles
real :: dry_depvel_so2_scaling = rmdi     ! Scaling factor for SO2 dry
                                          ! deposition velocity
real :: anth_so2_ems_scaling = rmdi       ! Scaling factor for anthropogenic
                                          ! SO2 emissions
real :: dry_depvel_acc_scaling = rmdi     ! Scaling factor for dry deposition
                                          ! velocity for the accumulation mode
real :: acc_cor_scav_scaling = rmdi       ! Scaling factor for scavenging
                                          ! parameters for the accumulation and
                                          ! coarse modes
real :: sigma_updraught_scaling = rmdi    ! Scaling factor for standard
                                          ! deviation of updraught velocities
real :: bc_refrac_im_scaling = rmdi       ! Scaling factor for the imaginary
                                          ! part of the BC refractive index
logical :: l_ukca_scale_ppe = .false.     ! Apply scaling to parameters used in
                                          ! perturbed parameter ensembles

! Define the RUN_UKCA namelist

namelist/run_ukca/ l_ukca, l_ukca_aie1, l_ukca_aie2,                           &
         i_ukca_chem, l_ukca_chem_aero, l_ukca_ageair,                         &
         i_ukca_photol,                                                        &
         l_ukca_mode,                                                          &
         l_ukca_dust,                                                          &
         l_ukca_qch4inter, l_ukca_emsdrvn_ch4,                                 &
         l_ukca_het_psc, l_ukca_sa_clim,                                       &
         l_ukca_h2o_feedback,                                                  &
         l_ukca_rado3, l_ukca_radch4, l_ukca_radn2o,                           &
         l_ukca_radf11, l_ukca_radf12, l_ukca_radf113,                         &
         l_ukca_radf22, l_ukca_radaer, i_ukca_tune_bc,                         &
         l_ukca_radaer_sustrat,                                                &
         i_ukca_radaer_prescribe_ssa,                                          &
         l_ukca_intdd, l_ukca_trophet, l_ukca_prescribech4,                    &
         l_ukca_set_trace_gases, l_ukca_use_background_aerosol,                &
         i_ukca_hetconfig, i_ukca_topboundary,                                 &
         l_ukca_asad_columns, l_ukca_asad_full, l_ukca_ro2_ntp,                &
         l_ukca_ro2_perm, l_ukca_intph, ph_fit_coeff_a, ph_fit_coeff_b,        &
         ph_fit_intercept, l_ukca_primsu, l_ukca_primss, i_primss_method,      &
         l_ukca_fine_no3_prod, l_ukca_coarse_no3_prod, l_no3_prod_in_aero_step,&
         l_ukca_mp_fragment, l_ukca_mp_fibre,                                  &
         hno3_uptake_coeff, l_dust_mp_slinn_impc_scav, l_ukca_primbcoc,        &
         l_ukca_prim_moc, l_ukca_primdu, l_dust_mp_ageing,                     &
         l_bcoc_ff, l_bcoc_bf, l_bcoc_bm, l_mode_bhn_on,                       &
         l_mode_bln_on, i_ukca_activation_scheme,                              &
         l_ukca_sfix, i_mode_setup, i_mode_nzts,                               &
         i_mode_bln_param_method, mode_parfrac,                                &
         mode_aitsol_cvscav, mode_activation_dryr,                             &
         mode_incld_so2_rfrac,                                                 &
         l_ukca_scale_biom_aer_ems, biom_aer_ems_scaling,                      &
         l_ukca_scale_soa_yield_mt, soa_yield_scaling_mt,                      &
         l_ukca_scale_soa_yield_isop, soa_yield_scaling_isop,                  &
         chem_timestep, dts0, nit, nrsteps,                                    &
         jvspec_dir, jvspec_file, jvscat_file, jvsolar_file,                   &
         phot2d_dir, fastjx_numwl, fastjx_mode,                                &
         fastjx_prescutoff, dir_strat_aer, file_strat_aer,                     &
         l_ukca_so2ems_plumeria, file_volc_so2,                                &
         i_ukca_scenario, ukca_RCPdir, ukca_RCPfile,                           &
         ukca_MeBrmmr, ukca_MeClmmr, ukca_CH2Br2mmr, ukca_H2mmr,               &
         ukca_N2mmr, ukca_CFC115mmr, ukca_CCl4mmr,                             &
         ukca_MeCCl3mmr, ukca_HCFC141bmmr, ukca_HCFC142bmmr,                   &
         ukca_H1211mmr, ukca_H1202mmr, ukca_H1301mmr,                          &
         ukca_H2402mmr, ukca_COSmmr,                                           &
         ukca_em_dir, ukca_em_files, tc_lbc_ukca,                              &
         i_ukca_conserve_method, i_ukca_hiorder_scheme,                        &
         L_ukca_src_in_conservation,                                           &
         i_ukca_dms_flux,                                                      &
         ukca_offline_dir, ukca_offline_files, l_ukca_ibvoc,                   &
         l_ukca_classic_hetchem, l_ukca_ddep_lev1,                             &
         ukcaaclw, ukcaacsw, ukcaanlw, ukcaansw, ukcacrlw, ukcacrsw,           &
         ukcacnlw, ukcacnsw, ukcasulw, ukcasusw, ukcaprec, lightnox_scale_fac, &
         L_ukca_so2ems_expvolc, i_ukca_solcyc, i_ukca_solcyc_start_year,       &
         l_ukca_scale_seadms_ems, seadms_ems_scaling,                          &
         l_ukca_quasinewton, i_ukca_quasinewton_start,                         &
         i_ukca_quasinewton_end,                                               &
         i_ageair_reset_method, max_ageair_reset_level,                        &
         max_ageair_reset_height,                                              &
         i_ukca_sad_months, i_ukca_sad_start_year,                             &
         l_ukca_limit_nat, l_ukca_linox_scaling, i_ukca_light_param,           &
         l_ukca_debug_asad, nerupt,                                            &
         l_ukca_inferno, l_ukca_inferno_ch4, i_inferno_emi,                    &
         l_ukca_ddepo3_ocean,                                                  &
         i_ukca_nwbins, i_ukca_chem_version,                                   &
         l_ukca_dry_dep_so2wet, ukca_radaer_dir,                               &
         ukca_radaer_swext_file, ukca_radaer_swabs_file,                       &
         ukca_radaer_lwext_file, ukca_radaer_lwabs_file,                       &
         l_environ_jo2, l_environ_jo2b,                                        &
         l_ukca_scale_sea_salt_ems, sea_salt_ems_scaling,                      &
         l_ukca_scale_marine_pom_ems, marine_pom_ems_scaling,                  &
         dry_depvel_so2_scaling, anth_so2_ems_scaling,                         &
         dry_depvel_acc_scaling, acc_cor_scav_scaling,                         &
         sigma_updraught_scaling, bc_refrac_im_scaling, l_ukca_scale_ppe

! -----------------------------------------------------------------------------
! These are set by UKCA via the 'atmos_ukca_setup' call after the namelist is
! read

integer :: ukca_int_method  = imdi   ! Defines chemical integration method
logical :: l_ukca_chem      =.false. ! True when UKCA chemistry is on
logical :: l_ukca_trop      =.false. ! True for tropospheric chemistry (B-E)
logical :: l_ukca_raq       =.false. ! True for regional air quality chem (B-E)
logical :: l_ukca_raqaero   =.false. ! True for regional air quality chem (B-E)
                                     !         with aerosols
logical :: l_ukca_offline_be=.false. ! True for offline oxidants chem. (B-E)
logical :: l_ukca_tropisop  =.false. ! True for trop chemistry + isoprene
logical :: l_ukca_strat     =.false. ! True for strat+reduced trop chemistry
logical :: l_ukca_strattrop =.false. ! True for std strat+trop chemistry
logical :: l_ukca_cristrat  =.false. ! True for CRI-strat chemistry N-R
logical :: l_ukca_offline   =.false. ! True for offline oxidants chemistry N-R

! These schemes are not yet included but logicals used in code so needed here
logical :: l_ukca_stratcfc  =.false. ! True for extended strat chemistry

! Logical array controlling tracers - set up in primary based
! on calls to tstmsk
logical :: tr_ukca_a (0:a_max_ukcavars) = .false.

! Controls whether UKCA tracers are conserved with same option as
! CLASSIC tracers. Avoids duplication when same (Priestley/ ADAS) scheme
! and hi-order interpolation method are being used.
! This affects whether UKCA tracers are lumped with CLASSIC ones
! in the original or Priestley scheme.
logical :: l_conserve_ukca_with_tr

!DrHook-related parameters
integer(kind=jpim), parameter, private :: zhook_in  = 0
integer(kind=jpim), parameter, private :: zhook_out = 1

character(len=*), parameter, private :: ModuleName='UKCA_OPTION_MOD'

contains

end module ukca_option_mod
