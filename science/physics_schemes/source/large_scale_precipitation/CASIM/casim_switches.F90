! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Cloud Aerosol Interacting Microphysics (CASIM)
! Module to hold the switches for CASIM- temporary until CASIM
! Microphysics is lodged on the trunk.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

module casim_switches

use missing_data_mod, only: imdi

implicit none

! Switches for the activated moments
logical :: l_mp_cloudnumber       = .false.
logical :: l_mp_rainnumber        = .false.
logical :: l_mp_rain3mom          = .false.
logical :: l_mp_icenumber         = .false.
logical :: l_mp_snownumber        = .false.
logical :: l_mp_snow3mom          = .false.
logical :: l_mp_graupnumber       = .false.
logical :: l_mp_graup3mom         = .false.

! Switches for the activated aerosols
logical :: l_mp_activesolliquid   = .false.
logical :: l_mp_activesolrain     = .false.
logical :: l_mp_activeinsolice    = .false.
logical :: l_mp_activesolice      = .false.
logical :: l_mp_activeinsolliquid = .false.
logical :: l_mp_activesolnumber   = .false.
logical :: l_mp_activeinsolnumber = .false.

! Switches to control the aerosol coupling with CASIM
logical :: l_fix_aerosol          = .false.
logical :: l_tracer_aerosol       = .false.
logical :: l_ukca_aerosol         = .false.
logical :: l_ukca_feeding_in      = .false.
logical :: l_ukca_feeding_out     = .false.

logical :: l_asci_progs           = .false.

! Switch for a warm rain only simulation
logical :: l_casim_warm_only      = .false.

! Switch for single moment simulations in the UM. This allows
! high-level UM routines to know that CASIM is running as
! a single-moment simulation and behave accordingly.
logical :: l_kfsm                 = .false.

! Number of extra CASIM prognostics
integer :: n_casim_progs          = imdi

! CASIM moments option
integer :: casim_moments_option   = imdi

! Additional Switches for the cloud scheme from Dan G (Leeds)
logical :: l_cfrac_casim_precip_frac = .false.
logical :: l_mphys_CF_1 = .false.
logical :: l_rad_CF_1 = .false.
logical :: l_div_cfrac_tidy_only=.false.
logical :: l_cfrac_casim_precip_frac_seeding = .false.
logical :: l_cfrac_casim_div_CF = .true.

! Logical switches to maintain bit-comparability
logical :: l_min_cdnc = .true.   !PRF try this - otherwise need to
! qtidy for droplet number and mass before entering radiation.
! This variable controls whether a minimum cloud drop number concentration
! is applied to the fields send to the radiation scheme. If it is, then
! this will likely loose bit-reproducibility with the UM trunk, but will
! ensure that versions of the UM running with the cloud scheme do not crash.

logical :: l_rm_bit_compare = .false.
! A general logical to remove bit-comparability in the UM compared to the
! trunk.


! Switch for a diagnostic cloud scheme (under development)
logical :: l_cfrac_casim_diag_scheme = .false.

! Switch for CASIM diagnostic call
logical :: l_casimmphys_diags = .false.

! determined in casim_ctl based on lam or global model
logical :: l_set_casim_lbc_number = .false.

! CASIM moments - defined in casim_set_dependent_switches and
! for use in the interface
integer :: cloud_mom
integer :: rain_mom
integer :: ice_mom
integer :: snow_mom
integer :: graup_mom

! Parameters for each moment type (single, double, triple)
integer, parameter :: no_moments    = 0
integer, parameter :: single_moment = 1
integer, parameter :: double_moment = 2
integer, parameter :: triple_moment = 3

! Parameters for casim_aerosol_option
integer, parameter :: no_aerosol_modes             = 0
integer, parameter :: soluble_accumulation_coarse  = 1
integer, parameter :: soluble_all_modes            = 2
integer, parameter :: soluble_insoluble_modes      = 3

! CASIM dimensions
integer :: its
integer :: ite
integer :: jts
integer :: jte
integer :: kts
integer :: kte

! First and last dimensions over which CASIM is actually run.
! If l_micro_in_rim = .true., they are equivalent to tdims
! except kle is tdims%k_end - 2
! If l_micro_in_rim = .false., they are the same as the
! rim dimensions below.
integer :: ils
integer :: ile
integer :: jls
integer :: jle
integer :: kls
integer :: kle

!----------------------------------------------------------
! Dimensions without the rim at the start and end.
! So for a LAM these would be tdims%i_start + rimwidth and
! tdims%i_end - rimwidth etc

! Used for CASIM unless l_micro_in_rim = .true.
integer :: irs
integer :: ire
integer :: jrs
integer :: jre
integer :: krs
integer :: kre

!--------------------------------------------------------
! Indexes for aerosol tracers used in the interface
!-------------------------------------------------------
! Soluble aerosol
integer, parameter :: i_AitkenSolMass    = 1
integer, parameter :: i_AitkenSolNumber  = 2
integer, parameter :: i_AccumSolMass     = 3
integer, parameter :: i_AccumSolNumber   = 4
integer, parameter :: i_CoarseSolMass    = 5
integer, parameter :: i_CoarseSolNumber  = 6

! Insoluble aerosol
integer, parameter :: i_CoarseDustMass   = 7
integer, parameter :: i_CoarseDustNumber = 8
integer, parameter :: i_AccumDustMass    = 9
integer, parameter :: i_AccumDustNumber  = 10

integer, parameter :: n_casim_tracers = 10

! Parameters for casim_aerosol_process_level
integer, parameter :: no_processing       = 0
integer, parameter :: passive_processing  = 1
integer, parameter :: full_processing     = 2
integer, parameter :: passive_ice_only    = 3
integer, parameter :: passive_liquid_only = 4

end module casim_switches
