! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module mphys_constants_mod

! Description:
! Holds constants required by the large-scale
! precipitation scheme

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

use missing_data_mod, only: rmdi

use um_types, only: real_umphys

implicit none

save

  ! Offset between ice types (crystals and aggregates)
integer, parameter :: ice_type_offset = 20

  ! Arrays cx and constp hold microphysics constants
  ! Set to missing data to avoid non-initialised data in
  ! the model
real(kind=real_umphys) :: cx(200)     = rmdi
real(kind=real_umphys) :: constp(200) = rmdi

  ! Saved value of m_ci_rp; used to check for changes from
  ! stochastic physics and if we need to rerun code
  ! to calculate microphysics constants
real(kind=real_umphys) :: m_ci_sav = 1.0

  ! air density * specific humidity * bulk velocity
  ! required for Abel and Shipway (diagnostic scheme)
real(kind=real_umphys) :: rho_q_veloc = 3.0e-4 ! This is a typical value

  ! Logical to determine whether to calculate microphysics
  ! constants (only when constants are undefined or input
  ! from stochastic physics change, causing the constants
  ! to need recalculating).
logical :: l_calc_mp_const = .true.

  ! logical to use inhomogeneity parametrization for autoconv
  ! and accretion as described in Boutle et al 2012 QJ
  ! hardwired to true, but left as an option for testing
logical, parameter :: l_inhomog = .true.

!-----------------------------------------------------------------------
! Tolerances for microphysics
!-----------------------------------------------------------------------
real(kind=real_umphys), parameter :: mprog_min = 1.0e-8
! Minimum value of a moist prognostic variable to be considered by the
! microphysics

real(kind=real_umphys), parameter :: mprog_abs = 0.01 * mprog_min
! Absolute value of a moist prognostic variable; below this threshold
! it should really be considered to be zero. Currently set to two orders
! of magnitude below mprog_min

real(kind=real_umphys), parameter :: tnuc = -10.00
! Maximum Temp for ice nuclei nucleation (deg C)
! Temperature at which heterogeneous nucleation of ice first starts to form
! in the model. See UMDP26. Typically minus 10 C.
! Although ice can exist at warmer temperatures, ice is not created until the
! temperature reaches this value. Changing to colder values will reduce amount
! of ice present; warmer values should increase the ice present in model run.

!-----------------------------------------------------------------------
! Parameters for autoconversion scheme if l_warm_new=.true.
! Currently these are for Khairoutdinov & Kogan (2000, MWR)
!-----------------------------------------------------------------------
    ! Pre-factor
real(kind=real_umphys), parameter :: acc_pref = 67.0
  ! qc power
real(kind=real_umphys), parameter :: acc_qc  = 1.15
  ! qr power
real(kind=real_umphys), parameter :: acc_qr  = 1.15

!------------------------------------

! Mass diameter relationship for ice crystals:  m(D) = AI D^BI
! Recommended values for the generic particle size distribution are
! from Brown and Francis and are identical to AI/BI
! AIC = AI = 1.85E-2, BIC = BI = 1.9.  Copy done in module lspcon_mod

!------------------------------------

! Best-Reynolds numbers for aggregates
real(kind=real_umphys), parameter :: lsp_ei      = 0.2072
real(kind=real_umphys), parameter :: lsp_fi      = 0.638

!------------------------------------

! Best-Reynolds numbers for ice crystals
real(kind=real_umphys), parameter :: lsp_eic     = 6.049000e-2
real(kind=real_umphys), parameter :: lsp_fic     = 8.3100e-1

!------------------------------------

! Ice intercept values
real(kind=real_umphys)            :: x1i         = 2.0e6
real(kind=real_umphys)            :: x1ic        = 4.0e7

!------------------------------------
! Rain particle size distribution
! values
real(kind=real_umphys), parameter :: x4r         = 0.00

!------------------------------------
!autoconversion efficiency
real(kind=real_umphys) :: ec_auto     = 0.55
!------------------------------------

!Droplet number over land
real(kind=real_umphys), parameter :: ntot_land   = 3.0e8
!------------------------------------

!Droplet number over sea
real(kind=real_umphys), parameter :: ntot_sea    = 1.0e8
!------------------------------------
! Maximum droplet number assumed through the profile for W-B
real(kind=real_umphys), parameter :: max_drop    = 375.0e6
! Maximum droplet number assumed through the profile for casim for murk
real(kind=real_umphys), parameter :: max_drop_casim = 250.0e6
! Minimum droplet number assumed through the profile for casim for murk
real(kind=real_umphys), parameter :: min_drop_casim = 1.0e3
!------------------------------------

! MURK aerosol relationships
! Scaling concentration (m-3) in aerosol number concentration
real(kind=real_umphys) :: n0_murk

! Scaling mass (kg/kg) in aerosol number calculation from aerosol mass
real(kind=real_umphys) :: m0_murk

! Conversion term from micro-g per kg to kg per kg:
real(kind=real_umphys), parameter :: mu_g_to_kg = 1.0e-9_real_umphys

! Minimum aerosol number
real(kind=real_umphys), parameter :: min_n_aer = 0.0001_real_umphys

! Max aerosol number
real(kind=real_umphys), parameter :: max_n_aer = 1.0e9_real_umphys

!-----------------------------------------------------------------
! Iterations of microphysics
!-----------------------------------------------------------------

! Height at which to start iterative melting

real(kind=real_umphys), parameter :: iter_z = 3380.0

! Corresponds to level 13 in 38 level data set
! (Change added for move to 70 levels)

! Maximum iterations for iterative melting

integer, parameter :: max_it_mlt = 5

! Timestep of each microphysics iteration (s) when iterating over columns
! This is calculated from the number of iterations over columns,
! which in turn can be calculated from a requested sub-step length
! Note, however, that this can differ from the requested sub-step
! length as it needs to divide exactly into the model timestep
real(kind=real_umphys) :: timestep_mp = rmdi

!-----------------------------------------------------------------
! Abel & Shipway (2007) Terms
!-----------------------------------------------------------------

! Specify the default Lambda (function of rain rate)
! This value is given assuming a rain rate of 10 mm per hour.


real(kind=real_umphys) :: lam_evap_enh = 343.9988

! Maximum possible enhancement of Abel & Shipway evaporation
! Largely to avoid issues with it having no upper bound
! and hence evaporating way too too much.

real(kind=real_umphys), parameter :: max_as_enh = 50.0

! Shape dependent riming rate parameters
!
real(kind=real_umphys) :: qclmin_rime = rmdi
real(kind=real_umphys) :: area_ratio_prefac = rmdi
real(kind=real_umphys) :: area_ratio_expn = rmdi

!-----------------------------------------------------------------
! Constants and switches for mixed-phase scheme
!-----------------------------------------------------------------

! Small value used in if tests; taken to be the same as in lsp_moments uses
! for consistency
real(kind=real_umphys), parameter :: mp_smallnum = 2.2e-14

real(kind=real_umphys), parameter :: mp_one_third = 1.0/3.0

end module mphys_constants_mod
