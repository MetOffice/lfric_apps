! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module lsp_autoc_consts_mod

! Description:
! Holds autoconversion constants required by the large-scale
! precipitation scheme
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

use mphys_constants_mod, only: max_drop

use um_types, only: real_umphys

implicit none

! ----------------------------------------------------------------------
!      AUTOCONVERSION TERMS
! ----------------------------------------------------------------------

!     logical, parameter :: L_AUTOCONV_MURK is set in the run_precip namelist.
! Set to .true. to calculate droplet concentration from MURK aerosol,
! which will override L_USE_SULPHATE_AUTOCONV (second indirect effect
! of sulphate aerosol). If both are .false., droplet concentrations
! from mphys_constants_mod are used to be consistent with the values
! used in the radiation scheme.

      ! This next set of parameters is to allow the 3B scheme to
      ! be replicated at 3C/3D
        ! Inhomogeneity factor for autoconversion rate
real(kind=real_umphys),parameter:: inhomog_rate=1.0

      ! Inhomogeneity factor for autoconversion limit
real(kind=real_umphys),parameter:: inhomog_lim=1.0

      ! Threshold droplet radius for autoconversion
real(kind=real_umphys),parameter:: r_thresh=7.0e-6
    ! End of 3B repeated code

    !Do not alter R_AUTO and N_AUTO since these values are effectively
    ! hard wired into a numerical approximation in the autoconversion
    ! code.

    ! Threshold radius for autoconversion
real(kind=real_umphys), parameter :: r_auto=20.0e-6

    ! Critical droplet number for autoconversion
real(kind=real_umphys), parameter :: n_auto=1000.0

    ! The autoconversion powers define the variation of the rate with
    ! liquid water content and droplet concentration. The following are
    ! from Tripoli and Cotton

    !  Dependency of autoconversion rate on droplet concentration
real(kind=real_umphys), parameter :: power_droplet_auto=-0.33333

    ! Dependency of autoconversion rate on water content
real(kind=real_umphys), parameter :: power_qcl_auto=2.33333

    ! Dependency of autoconversion rate on air density
real(kind=real_umphys), parameter :: power_rho_auto=1.33333

    ! CONSTS_AUTO = (4 pi)/( 18 (4 pi/3)^(4/3)) g /  mu (rho_w)^(1/3)
    ! See UM documentation paper 26, equation P26.132

    ! Combination of physical constants
real(kind=real_umphys), parameter :: consts_auto=5907.24

    ! Quantites for calculation of drop number by aerosols.
    ! Need only set if L_AUTOCONV_MURK=.true.  See file C_VISBTY

    ! Scaling concentration (m-3) in aerosol number concentration
    ! Value from Clark et al (2008) QJRMS paper
real(kind=real_umphys), parameter :: n0_clark   = 500.0e6

    ! Scaling mass (kg/kg) in aerosol number calculation from aersol mass
    ! M0_MURK = 4/3 * pi * N0_MURK * (r0)^3 * (rho_aerosol/rho_air)
    ! r0 = 0.11E-6 in Haywood et al (2008) Scheme
    ! rho_aerosol = 1700.0 in Haywood et al (2008) scheme
    ! rho_air = 1.0 in Haywood et al (2008) scheme
    ! Value from Clark et al (2008) QJRMS paper
real(kind=real_umphys), parameter :: m0_clark   = 1.4584e-8

    ! Power in droplet number calculation from aerosols
real(kind=real_umphys), parameter :: power_murk=0.5

    ! Ice water content threshold for graupel autoconversion (kg/m^3)
real(kind=real_umphys), parameter :: auto_graup_qcf_thresh = 3.0e-4

    ! Temperature threshold for graupel autoconversion (degC)
real(kind=real_umphys), parameter :: auto_graup_t_thresh = -4.0

    ! Temperature threshold for graupel autoconversion
real(kind=real_umphys), parameter :: auto_graup_coeff = 0.5

!-----------------------------------------------------------------------
! New parameters added for the droplet tapering process
! Put in this module for now, eventually may pass from the run_precip
! namelist if there is a user requirement for these.
!-----------------------------------------------------------------------

      ! Altitude (m) of where the droplet number goes to its minimum
      ! value
real(kind=real_umphys), parameter :: z_low_nd     = 2000.0

    ! Droplet number at altitudes of z_low_nd and above
real(kind=real_umphys), parameter :: min_drop_alt = 100.0e6

      ! For cosine function, this is equivalent to
      ! ( 375.0E6 m 3 - min_drop_alt ) / 2.0
      ! Where 375.0E6 / m 3 is the maximum droplet number assigned
      ! by the Jones et al (1994) paper.
real(kind=real_umphys)                                                         &
    , parameter :: half_range = (max_drop - min_drop_alt) / 2.0

! Eta value below which to taper the droplet number
real(kind=real_umphys) :: eta_peak

! Eta value above which to use low droplet number
real(kind=real_umphys) :: eta_low_nd

! Eta value at the grid box immediately above the tapering level
real(kind=real_umphys) :: eta_before_taper

integer :: level_peak           ! level of eta_peak
integer :: level_surf           ! level of eta_surf

real(kind=real_umphys) :: vala_fac1, vala_fac2
                                ! constants for drop taper formula

!-----------------------------------------------------------------------
! Parameters for autoconversion scheme if l_warm_new=.true.
! Currently these are for Khairoutdinov & Kogan (2000, MWR)
!-----------------------------------------------------------------------
      ! Pre-factor
real(kind=real_umphys), parameter :: aut_pref = 1350.0
    ! n_drop power
real(kind=real_umphys), parameter :: aut_nc   = -1.79

end module lsp_autoc_consts_mod
