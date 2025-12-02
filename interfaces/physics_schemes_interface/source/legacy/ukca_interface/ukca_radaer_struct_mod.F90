! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose:
! Defines maximum dimensions
! Defines type ukca_radaer_struct, the structure used by UKCA_RADAER
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA_UM
!
! Contained subroutines:
!      allocate_radaer_struct
!      allocate_ukca_radaer_fields
!      deallocate_ukca_radaer_fields
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------

module ukca_radaer_struct_mod

implicit none

character(len=*), parameter, private :: ModuleName='UKCA_RADAER_STRUCT_MOD'

! Internal IDs for the types of UKCA aerosol components.
!
! Note: When adding an aerosol component, also increase
!       npd_ukca_maxcomptype in module ukca_radaer_precalc_mod and
!       update the pre-computed file read by ukca_radaer_read_precalc.
!       Any new aerosol component should be added before ip_ukca_water.
!
! Water is not an aerosol species, but is included here
! as it behaves like one in ukca_radaer_band_average().
! However, no aerosol component should be of type ip_ukca_water.

! cp_su=1  from ukca_mode_setup
! cp_bc=2  from ukca_mode_setup
! cp_oc=3  from ukca_mode_setup
! cp_cl=4  from ukca_mode_setup
! cp_du=5  from ukca_mode_setup
! cp_so=6  from ukca_mode_setup
! Note that cp_no3=7 from ukca_mode_setup is not being used by RADAER
! Note that cp_nh4=8 from ukca_mode_setup is not being used by RADAER
integer, parameter :: ip_ukca_h2so4 =  9
! cp_mp=10 from ukca_mode_setup
integer, parameter :: ip_ukca_water = 11

integer, save :: npd_ukca_cpnt       ! nmodes*ncp

integer, save :: ncp_max_x_nmodes    ! nmodes * ncp_max

! Thresholds on the modal mass-mixing ratio and modal number
! concentrations above which aerosol optical properties are to be
! computed. Placed here as used by multiple subroutines.
real, parameter :: threshold_mmr = 1.0e-12 ! kg/kg
! Corresponds to burden of 0.01 mg/m2 if mmr=1.e-12 everywhere.

real, parameter :: threshold_vol = 1.0e-25 ! m3/m3
! Corresponds to particle diameter < 10nm

real, parameter :: threshold_nbr = 1.0e+00 ! m-3
! A single coarse-mode particle with d=10um per would give
! a mixing ratio of ~1.e-12kg/kg in the lower troposphere

! Main structure holding all the variables needed for
! interacting UKCA aerosols with radiation.

type :: ukca_radaer_struct

  !
  ! Information about UKCA aerosol modes
  !

  ! Actual number of modes, with a default value for
  ! minimising array dimensions.
  integer :: n_mode = 1

  ! Type of mode (i.e. nucleation, Aitken, accum, or coarse)
  integer, allocatable :: i_mode_type(:)

  ! Solubility of mode (soluble if true)
  logical, allocatable :: l_soluble(:)

  ! Lower and upper limits on the geometric mean diameter (m)
  ! in each mode
  real, allocatable :: d0low(:)
  real, allocatable :: d0up(:)

  ! Geometric standard deviation in this mode
  real, allocatable :: sigma(:)

  ! STASH code, D1 address, number of levels and total length
  ! of the D1 fields corresponding to modal dry diameter.
  integer, allocatable :: stashcode_dry(:)
  integer, allocatable :: d1_address_dry(:)
  integer, allocatable :: d1_nlevs_dry(:)
  integer, allocatable :: d1_length_dry(:)

  ! STASH code, D1 address, number of levels and total length
  ! of the D1 fields corresponding to modal wet diameter.
  integer, allocatable :: stashcode_wet(:)
  integer, allocatable :: d1_address_wet(:)
  integer, allocatable :: d1_nlevs_wet(:)
  integer, allocatable :: d1_length_wet(:)

  ! STASH code, D1 address, number of levels and total length
  ! of the D1 fields corresponding to modal density.
  integer, allocatable :: stashcode_rho(:)
  integer, allocatable :: d1_address_rho(:)
  integer, allocatable :: d1_nlevs_rho(:)
  integer, allocatable :: d1_length_rho(:)

  ! STASH code, D1 address, number of levels and total length
  ! of the D1 fields corresponding to water volume in each mode.
  integer, allocatable :: stashcode_wtv(:)
  integer, allocatable :: d1_address_wtv(:)
  integer, allocatable :: d1_nlevs_wtv(:)
  integer, allocatable :: d1_length_wtv(:)

  ! STASH code, D1 address, number of levels, total length
  ! and halo type of the D1 fields corresponding to
  ! modal number concentrations.
  integer, allocatable :: stashcode_nbr(:)
  integer, allocatable :: d1_address_nbr(:)
  integer, allocatable :: d1_nlevs_nbr(:)
  integer, allocatable :: d1_length_nbr(:)
  integer, allocatable :: d1_halo_type_nbr(:)

  ! Number of components in each mode and index of each
  ! component in array ukca_cpnt_info
  integer, allocatable :: n_cpnt_in_mode(:)
  integer, allocatable :: i_cpnt_index(:,:)

  ! Modal diameter of the dry aerosol (m)
  real,    allocatable :: dry_diam(:, :, :, :)

  ! Modal diameter of the wet aerosol (m)
  real,    allocatable :: wet_diam(:, :, :, :)

  ! Modal densities (kg/m3)
  real,    allocatable :: modal_rho(:, :, :, :)

  ! Modal volumes (including water for soluble modes)
  real,    allocatable :: modal_vol(:, :, :, :)

  ! Fractional volume of water in each mode
  real,    allocatable :: modal_wtv(:, :, :, :)

  ! Modal number concentrations
  real,    allocatable :: modal_nbr(:, :, :, :)

  !
  ! Information about UKCA aerosol components
  !

  ! Actual number of components, with a default value for
  ! minimising array dimensions.
  integer :: n_cpnt = 1

  ! Size of ncp_max
  integer :: ncp_max

  ! Size of ncp
  integer :: ncp

  ! Size of nmodes
  integer :: nmodes

  ! Size of ncp_max_x_nmodes
  integer :: ncp_max_x_nmodes

  ! Type of component (e.g. sulphate)
  integer, allocatable :: i_cpnt_type(:)

  ! Mass density of each component (kg/m3)
  real,    allocatable :: density(:)

  ! Array index of the mode this component belongs to
  integer, allocatable :: i_mode(:)

  ! Array index for translating component and mode
  ! to ukca_radaer index
  integer, allocatable :: idx_cpnt_mode(:,:)

  ! STASH code, D1 address, number of levels, total length
  ! and halo type of the D1 fields corresponding to the
  ! mass-mixing ratio of each component.
  integer, allocatable :: stashcode_mmr(:)
  integer, allocatable :: d1_address_mmr(:)
  integer, allocatable :: d1_nlevs_mmr(:)
  integer, allocatable :: d1_length_mmr(:)
  integer, allocatable :: d1_halo_type_mmr(:)

  ! STASH code, D1 address, number of levels and total length
  ! of the D1 fields corresponding to the fractional volume of
  ! each component.
  integer, allocatable :: stashcode_cvl(:)
  integer, allocatable :: d1_address_cvl(:)
  integer, allocatable :: d1_nlevs_cvl(:)
  integer, allocatable :: d1_length_cvl(:)

  ! Component mass-mixing ratio (kg/kg)
  real,    allocatable :: mix_ratio(:, :, :, :)

  ! Component volumes
  real,    allocatable :: comp_vol(:, :, :, :)

  ! Switch: if true, use sulphuric acid optical properties in the
  ! stratosphere, instead of ammonium sulphate (if false).
  ! Has the same default as run_ukca:l_ukca_radaer_sustrat or
  ! run_glomap_aeroclim:l_glomap_clim_radaer_sustrat
  ! from which is it assigned.
  logical :: l_sustrat = .false.

  ! Switch: if true, nitrate scheme is on and sodium nitrate
  ! refractive indices are included in the pcalc file
  logical :: l_nitrate = .false.

  ! Switch: if true, use the narrow coarse mode LUTs for the
  ! coarse insoluble mode. This is the case if the super coarse
  ! insoluble mode is used
  logical :: l_cornarrow_ins = .false.

end type ukca_radaer_struct

public :: allocate_radaer_struct,                                              &
          allocate_ukca_radaer_fields,                                         &
          deallocate_ukca_radaer_fields

contains

! #############################################################################

subroutine allocate_radaer_struct(ukca_radaer, glomap_variables)

! To allocate arrays in the structure ukca_radaer using the number of modes
! and number of components configured in the GLOMAP setup routine.

use ukca_mode_setup, only: nmodes, ncp_max, cp_no3, glomap_variables_type
use parkind1,        only: jprb, jpim
use yomhook,         only: lhook, dr_hook
implicit none

type(ukca_radaer_struct),    intent(in out) :: ukca_radaer
type(glomap_variables_type), intent(in)     :: glomap_variables

integer :: ncp

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='ALLOCATE_RADAER_STRUCT'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncp                          = glomap_variables%ncp
npd_ukca_cpnt                = ncp * nmodes
ncp_max_x_nmodes             = ncp_max * nmodes

ukca_radaer%ncp              = ncp
ukca_radaer%ncp_max          = ncp_max
ukca_radaer%nmodes           = nmodes
ukca_radaer%ncp_max_x_nmodes = ncp_max_x_nmodes

! Amend the indices of H2SO4 and water for sodium nitrate inclusion
if (ncp >= cp_no3) then
  if (any(glomap_variables%component (:,cp_no3))) then
    ukca_radaer%l_nitrate = .true.
  end if
end if

if (.not. allocated(ukca_radaer%i_mode_type))                                  &
      allocate(ukca_radaer%i_mode_type(nmodes))
if (.not. allocated(ukca_radaer%l_soluble))                                    &
      allocate(ukca_radaer%l_soluble(nmodes))
if (.not. allocated(ukca_radaer%d0low))                                        &
      allocate(ukca_radaer%d0low(nmodes))
if (.not. allocated(ukca_radaer%d0up))                                         &
      allocate(ukca_radaer%d0up(nmodes))
if (.not. allocated(ukca_radaer%sigma))                                        &
      allocate(ukca_radaer%sigma(nmodes))
if (.not. allocated(ukca_radaer%stashcode_dry))                                &
      allocate(ukca_radaer%stashcode_dry(nmodes))
if (.not. allocated(ukca_radaer%d1_address_dry))                               &
      allocate(ukca_radaer%d1_address_dry(nmodes))
if (.not. allocated(ukca_radaer%d1_nlevs_dry))                                 &
      allocate(ukca_radaer%d1_nlevs_dry(nmodes))
if (.not. allocated(ukca_radaer%d1_length_dry))                                &
      allocate(ukca_radaer%d1_length_dry(nmodes))
if (.not. allocated(ukca_radaer%stashcode_wet))                                &
      allocate(ukca_radaer%stashcode_wet(nmodes))
if (.not. allocated(ukca_radaer%d1_address_wet))                               &
      allocate(ukca_radaer%d1_address_wet(nmodes))
if (.not. allocated(ukca_radaer%d1_nlevs_wet))                                 &
      allocate(ukca_radaer%d1_nlevs_wet(nmodes))
if (.not. allocated(ukca_radaer%d1_length_wet))                                &
      allocate(ukca_radaer%d1_length_wet(nmodes))
if (.not. allocated(ukca_radaer%stashcode_rho))                                &
      allocate(ukca_radaer%stashcode_rho(nmodes))
if (.not. allocated(ukca_radaer%d1_address_rho))                               &
      allocate(ukca_radaer%d1_address_rho(nmodes))
if (.not. allocated(ukca_radaer%d1_nlevs_rho))                                 &
      allocate(ukca_radaer%d1_nlevs_rho(nmodes))
if (.not. allocated(ukca_radaer%d1_length_rho))                                &
      allocate(ukca_radaer%d1_length_rho(nmodes))
if (.not. allocated(ukca_radaer%d1_length_wet))                                &
      allocate(ukca_radaer%d1_length_wet(nmodes))
if (.not. allocated(ukca_radaer%stashcode_wtv))                                &
      allocate(ukca_radaer%stashcode_wtv(nmodes))
if (.not. allocated(ukca_radaer%d1_address_wtv))                               &
      allocate(ukca_radaer%d1_address_wtv(nmodes))
if (.not. allocated(ukca_radaer%d1_nlevs_wtv))                                 &
      allocate(ukca_radaer%d1_nlevs_wtv(nmodes))
if (.not. allocated(ukca_radaer%d1_length_wtv))                                &
      allocate(ukca_radaer%d1_length_wtv(nmodes))
if (.not. allocated(ukca_radaer%stashcode_nbr))                                &
      allocate(ukca_radaer%stashcode_nbr(nmodes))
if (.not. allocated(ukca_radaer%d1_address_nbr))                               &
      allocate(ukca_radaer%d1_address_nbr(nmodes))
if (.not. allocated(ukca_radaer%d1_nlevs_nbr))                                 &
      allocate(ukca_radaer%d1_nlevs_nbr(nmodes))
if (.not. allocated(ukca_radaer%d1_length_nbr))                                &
      allocate(ukca_radaer%d1_length_nbr(nmodes))
if (.not. allocated(ukca_radaer%d1_halo_type_nbr))                             &
      allocate(ukca_radaer%d1_halo_type_nbr(nmodes))
if (.not. allocated(ukca_radaer%n_cpnt_in_mode))                               &
      allocate(ukca_radaer%n_cpnt_in_mode(nmodes))
if (.not. allocated(ukca_radaer%i_cpnt_index))                                 &
      allocate(ukca_radaer%i_cpnt_index(ncp_max,nmodes))
if (.not. allocated(ukca_radaer%i_cpnt_type))                                  &
      allocate(ukca_radaer%i_cpnt_type(ncp_max_x_nmodes))
if (.not. allocated(ukca_radaer%density))                                      &
      allocate(ukca_radaer%density(npd_ukca_cpnt))
if (.not. allocated(ukca_radaer%i_mode))                                       &
      allocate(ukca_radaer%i_mode(npd_ukca_cpnt))
if (.not. allocated(ukca_radaer%idx_cpnt_mode))                                &
      allocate(ukca_radaer%idx_cpnt_mode(nmodes,ncp_max))
if (.not. allocated(ukca_radaer%stashcode_mmr))                                &
      allocate(ukca_radaer%stashcode_mmr(npd_ukca_cpnt))
if (.not. allocated(ukca_radaer%d1_address_mmr))                               &
      allocate(ukca_radaer%d1_address_mmr(npd_ukca_cpnt))
if (.not. allocated(ukca_radaer%d1_nlevs_mmr))                                 &
      allocate(ukca_radaer%d1_nlevs_mmr(npd_ukca_cpnt))
if (.not. allocated(ukca_radaer%d1_length_mmr))                                &
      allocate(ukca_radaer%d1_length_mmr(npd_ukca_cpnt))
if (.not. allocated(ukca_radaer%d1_halo_type_mmr))                             &
      allocate(ukca_radaer%d1_halo_type_mmr(npd_ukca_cpnt))
if (.not. allocated(ukca_radaer%stashcode_cvl))                                &
      allocate(ukca_radaer%stashcode_cvl(npd_ukca_cpnt))
if (.not. allocated(ukca_radaer%d1_address_cvl))                               &
      allocate(ukca_radaer%d1_address_cvl(npd_ukca_cpnt))
if (.not. allocated(ukca_radaer%d1_nlevs_cvl))                                 &
      allocate(ukca_radaer%d1_nlevs_cvl(npd_ukca_cpnt))
if (.not. allocated(ukca_radaer%d1_length_cvl))                                &
      allocate(ukca_radaer%d1_length_cvl(npd_ukca_cpnt))

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine allocate_radaer_struct

! =============================================================================
! Allocate ukca_radaer fields
!  Full ukca_radaer fields are allocated on radiation time steps
!  Minimum ukca_radaer fields are allocated on other time steps
! =============================================================================
subroutine allocate_ukca_radaer_fields( row_length_in,                         &
                                        rows_in,                               &
                                        model_levels_in,                       &
                                        n_cpnt_in,                             &
                                        n_mode_in,                             &
                                        ukca_radaer )

use parkind1,             only:                                                &
    jpim,                                                                      &
    jprb

use yomhook,              only:                                                &
    lhook,                                                                     &
    dr_hook

implicit none

! Arguments

integer, intent(in) :: row_length_in    ! Either row_length or 1
integer, intent(in) :: rows_in          ! Either rows or 1
integer, intent(in) :: model_levels_in  ! Either model_levels or 1
integer, intent(in) :: n_cpnt_in        ! Either ukca_radaer%n_cpnt or 1
integer, intent(in) :: n_mode_in        ! Either ukca_radaer%n_mode or 1

type(ukca_radaer_struct), intent(in out) :: ukca_radaer

! Local variables

character(len=*),   parameter :: RoutineName='ALLOCATE_UKCA_RADAER_FIELDS'
integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

if (lhook) call dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

if (.not. allocated( ukca_radaer%mix_ratio ) )                                 &
          allocate(  ukca_radaer%mix_ratio ( row_length_in,                    &
                                             rows_in,                          &
                                             model_levels_in,                  &
                                             n_cpnt_in ) )

if (.not. allocated( ukca_radaer%comp_vol ) )                                  &
          allocate(  ukca_radaer%comp_vol (  row_length_in,                    &
                                             rows_in,                          &
                                             model_levels_in,                  &
                                             n_cpnt_in ) )

if (.not. allocated( ukca_radaer%dry_diam ) )                                  &
          allocate(  ukca_radaer%dry_diam (  row_length_in,                    &
                                             rows_in,                          &
                                             model_levels_in,                  &
                                             n_mode_in ) )

if (.not. allocated( ukca_radaer%wet_diam ) )                                  &
          allocate(  ukca_radaer%wet_diam (  row_length_in,                    &
                                             rows_in,                          &
                                             model_levels_in,                  &
                                             n_mode_in ) )

if (.not. allocated( ukca_radaer%modal_rho ) )                                 &
          allocate(  ukca_radaer%modal_rho ( row_length_in,                    &
                                             rows_in,                          &
                                             model_levels_in,                  &
                                             n_mode_in ) )

if (.not. allocated( ukca_radaer%modal_wtv ) )                                 &
          allocate(  ukca_radaer%modal_wtv ( row_length_in,                    &
                                             rows_in,                          &
                                             model_levels_in,                  &
                                             n_mode_in ) )

if (.not. allocated( ukca_radaer%modal_vol ) )                                 &
          allocate(  ukca_radaer%modal_vol ( row_length_in,                    &
                                             rows_in,                          &
                                             model_levels_in,                  &
                                             n_mode_in ) )

if (.not. allocated( ukca_radaer%modal_nbr ) )                                 &
          allocate(  ukca_radaer%modal_nbr ( row_length_in,                    &
                                             rows_in,                          &
                                             model_levels_in,                  &
                                             n_mode_in ) )

if (lhook) call dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
return
end subroutine allocate_ukca_radaer_fields

! =============================================================================
! Deallocate ukca_radaer fields
! =============================================================================
subroutine deallocate_ukca_radaer_fields( ukca_radaer )

use parkind1, only:                                                            &
    jpim,                                                                      &
    jprb

use yomhook,  only:                                                            &
    lhook,                                                                     &
    dr_hook

implicit none

! Arguments

type(ukca_radaer_struct), intent(in out) :: ukca_radaer

! Local variables

character(len=*),   parameter :: RoutineName='DEALLOCATE_UKCA_RADAER_FIELDS'
integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

if (lhook) call dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

if (allocated( ukca_radaer%modal_nbr ))  deallocate( ukca_radaer%modal_nbr )
if (allocated( ukca_radaer%modal_vol ))  deallocate( ukca_radaer%modal_vol )
if (allocated( ukca_radaer%modal_wtv ))  deallocate( ukca_radaer%modal_wtv )
if (allocated( ukca_radaer%modal_rho ))  deallocate( ukca_radaer%modal_rho )
if (allocated( ukca_radaer%wet_diam  ))  deallocate( ukca_radaer%wet_diam  )
if (allocated( ukca_radaer%dry_diam  ))  deallocate( ukca_radaer%dry_diam  )
if (allocated( ukca_radaer%comp_vol  ))  deallocate( ukca_radaer%comp_vol  )
if (allocated( ukca_radaer%mix_ratio ))  deallocate( ukca_radaer%mix_ratio )

if (lhook) call dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
end subroutine deallocate_ukca_radaer_fields

end module ukca_radaer_struct_mod
