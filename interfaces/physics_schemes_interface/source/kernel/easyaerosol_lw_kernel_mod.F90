!-------------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Populate socrates arrays in the radiation_fields collection
!>        with the easyaerosol longwave optical property variables.
!>        Also, dividing absorption and extinction by air density to convert
!>        to suitable units for socrates (from m-1 to m2 kg-1).

module easyaerosol_lw_kernel_mod
use argument_mod,      only: arg_type,                                         &
                             GH_FIELD, GH_REAL, GH_READ,                       &
                             GH_SCALAR, GH_INTEGER,                            &
                             GH_WRITE, CELL_COLUMN,                            &
                             ANY_DISCONTINUOUS_SPACE_1,                        &
                             ANY_DISCONTINUOUS_SPACE_2,                        &
                             ANY_DISCONTINUOUS_SPACE_3

use kernel_mod,        only: kernel_type

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel.
!> Contains the metadata needed by the Psy layer

type, public, extends(kernel_type) :: easyaerosol_lw_kernel_type
  private
  type(arg_type) :: meta_args(9) = (/                                          &
       ! n_radaer_mode
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                               &
       ! n_lw_band
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                               &
       ! aer_lw_absorption
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),       &
       ! aer_lw_scattering
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),       &
       ! aer_lw_asymmetry
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),       &
       ! easy_absorption_lw
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2),       &
       ! easy_extinction_lw
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2),       &
       ! easy_asymmetry_lw
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2),       &
       ! rho_wtheta
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3)        &
       /)

  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: easyaerosol_lw_code
end type

public :: easyaerosol_lw_code

contains

!> @brief         Easy Aerosol longwave kernel required to populate socrates
!>                arrays with optical property values from the climatology.
!> @param[in]     nlayers             The number of layers
!> @param[in,out] aer_lw_absorption   aerosol LW absorption
!> @param[in,out] aer_lw_scattering   aerosol LW scattering
!> @param[in,out] aer_lw_asymmetry    aerosol LW asymmetry
!> @param[in]     easy_absorption_lw  easyaerosol LW absorption
!> @param[in]     easy_extinction_lw  easyaerosol LW scattering
!> @param[in]     easy_asymmetry_lw   easyaerosol LW asymmetry
!> @param[in]     rho_wtheta          Air density on theta levels
!> @param[in]     ndf_easy_lw         No. of DOFs per cell for aer_lw space
!> @param[in]     undf_easy_lw        No. unique of DOFs for aer_lw space
!> @param[in]     map_easy_lw         Dofmap for aer_lw space column base cell
!> @param[in]     ndf_aer_lw          No. of DOFs per cell for aer_lw space
!> @param[in]     undf_aer_lw         No. unique of DOFs for aer_lw space
!> @param[in]     map_aer_lw          Dofmap for aer_lw space column base cell
!> @param[in]     ndf_rho             No. DOFs per cell for 2D space
!> @param[in]     undf_rho            No. unique DOFs for rho_wtheta space
!> @param[in]     map_rho             Dofmap for rho_wtheta space

subroutine easyaerosol_lw_code( nlayers,                                       &
                                n_radaer_mode,                                 &
                                n_lw_band,                                     &
                                aer_lw_absorption,                             &
                                aer_lw_scattering,                             &
                                aer_lw_asymmetry,                              &
                                easy_absorption_lw,                            &
                                easy_extinction_lw,                            &
                                easy_asymmetry_lw,                             &
                                rho_wtheta,                                    &
                                ndf_aer_lw, undf_aer_lw, map_aer_lw,           &
                                ndf_easy_lw, undf_easy_lw, map_easy_lw,        &
                                ndf_rho,  undf_rho,  map_rho )

  use constants_mod,                     only: r_def, i_def, i_um

  implicit none

  ! Arguments

  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: n_radaer_mode
  integer(kind=i_def), intent(in) :: n_lw_band

  integer(kind=i_def), intent(in) :: ndf_aer_lw
  integer(kind=i_def), intent(in) :: undf_aer_lw
  integer(kind=i_def), dimension(ndf_aer_lw), intent(in) :: map_aer_lw

  integer(kind=i_def), intent(in) :: ndf_easy_lw
  integer(kind=i_def), intent(in) :: undf_easy_lw
  integer(kind=i_def), dimension(ndf_easy_lw), intent(in) :: map_easy_lw

  integer(kind=i_def), intent(in) :: ndf_rho
  integer(kind=i_def), intent(in) :: undf_rho
  integer(kind=i_def), dimension(ndf_rho), intent(in) :: map_rho

  real(kind=r_def), intent(inout), dimension(undf_aer_lw)  :: aer_lw_absorption
  real(kind=r_def), intent(inout), dimension(undf_aer_lw)  :: aer_lw_scattering
  real(kind=r_def), intent(inout), dimension(undf_aer_lw)  :: aer_lw_asymmetry
  real(kind=r_def), intent(in),    dimension(undf_easy_lw) :: easy_absorption_lw
  real(kind=r_def), intent(in),    dimension(undf_easy_lw) :: easy_extinction_lw
  real(kind=r_def), intent(in),    dimension(undf_easy_lw) :: easy_asymmetry_lw
  real(kind=r_def), intent(in),    dimension(undf_rho)     :: rho_wtheta

  ! Loop counters / array index
  integer(i_um) :: k, i_band, i_bnd_mode, i_easy, i_aer, i_rho

  !-----------------------------------------------------------------------
  ! LW aerosol optical properties are copied from
  ! the aerosol_fields arrays to radiation_fields arrays.
  !
  i_bnd_mode = 0
  do i_band = 1, n_lw_band

    ! The socrates arrays will be populated with radaer outputs
    ! for modes 1 to n_radaer_mode so need to advance i_bnd_mode by
    ! n_radaer_mode + 1 each time this loop increments i_band
    !
    i_bnd_mode = i_bnd_mode + n_radaer_mode + 1

    ! The first model level is the surface where aerosols are zero
    k=0
    i_aer = map_aer_lw (1) + ( (i_bnd_mode-1) * (nlayers+1) ) + k
    aer_lw_absorption (i_aer) = 0.0_r_def
    aer_lw_scattering (i_aer) = 0.0_r_def
    aer_lw_asymmetry  (i_aer) = 0.0_r_def

    do k = 1, nlayers

      ! The easyaerosol arrays are arranged differently with only one mode
      ! and with the data ordered so that the waveband advances from one
      ! element to the next and data advances through the model levels in
      ! strides of n_lw_band

      i_easy = map_easy_lw(1) + (k * n_lw_band) + i_band-1
      i_aer  = map_aer_lw(1)  + ((i_bnd_mode-1) * (nlayers+1)) + k
      i_rho  = map_rho(1)     + k

      ! There is no need to divide absorption and scattering by aerosol
      ! mixing ratio here as mmr is unknown and was set to dummy values of
      ! 1.0 in easyaerosol_aer_kernel_mod.
      !
      aer_lw_absorption (i_aer) =   easy_absorption_lw(i_easy)                 &
                                    / rho_wtheta(i_rho)
      aer_lw_scattering (i_aer) =  (easy_extinction_lw(i_easy) -               &
                                    easy_absorption_lw(i_easy))                &
                                    / rho_wtheta(i_rho)
      aer_lw_asymmetry  (i_aer) =  easy_asymmetry_lw (i_easy)

    end do
  end do

end subroutine easyaerosol_lw_code
end module easyaerosol_lw_kernel_mod
