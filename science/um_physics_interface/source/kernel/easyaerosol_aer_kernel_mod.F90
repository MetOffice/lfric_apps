!-------------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Populate the socrates aer_mix_ratio arrays in the radiation_fields
!> collection with dummy values of 1.0 for the elements corresponding to
!> the easyaerosol mode.
!>
!> Socrates assumes that the aerosol mixing ratio is known for every aerosol
!> mode and will multiply absorption and scattering coefficients by the
!> aer_mix_ratio array. However, the optical properties read in from the ancils
!> are effectively already the product of the aerosol absorption or extinction
!> efficiency multiplied by aer_mix_ratio so don't need this multiplication step.
!>
!> Rather than altering socrates to exclude the multipication step for
!> the special case when optical properties originate from easyaerosol, we can
!> set the easyaerosol aer_mix_ratio values to 1.0 here.
!>
!> A separate kernel has been created to do this since only the last seventh of
!> the elements in the column corresponds to easyaerosol and a loop is written to
!> locate these.


module easyaerosol_aer_kernel_mod
use argument_mod,      only: arg_type,                                         &
                             GH_FIELD, GH_REAL,                                &
                             GH_SCALAR, GH_INTEGER, GH_READ,                   &
                             GH_WRITE, CELL_COLUMN,                            &
                             ANY_DISCONTINUOUS_SPACE_1

use kernel_mod,        only: kernel_type

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel.
!> Contains the metadata needed by the Psy layer

type, public, extends(kernel_type) :: easyaerosol_aer_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                                          &
      ! mode_dimen
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                               &
      ! aer_mix_ratio
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1) /)

  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: easyaerosol_aer_code
end type

public :: easyaerosol_aer_code

contains

!> @brief  Locate aer_mix_ratio values associated with easyaerosol and set to 1.0
!> @param[in]     nlayers             The number of layers
!> @param[in]     mode_dimen          The number of modes (RADAER + EasyAerosol)
!> @param[in,out] aer_mix_ratio       aerosol mixing ratios
!> @param[in]     ndf_mode            No. of DOFs per cell for mode space
!> @param[in]     undf_mode           No. unique of DOFs for mode space
!> @param[in]     map_mode            Dofmap for mode space column base cell

subroutine easyaerosol_aer_code( nlayers,                                      &
                                 mode_dimen,                                   &
                                 aer_mix_ratio,                                &
                                 ndf_mode, undf_mode, map_mode)

  use constants_mod,       only: r_def, i_def, i_um

  implicit none

  ! Arguments

  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: mode_dimen
  integer(kind=i_def), intent(in) :: ndf_mode
  integer(kind=i_def), intent(in) :: undf_mode
  integer(kind=i_def), dimension(ndf_mode), intent(in) :: map_mode

  real(kind=r_def), intent(inout), dimension(undf_mode) :: aer_mix_ratio

  ! Loop counters / array index
  integer(i_um) :: k, i_easy

  !-----------------------------------------------------------------------

  ! Easyaerosol data will occupy the final mode of the socrates array.
  ! The total number of modes is mode_dimen = n_radaer_modes + 1

  ! Dummy aerosol mixing ratios for easyaerosols
  ! Set to 0.0 for surface layer and 1.0 for all other levels
  !
  k = 0
  i_easy = map_mode(1) + ((mode_dimen-1) * (nlayers+1)) + k
  aer_mix_ratio(i_easy) = 0.0_r_def

  do k = 1, nlayers
    i_easy = map_mode(1) + ((mode_dimen-1) * (nlayers+1)) + k
    aer_mix_ratio(i_easy) = 1.0_r_def
  end do

end subroutine easyaerosol_aer_code

end module easyaerosol_aer_kernel_mod
