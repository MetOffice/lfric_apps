!-------------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Calculate fields for easyaerosol diagnostics

module easyaerosol_fld2diag_kernel_mod
use argument_mod,      only: arg_type,                                         &
                             GH_FIELD, GH_REAL, GH_READ,                       &
                             GH_INTEGER, GH_SCALAR,                            &
                             GH_WRITE, CELL_COLUMN,                            &
                             ANY_DISCONTINUOUS_SPACE_1,                        &
                             ANY_DISCONTINUOUS_SPACE_2

use kernel_mod,        only: kernel_type

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel.
!> Contains the metadata needed by the Psy layer

type, public, extends(kernel_type) :: easyaerosol_fld2diag_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                                          &
       ! n_band
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                               &
       ! aer_prop
       arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1),        &
       ! aer_prop_diag
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_2)        &
       /)

  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: easyaerosol_fld2diag_code
end type

public :: easyaerosol_fld2diag_code

contains

!> @brief         Transfer data from easyaerosol fields into diagnostic arrays.
!>                This is not totally straightforward as the diag arrays have
!>                a different shape to those containing the easyaerosol
!>                climatology.
!> @param[in]     nlayers          The number of layers
!> @param[in]     n_band           Number of spectral bands (sw or lw)
!> @param[in]     easy_prop        Aerosol property from easyaerosol
!> @param[in,out] easy_prop_diag   Aerosol property diagnostic
!> @param[in]     ndf_easy         No. of DOFs per cell for easy space
!> @param[in]     undf_easy        No. unique of DOFs for easy space
!> @param[in]     map_easy         Dofmap for aer space column base cell
!> @param[in]     ndf_diag         No. of DOFs per cell for diag space
!> @param[in]     undf_diag        No. unique of DOFs for diag space
!> @param[in]     map_diag         Dofmap for diag ce column base cell
!
subroutine easyaerosol_fld2diag_code( nlayers, n_band,                         &
                                      easy_prop,                               &
                                      easy_diag,                               &
                                      ndf_easy, undf_easy, map_easy,           &
                                      ndf_diag, undf_diag, map_diag)

  use constants_mod, only: r_def, i_def, i_um

  implicit none

  ! Arguments

  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: n_band

  integer(kind=i_def), intent(in) :: ndf_easy
  integer(kind=i_def), intent(in) :: undf_easy
  integer(kind=i_def), dimension(ndf_easy), intent(in) :: map_easy

  integer(kind=i_def), intent(in) :: ndf_diag
  integer(kind=i_def), intent(in) :: undf_diag
  integer(kind=i_def), dimension(ndf_diag), intent(in) :: map_diag

  real(kind=r_def), intent(in),    dimension(undf_easy) :: easy_prop
  real(kind=r_def), intent(inout), dimension(undf_diag) :: easy_diag

  ! Loop counters / array index
  integer(i_um) :: k, i_band, i_easy, i_diag

  do i_band = 1, n_band
    k = 0
    i_diag  = map_easy(1) + ((i_band-1) * (nlayers+1)) + k
    easy_diag (i_diag) = 0.0_r_def

    do k = 1, nlayers
      !
      i_easy = map_easy(1) + (k * n_band) + i_band - 1
      i_diag = map_diag(1) + ((i_band-1) * (nlayers+1)) + k

      easy_diag (i_diag) =  easy_prop(i_easy)

    end do
  end do

end subroutine easyaerosol_fld2diag_code

end module easyaerosol_fld2diag_kernel_mod
