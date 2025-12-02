!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Computes coefficients for third-order FFSL edge reconstruction.
!!        Note that this kernel only works when field is a W3 field at lowest
!!        order, since it is assumed that ndf_w3 = 1.

module ffsl_third_dldz_kernel_mod

use argument_mod,                   only : arg_type,              &
                                           GH_FIELD, GH_REAL,     &
                                           GH_READ, GH_WRITE,     &
                                           CELL_COLUMN
use fs_continuity_mod,              only : W3, W2v
use constants_mod,                  only : r_tran, i_def, l_def, EPS_R_TRAN
use kernel_mod,                     only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: ffsl_third_dldz_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                &
       arg_type(GH_FIELD*3, GH_REAL, GH_WRITE, W3),  & ! dla_dz
       arg_type(GH_FIELD*3, GH_REAL, GH_WRITE, W3),  & ! dlb_dz
       arg_type(GH_FIELD,   GH_REAL, GH_READ,  W3)   & ! dz
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: ffsl_third_dldz_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: ffsl_third_dldz_code

contains

!> @brief Computes coefficients for third-order FFSL edge reconstruction.
!> @param[in]     nlayers   Number of layers
!> @param[in,out] dla_dz_1  The first set of coefficients to be calculated, for
!!                          reconstruction of the edge above each cell
!> @param[in,out] dla_dz_2  The second set of coefficients to be calculated, for
!!                          reconstruction of the edge above each cell
!> @param[in,out] dla_dz_3  The third set of coefficients to be calculated, for
!!                          reconstruction of the edge above each cell
!> @param[in,out] dlb_dz_1  The first set of coefficients to be calculated, for
!!                          reconstruction of the edge below each cell
!> @param[in,out] dlb_dz_2  The second set of coefficients to be calculated, for
!!                          reconstruction of the edge below each cell
!> @param[in,out] dlb_dz_3  The third set of coefficients to be calculated, for
!!                          reconstruction of the edge below each cell
!> @param[in]     dz        Vertical length of the W3 cell
!> @param[in]     ndf_w3    Number of degrees of freedom for W3 per cell
!> @param[in]     undf_w3   Number of unique degrees of freedom for W3
!> @param[in]     map_w3    The dofmap for the cell at the base of the column
subroutine ffsl_third_dldz_code( nlayers,   &
                                 dla_dz_1,   &
                                 dla_dz_2,   &
                                 dla_dz_3,   &
                                 dlb_dz_1,   &
                                 dlb_dz_2,   &
                                 dlb_dz_3,   &
                                 dz,         &
                                 ndf_w3,     &
                                 undf_w3,    &
                                 map_w3 )

  use subgrid_vertical_support_mod, only: third_precompute_height

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: undf_w3
  integer(kind=i_def), intent(in)    :: ndf_w3
  real(kind=r_tran),   intent(inout) :: dla_dz_1(undf_w3)
  real(kind=r_tran),   intent(inout) :: dla_dz_2(undf_w3)
  real(kind=r_tran),   intent(inout) :: dla_dz_3(undf_w3)
  real(kind=r_tran),   intent(inout) :: dlb_dz_1(undf_w3)
  real(kind=r_tran),   intent(inout) :: dlb_dz_2(undf_w3)
  real(kind=r_tran),   intent(inout) :: dlb_dz_3(undf_w3)
  real(kind=r_tran),   intent(in)    :: dz(undf_w3)
  integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)

  ! Internal variables
  integer(kind=i_def) :: k, w3_idx, local_cell_idx

  w3_idx = map_w3(1)

  ! Bottom layer ---------------------------------------------------------------
  local_cell_idx = 1
  call third_precompute_height(                                                &
          dla_dz_1(w3_idx+0), dla_dz_2(w3_idx+0), dla_dz_3(w3_idx+0),          &
          dlb_dz_1(w3_idx+0), dlb_dz_2(w3_idx+0), dlb_dz_3(w3_idx+0),          &
          dz(w3_idx+0 : w3_idx+2), local_cell_idx                              &
  )

  ! Centre layers --------------------------------------------------------------
  local_cell_idx = 2
  do k = 1, nlayers - 2
    call third_precompute_height(                                              &
            dla_dz_1(w3_idx+k), dla_dz_2(w3_idx+k), dla_dz_3(w3_idx+k),        &
            dlb_dz_1(w3_idx+k), dlb_dz_2(w3_idx+k), dlb_dz_3(w3_idx+k),        &
            dz(w3_idx+k-1 : w3_idx+k+1), local_cell_idx                        &
    )
  end do

  ! Top layer ------------------------------------------------------------------
  local_cell_idx = 3
  k = nlayers - 1
  call third_precompute_height(                                                &
          dla_dz_1(w3_idx+k), dla_dz_2(w3_idx+k), dla_dz_3(w3_idx+k),          &
          dlb_dz_1(w3_idx+k), dlb_dz_2(w3_idx+k), dlb_dz_3(w3_idx+k),          &
          dz(w3_idx+k-2 : w3_idx+k), local_cell_idx                            &
  )

end subroutine ffsl_third_dldz_code

end module ffsl_third_dldz_kernel_mod
