!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Computes coefficients for fourth-order FFSL edge reconstruction.
!!        Note that this kernel only works when field is a W3 field at lowest
!!        order, since it is assumed that ndf_w3 = 1.

module ffsl_fourth_dldz_kernel_mod

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
type, public, extends(kernel_type) :: ffsl_fourth_dldz_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                &
       arg_type(GH_FIELD*4, GH_REAL, GH_WRITE, W2v), & ! dl_dz
       arg_type(GH_FIELD,   GH_REAL, GH_READ,  W3)   & ! dz
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: ffsl_fourth_dldz_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: ffsl_fourth_dldz_code

contains

!> @brief Computes coefficients for fourth-order FFSL edge reconstruction.
!> @param[in]     nlayers   Number of layers
!> @param[in,out] dl_dz_1   The first set of coefficients to be calculated
!> @param[in,out] dl_dz_2   The second set of coefficients to be calculated
!> @param[in,out] dl_dz_3   The third set of coefficients to be calculated
!> @param[in,out] dl_dz_4   The fourth set of coefficients to be calculated
!> @param[in]     dz        Vertical length of the W3 cell
!> @param[in]     ndf_w2v   Number of degrees of freedom for W2v per cell
!> @param[in]     undf_w2v  Number of unique degrees of freedom for W2v
!> @param[in]     map_w2v   The dofmap for the W2v cell at the base of the column
!> @param[in]     ndf_w3    Number of degrees of freedom for W3 per cell
!> @param[in]     undf_w3   Number of unique degrees of freedom for W3
!> @param[in]     map_w3    The dofmap for the cell at the base of the column
subroutine ffsl_fourth_dldz_code( nlayers,   &
                                  dl_dz_1,   &
                                  dl_dz_2,   &
                                  dl_dz_3,   &
                                  dl_dz_4,   &
                                  dz,        &
                                  ndf_w2v,   &
                                  undf_w2v,  &
                                  map_w2v,   &
                                  ndf_w3,    &
                                  undf_w3,   &
                                  map_w3 )

  use subgrid_vertical_support_mod, only: fourth_precompute_height

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: undf_w2v
  integer(kind=i_def), intent(in)    :: ndf_w2v
  integer(kind=i_def), intent(in)    :: undf_w3
  integer(kind=i_def), intent(in)    :: ndf_w3
  real(kind=r_tran),   intent(inout) :: dl_dz_1(undf_w2v)
  real(kind=r_tran),   intent(inout) :: dl_dz_2(undf_w2v)
  real(kind=r_tran),   intent(inout) :: dl_dz_3(undf_w2v)
  real(kind=r_tran),   intent(inout) :: dl_dz_4(undf_w2v)
  real(kind=r_tran),   intent(in)    :: dz(undf_w3)
  integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)
  integer(kind=i_def), intent(in)    :: map_w2v(ndf_w2v)

  ! Internal variables
  integer(kind=i_def) :: k, w2v_idx, w3_idx, local_edge_idx

  w2v_idx = map_w2v(1)
  w3_idx = map_w3(1)

  ! Lowest cells ---------------------------------------------------------------
  local_edge_idx = 0
  call fourth_precompute_height(                                               &
          dl_dz_1(w2v_idx+0), dl_dz_2(w2v_idx+0),                              &
          dl_dz_3(w2v_idx+0), dl_dz_4(w2v_idx+0),                              &
          dz(w3_idx+0 : w3_idx+3), local_edge_idx                              &
  )

  local_edge_idx = 1
  call fourth_precompute_height(                                               &
          dl_dz_1(w2v_idx+1), dl_dz_2(w2v_idx+1),                              &
          dl_dz_3(w2v_idx+1), dl_dz_4(w2v_idx+1),                              &
          dz(w3_idx+0 : w3_idx+3), local_edge_idx                              &
  )

  ! Interior cells -------------------------------------------------------------
  local_edge_idx = 2
  do k = 2, nlayers - 2
    call fourth_precompute_height(                                             &
            dl_dz_1(w2v_idx+k), dl_dz_2(w2v_idx+k),                            &
            dl_dz_3(w2v_idx+k), dl_dz_4(w2v_idx+k),                            &
            dz(w3_idx+k-2 : w3_idx+k+1), local_edge_idx                        &
    )
  end do

  ! Highest cells --------------------------------------------------------------
  local_edge_idx = 3
  k = nlayers - 1
  call fourth_precompute_height(                                               &
          dl_dz_1(w2v_idx+k), dl_dz_2(w2v_idx+k),                              &
          dl_dz_3(w2v_idx+k), dl_dz_4(w2v_idx+k),                              &
          dz(w3_idx+k-3 : w3_idx+k), local_edge_idx                            &
  )

  local_edge_idx = 4
  k = nlayers
  call fourth_precompute_height(                                               &
          dl_dz_1(w2v_idx+k), dl_dz_2(w2v_idx+k),                              &
          dl_dz_3(w2v_idx+k), dl_dz_4(w2v_idx+k),                              &
          dz(w3_idx+k-4 : w3_idx+k-1), local_edge_idx                          &
  )


end subroutine ffsl_fourth_dldz_code

end module ffsl_fourth_dldz_kernel_mod
