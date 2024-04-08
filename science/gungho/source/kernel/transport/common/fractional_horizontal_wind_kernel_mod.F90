!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes the fractional part of the horizontal wind.
!> @details Computes the fractional part of the wind from the horizontal
!!          departure distance, and the volume of the upwind cell.

module fractional_horizontal_wind_kernel_mod

use argument_mod,       only : arg_type, GH_INTEGER,    &
                               GH_FIELD, GH_REAL,       &
                               CELL_COLUMN, GH_WRITE,   &
                               GH_READ, STENCIL, CROSS, &
                               ANY_DISCONTINUOUS_SPACE_3
use constants_mod,      only : i_def, r_tran
use fs_continuity_mod,  only : W3, W2h
use kernel_mod,         only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
type, public, extends(kernel_type) :: fractional_horizontal_wind_kernel_type
  private
  type(arg_type) :: meta_args(5) = (/                                       &
      arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, W2h),                       & ! frac_wind_xy
      arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2h),                       & ! dep_dist_xy
      arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3, STENCIL(CROSS)),        & ! detj
      arg_type(GH_FIELD,  GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! face_selector_ew
      arg_type(GH_FIELD,  GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3)  & ! face_selector_ns
      /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: fractional_horizontal_wind_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: fractional_horizontal_wind_code

contains

!> @brief Computes the fractional part of the horizontal wind.
!> @param[in]     nlayers           Number of layers
!> @param[in,out] frac_wind_xy      Fractional wind
!> @param[in]     dep_dist_xy       Departure distance
!> @param[in]     detj              Volume factor at W3
!> @param[in]     stencil_size      Local length of Det(J) at W3 stencil
!> @param[in]     stencil_map       Dofmap for the Det(J) at W3 stencil
!> @param[in]     face_selector_ew  2D field indicating which faces to loop over
!> @param[in]     face_selector_ns  2D field indicating which faces to loop over
!> @param[in]     ndf_w2h           Number of DoFs per cell for W2H
!> @param[in]     undf_w2h          Num of W2H DoFs in memory for this partition
!> @param[in]     map_w2h           Map of lowest-cell W2H DoFs
!> @param[in]     ndf_w3            Number of DoFs per cell for W3
!> @param[in]     undf_w3           Num of W3 DoFs in memory for this partition
!> @param[in]     map_w3            Map of lowest-cell W3 DoFs
!> @param[in]     ndf_w3_2d         Number of DoFs for 2D W3 per cell
!> @param[in]     undf_w3_2d        Number of DoFs for this partition for 2D W3
!> @param[in]     map_w3_2d         Map for 2D W3
subroutine fractional_horizontal_wind_code( nlayers,          &
                                            frac_wind_xy,     &
                                            dep_dist_xy ,     &
                                            detj,             &
                                            stencil_size,     &
                                            stencil_map,      &
                                            face_selector_ew, &
                                            face_selector_ns, &
                                            ndf_w2h,          &
                                            undf_w2h,         &
                                            map_w2h,          &
                                            ndf_w3,           &
                                            undf_w3,          &
                                            map_w3,           &
                                            ndf_w3_2d,        &
                                            undf_w3_2d,       &
                                            map_w3_2d )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: undf_w3
  integer(kind=i_def), intent(in) :: ndf_w3
  integer(kind=i_def), intent(in) :: undf_w3_2d
  integer(kind=i_def), intent(in) :: ndf_w3_2d
  integer(kind=i_def), intent(in) :: undf_w2h
  integer(kind=i_def), intent(in) :: ndf_w2h
  integer(kind=i_def), intent(in) :: stencil_size

  ! Arguments: Maps
  integer(kind=i_def), dimension(ndf_w3),    intent(in) :: map_w3
  integer(kind=i_def), dimension(ndf_w2h),   intent(in) :: map_w2h
  integer(kind=i_def), dimension(ndf_w3_2d), intent(in) :: map_w3_2d
  integer(kind=i_def), dimension(ndf_w3,stencil_size), intent(in) :: stencil_map

  ! Arguments: Fields
  real(kind=r_tran),   dimension(undf_w2h),   intent(inout) :: frac_wind_xy
  real(kind=r_tran),   dimension(undf_w2h),   intent(in)    :: dep_dist_xy
  real(kind=r_tran),   dimension(undf_w3),    intent(in)    :: detj
  integer(kind=i_def), dimension(undf_w3_2d), intent(in)    :: face_selector_ew
  integer(kind=i_def), dimension(undf_w3_2d), intent(in)    :: face_selector_ns

  ! Variables
  integer(kind=i_def) :: k, df, w2h_idx, df_idx
  integer(kind=i_def) :: pos_stencil_idx, neg_stencil_idx
  real(kind=r_tran)   :: frac_dist, dep_dist
  ! The stencil index of the upwind cell, for positive and negative winds
  ! The index of these arrays is the DoF (W, S, E, N)
  integer(kind=i_def), parameter :: pos_wind_cells(4) = (/ 2, 1, 1, 5 /)
  integer(kind=i_def), parameter :: neg_wind_cells(4) = (/ 1, 3, 4, 1 /)
  integer(kind=i_def), parameter :: x_dofs(2) = (/ 1, 3 /)
  integer(kind=i_def), parameter :: y_dofs(2) = (/ 2, 4 /)

  ! Cross stencil extent 1 has form
  !     | 5 |
  ! | 2 | 1 | 4 |
  !     | 3 |

  ! W2 has horizontal DOF map
  !       4
  !     1   3
  !       2

  ! Check if at edge of LAM: if so, do nothing
  ! The frac_wind fields should already have been initialised to 0
  if (stencil_size == 5) then

    ! Loop through E/W DoFs
    do df_idx = 1, face_selector_ew(map_w3_2d(1))
      df = x_dofs(df_idx)
      w2h_idx = map_w2h(df)
      pos_stencil_idx = stencil_map(1, pos_wind_cells(df))
      neg_stencil_idx = stencil_map(1, neg_wind_cells(df))

      do k = 0, nlayers-1
        dep_dist = dep_dist_xy(w2h_idx + k)
        frac_dist = dep_dist - REAL(INT(dep_dist, i_def), r_tran)

        ! Add both contributions, but such that only the upwind one is non-zero
        frac_wind_xy(w2h_idx + k) = frac_dist *                                &
          ((0.5_r_tran + SIGN(0.5_r_tran, dep_dist))*detj(pos_stencil_idx + k) &
           + (0.5_r_tran - SIGN(0.5_r_tran, dep_dist))*detj(neg_stencil_idx + k))
      end do
    end do

    ! Loop through N/S DoFs
    do df_idx = 1, face_selector_ns(map_w3_2d(1))
      df = y_dofs(df_idx)
      w2h_idx = map_w2h(df)
      pos_stencil_idx = stencil_map(1, pos_wind_cells(df))
      neg_stencil_idx = stencil_map(1, neg_wind_cells(df))

      do k = 0, nlayers-1
        dep_dist = dep_dist_xy(w2h_idx + k)
        frac_dist = dep_dist - REAL(INT(dep_dist, i_def), r_tran)

        ! Add both contributions, but such that only the upwind one is non-zero
        frac_wind_xy(w2h_idx + k) = frac_dist *                                &
          ((0.5_r_tran + SIGN(0.5_r_tran, dep_dist))*detj(pos_stencil_idx + k) &
           + (0.5_r_tran - SIGN(0.5_r_tran, dep_dist))*detj(neg_stencil_idx + k))
      end do
    end do

  end if  ! Not at edge of LAM

end subroutine fractional_horizontal_wind_code

end module fractional_horizontal_wind_kernel_mod

