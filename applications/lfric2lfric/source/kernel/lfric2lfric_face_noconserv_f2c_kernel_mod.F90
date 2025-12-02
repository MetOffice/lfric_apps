!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief The weighted non-conservative restriction operation for scalar face fields.
!> @details Restrict a scalar field on a fine grid to a coarse grid field.
!!          Weights are provided, for instance to give conservative restriction.
!!          This kernel only works for the lowest-order Wtheta spaces.
!!          It is assumed that the non-conservative weights will be 2D weights.
!!          This is a copy of restric_scalar_weighted_kernel_mod, with a
!!          modification to take into account multidata values in grid points

module lfric2lfric_face_noconserv_f2c_kernel_mod

use constants_mod,           only: i_def, r_double, r_single
use kernel_mod,              only: kernel_type
use argument_mod,            only: arg_type,                                  &
                                   GH_FIELD, GH_REAL,                         &
                                   GH_READ, GH_WRITE,                         &
                                   ANY_DISCONTINUOUS_SPACE_1,                 &
                                   ANY_DISCONTINUOUS_SPACE_2,                 &
                                   GH_COARSE, GH_FINE, CELL_COLUMN

implicit none

private

type, public, extends(kernel_type) :: lfric2lfric_face_noconserv_f2c_kernel_type
   private
   type(arg_type) :: meta_args(3) = (/                                        &
        arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1,      &
                                              mesh_arg=GH_COARSE),            &
        arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2,      &
                                              mesh_arg=GH_FINE  ),            &
        arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2,      &
                                              mesh_arg=GH_FINE  )             &
        /)
  integer :: operates_on = CELL_COLUMN
end type lfric2lfric_face_noconserv_f2c_kernel_type

public :: lfric2lfric_face_noconserv_f2c_kernel_code

  ! Generic interface for real32 and real64 types
  interface lfric2lfric_face_noconserv_f2c_kernel_code
    module procedure                                                          &
      lfric2lfric_face_noconserv_f2c_code_r_single,                           &
      lfric2lfric_face_noconserv_f2c_code_r_double
  end interface

contains

  !> @brief The weighted restriction operation for scalar fields.
  !> @param[in]     nlayers                  Number of layers in a model column
  !> @param[in]     cell_map                 A 2D index map of which fine grid
  !!                                         cells lie in the coarse grid cell
  !> @param[in]     ncell_fine_per_coarse_x  Number of fine cells per coarse
  !!                                         cell in the horizontal x-direction
  !> @param[in]     ncell_fine_per_coarse_y  Number of fine cells per coarse
  !!                                         cell in the horizontal y-direction
  !> @param[in]     ncell_fine               Number of cells in the partition
  !!                                         for the fine grid
  !> @param[in,out] coarse_field             Coarse grid field to write to
  !> @param[in]     fine_field               The fine grid field to restrict
  !> @param[in]     weights                  Fine grid field containing weights
  !!                                         for the prolongation. The inverse
  !!                                         of these weights are used for the
  !!                                         restriction. The actual dimension
  !!                                         of this field will be smaller than
  !!                                         the default and care needs to be
  !!                                         taken not to go out of bounds
  !> @param[in]     undf_coarse              Total num of DoFs on the coarse
  !!                                         grid for this mesh partition
  !> @param[in]     map_coarse               DoFmap of cells on the coarse grid
  !> @param[in]     ndf                      Num of DoFs per cell on both grids
  !> @param[in]     undf_fine                Total num of DoFs on the fine grid
  !!                                         for this mesh partition
  !> @param[in]     map_fine                 DoFmap of cells on the fine grid

  ! R_SINGLE PRECISION
  ! ==================
  subroutine lfric2lfric_face_noconserv_f2c_code_r_single(                    &
                                                  nlayers,                    &
                                                  cell_map,                   &
                                                  ncell_fine_per_coarse_x,    &
                                                  ncell_fine_per_coarse_y,    &
                                                  ncell_fine,                 &
                                                  coarse_field,               &
                                                  fine_field,                 &
                                                  weights,                    &
                                                  undf_coarse,                &
                                                  map_coarse,                 &
                                                  ndf,                        &
                                                  undf_fine,                  &
                                                  map_fine)

    implicit none

    integer(kind=i_def), intent(in)    :: nlayers
    integer(kind=i_def), intent(in)    :: ncell_fine_per_coarse_x
    integer(kind=i_def), intent(in)    :: ncell_fine_per_coarse_y
    integer(kind=i_def), intent(in)    :: cell_map(ncell_fine_per_coarse_x, ncell_fine_per_coarse_y)
    integer(kind=i_def), intent(in)    :: ncell_fine
    integer(kind=i_def), intent(in)    :: ndf
    integer(kind=i_def), intent(in)    :: map_fine(ndf, ncell_fine)
    integer(kind=i_def), intent(in)    :: map_coarse(ndf)
    integer(kind=i_def), intent(in)    :: undf_fine, undf_coarse
    real(kind=r_single), intent(inout) :: coarse_field(undf_coarse)
    real(kind=r_single), intent(in)    :: fine_field(undf_fine)
    real(kind=r_single), intent(in)    :: weights(ncell_fine *             &
                                                 (nlayers -1 + ndf))

    integer(kind=i_def) :: multidata, df, k, m, x_idx, y_idx, top_df
    real(kind=r_single) :: denom

    denom = 1.0_r_single/real(ncell_fine_per_coarse_x*ncell_fine_per_coarse_y, kind=r_single)

    ! Assume lowest order W3 or Wtheta space
    df = 1
    ! Loop is 0 -> nlayers-1 for W3 fields, but 0 -> nlayers for Wtheta fields
    top_df = nlayers - 2 + ndf
    ! Number of multidata values per grid cell
    multidata = undf_fine/((top_df+1)*ncell_fine) - 1

    ! Build up 1D array of new coarse values for this column
    do y_idx = 1, ncell_fine_per_coarse_y
      do x_idx = 1, ncell_fine_per_coarse_x
        do m = 0, multidata
          do k = 0, top_df
            coarse_field(map_coarse(df)+m*(top_df+1)+k) = &
              coarse_field(map_coarse(df)+m*(top_df+1)+k) +         &
              fine_field(map_fine(df,cell_map(x_idx,y_idx))+m*(top_df+1)+k) * &
              denom * weights((map_fine(df,cell_map(x_idx,y_idx))-1)/((top_df+1)*(multidata+1))+1)
          end do
        end do
      end do
    end do

  end subroutine lfric2lfric_face_noconserv_f2c_code_r_single

  ! R_DOUBLE PRECISION
  ! ==================
  subroutine lfric2lfric_face_noconserv_f2c_code_r_double(                    &
                                                  nlayers,                    &
                                                  cell_map,                   &
                                                  ncell_fine_per_coarse_x,    &
                                                  ncell_fine_per_coarse_y,    &
                                                  ncell_fine,                 &
                                                  coarse_field,               &
                                                  fine_field,                 &
                                                  weights,                    &
                                                  undf_coarse,                &
                                                  map_coarse,                 &
                                                  ndf,                        &
                                                  undf_fine,                  &
                                                  map_fine)

    implicit none

    integer(kind=i_def), intent(in)    :: nlayers
    integer(kind=i_def), intent(in)    :: ncell_fine_per_coarse_x
    integer(kind=i_def), intent(in)    :: ncell_fine_per_coarse_y
    integer(kind=i_def), intent(in)    :: cell_map(ncell_fine_per_coarse_x,   &
                                                   ncell_fine_per_coarse_y)
    integer(kind=i_def), intent(in)    :: ncell_fine
    integer(kind=i_def), intent(in)    :: ndf
    integer(kind=i_def), intent(in)    :: map_fine(ndf, ncell_fine)
    integer(kind=i_def), intent(in)    :: map_coarse(ndf)
    integer(kind=i_def), intent(in)    :: undf_fine, undf_coarse
    real(kind=r_double), intent(inout) :: coarse_field(undf_coarse)
    real(kind=r_double), intent(in)    :: fine_field(undf_fine)
    real(kind=r_double), intent(in)    :: weights(ncell_fine *             &
                                                 (nlayers -1 + ndf))

    integer(kind=i_def) :: multidata, df, k, m, x_idx, y_idx, top_df
    real(kind=r_double) :: denom

    denom = 1.0_r_double/real(ncell_fine_per_coarse_x*ncell_fine_per_coarse_y, kind=r_double)

    ! Assume lowest order W3 or Wtheta space
    df = 1
    ! Loop is 0 -> nlayers-1 for W3 fields, but 0 -> nlayers for Wtheta fields
    top_df = nlayers - 2 + ndf
    ! Number of multidata values per grid cell
    multidata = undf_fine/((top_df+1)*ncell_fine) - 1

    ! Build up 1D array of new coarse values for this column
    do y_idx = 1, ncell_fine_per_coarse_y
      do x_idx = 1, ncell_fine_per_coarse_x
        do m = 0, multidata
          do k = 0, top_df
            coarse_field(map_coarse(df)+m*(top_df+1)+k) = &
              coarse_field(map_coarse(df)+m*(top_df+1)+k) +         &
              fine_field(map_fine(df,cell_map(x_idx,y_idx))+m*(top_df+1)+k) * &
              denom * weights((map_fine(df,cell_map(x_idx,y_idx))-1)/((top_df+1)*(multidata+1))+1)
          end do
        end do
      end do
    end do

  end subroutine lfric2lfric_face_noconserv_f2c_code_r_double

end module lfric2lfric_face_noconserv_f2c_kernel_mod
