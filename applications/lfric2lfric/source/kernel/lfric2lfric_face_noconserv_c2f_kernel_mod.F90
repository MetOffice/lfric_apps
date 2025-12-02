!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Perform the non-conservative weighted injection-prolongation operation
!!        of a coarse grid scalar to a fine grid face
!> @details Prolong the coarse grid correction into all cells on the fine grid
!!          that are contained in that coarse grid cell. This is not an injection
!!          on the "mass" field, so that the prolongation is non-conservative.
!!          This kernel only works for the lowest-order face spaces.
!!          It is assumed that the non-conservative weights will be 2D weights.
!!          This is a copy of prolong_scalar_weighted_kernel_mod, with a
!!          modification to take into account multidata values in grid points

module lfric2lfric_face_noconserv_c2f_kernel_mod

use constants_mod,           only: i_def, r_double, r_single
use kernel_mod,              only: kernel_type
use argument_mod,            only: arg_type,                                  &
                                   GH_FIELD, GH_REAL,                         &
                                   GH_READ, GH_READWRITE,                     &
                                   ANY_DISCONTINUOUS_SPACE_1,                 &
                                   ANY_DISCONTINUOUS_SPACE_2,                 &
                                   GH_COARSE, GH_FINE, CELL_COLUMN

implicit none

private

type, public, extends(kernel_type) :: lfric2lfric_face_noconserv_c2f_kernel_type
   private
   type(arg_type) :: meta_args(3) = (/                                        &
        arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1,  &
                                                  mesh_arg=GH_FINE),          &
        arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_2,  &
                                                  mesh_arg=GH_COARSE ),       &
        arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1,  &
                                                  mesh_arg=GH_FINE )          &
        /)
  integer :: operates_on = CELL_COLUMN
end type lfric2lfric_face_noconserv_c2f_kernel_type

public :: lfric2lfric_face_noconserv_c2f_kernel_code

  ! Generic interface for real32 and real64 types
  interface lfric2lfric_face_noconserv_c2f_kernel_code
    module procedure  &
      lfric2lfric_face_noconserv_c2f_code_r_single,                           &
      lfric2lfric_face_noconserv_c2f_code_r_double
  end interface

contains

  !> @brief Performs the weighted injection-prolongation for scalar fields
  !> @param[in]     nlayers                  Number of layers in a model column
  !> @param[in]     cell_map                 A 2D index map of which fine grid
  !!                                         cells lie in the coarse grid cell
  !> @param[in]     ncell_fine_per_coarse_x  Number of fine cells per coarse
  !!                                         cell in the horizontal x-direction
  !> @param[in]     ncell_fine_per_coarse_y  Number of fine cells per coarse
  !!                                         cell in the horizontal y-direction
  !> @param[in]     ncell_fine               Number of cells in the partition
  !!                                         for the fine grid
  !> @param[in,out] fine_field               The fine grid field to write to
  !> @param[in]     coarse_field             Coarse grid field to prolong
  !> @param[in]     weights                  Fine grid field containing weights
  !!                                         for the prolongation. The actual
  !!                                         dimension of this field will be
  !!                                         smaller than the default and care
  !!                                         needs to be taken not go go out
  !!                                         of bounds
  !> @param[in]     ndf                      Num of DoFs per cell on both grids
  !> @param[in]     undf_fine                Total num of DoFs on the fine grid
  !!                                         for this mesh partition
  !> @param[in]     map_fine                 DoFmap of cells on the fine grid
  !> @param[in]     undf_coarse              Total num of DoFs on the coarse
  !!                                         grid for this mesh partition
  !> @param[in]     map_coarse               DoFmap of cells on the coarse grid

  ! R_SINGLE PRECISION
  ! ==================
  subroutine lfric2lfric_face_noconserv_c2f_code_r_single(                    &
                                          nlayers,                            &
                                          cell_map,                           &
                                          ncell_fine_per_coarse_x,            &
                                          ncell_fine_per_coarse_y,            &
                                          ncell_fine,                         &
                                          fine_field,                         &
                                          coarse_field,                       &
                                          weights,                            &
                                          ndf,                                &
                                          undf_fine,                          &
                                          map_fine,                           &
                                          undf_coarse,                        &
                                          map_coarse)

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
    real(kind=r_single), intent(in)    :: coarse_field(undf_coarse)
    real(kind=r_single), intent(inout) :: fine_field(undf_fine)
    real(kind=r_single), intent(in)    :: weights(ncell_fine)

    integer(kind=i_def) :: multidata, df, k, m, x_idx, y_idx, top_df

    ! Assume lowest order W3 or Wtheta space
    df = 1
    ! Loop is 0 -> nlayers-1 for W3 fields, but 0 -> nlayers for Wtheta fields
    top_df = nlayers - 2 + ndf
    ! Number of multidata values per grid cell
    multidata = undf_fine/((top_df+1)*ncell_fine) - 1

    do y_idx = 1, ncell_fine_per_coarse_y
      do x_idx = 1, ncell_fine_per_coarse_x
        do m = 0, multidata
          do k = 0, top_df
            fine_field(map_fine(df,cell_map(x_idx,y_idx))+m*(top_df+1)+k) =   &
              fine_field(map_fine(df,cell_map(x_idx,y_idx))+m*(top_df+1)+k) + &
              weights((map_fine(df,cell_map(x_idx,y_idx))-1)/((top_df+1)*(multidata+1))+1) * &
              coarse_field(map_coarse(df)+m*(top_df+1)+k)
          end do
        end do
      end do
    end do

  end subroutine lfric2lfric_face_noconserv_c2f_code_r_single

  ! R_DOUBLE PRECISION
  ! ==================
  subroutine lfric2lfric_face_noconserv_c2f_code_r_double(                    &
                                          nlayers,                            &
                                          cell_map,                           &
                                          ncell_fine_per_coarse_x,            &
                                          ncell_fine_per_coarse_y,            &
                                          ncell_fine,                         &
                                          fine_field,                         &
                                          coarse_field,                       &
                                          weights,                            &
                                          ndf,                                &
                                          undf_fine,                          &
                                          map_fine,                           &
                                          undf_coarse,                        &
                                          map_coarse)

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
    real(kind=r_double), intent(in)    :: coarse_field(undf_coarse)
    real(kind=r_double), intent(inout) :: fine_field(undf_fine)
    real(kind=r_double), intent(inout) :: weights(ncell_fine)

    integer(kind=i_def) :: multidata, df, k, m, x_idx, y_idx, top_df

    ! Assume lowest order W3 or Wtheta space
    df = 1
    ! Loop is 0 -> nlayers-1 for W3 fields, but 0 -> nlayers for Wtheta fields
    top_df = nlayers - 2 + ndf
    ! Number of multidata values per grid cell
    multidata = undf_fine/((top_df+1)*ncell_fine) - 1

    do y_idx = 1, ncell_fine_per_coarse_y
      do x_idx = 1, ncell_fine_per_coarse_x
        do m = 0, multidata
          do k = 0, top_df
            fine_field(map_fine(df,cell_map(x_idx,y_idx))+m*(top_df+1)+k) =   &
              fine_field(map_fine(df,cell_map(x_idx,y_idx))+m*(top_df+1)+k) + &
              weights((map_fine(df,cell_map(x_idx,y_idx))-1)/((top_df+1)*(multidata+1))+1) * &
              coarse_field(map_coarse(df)+m*(top_df+1)+k)
          end do
        end do
      end do
    end do

  end subroutine lfric2lfric_face_noconserv_c2f_code_r_double

end module lfric2lfric_face_noconserv_c2f_kernel_mod
