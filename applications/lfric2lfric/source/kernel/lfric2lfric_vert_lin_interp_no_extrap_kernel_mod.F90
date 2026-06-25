!-----------------------------------------------------------------------------
! (C) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief SHARKS Perform the conservative weighted injection-prolongation operation
!!        of a coarse grid scalar to a fine grid face field
!> @details SHARKS Prolong the coarse grid correction into all cells on the fine grid
!!          that are contained in that coarse grid cell. This is an injection
!!          on the "mass" field, so that the prolongation is conservative.
!!          This kernel only works for the lowest-order face spaces.
!!          This is a copy of prolong_scalar_weighted_kernel_mod, with a
!!          modification to take into account multidata values in grid points

module lfric2lfric_vert_lin_interp_no_extrap_kernel_mod

use constants_mod,           only: i_def, r_double, r_single
use kernel_mod,              only: kernel_type
use argument_mod,            only: arg_type, CELL_COLUMN,                     &
                                   GH_FIELD, GH_INTEGER, GH_REAL,             &
                                   GH_READ, GH_READWRITE, GH_SCALAR,          &
                                   ANY_DISCONTINUOUS_SPACE_1,                 &
                                   ANY_DISCONTINUOUS_SPACE_2
use fs_continuity_mod,       only : W3

implicit none

private

type, public, extends(kernel_type) :: lfric2lfric_vert_lin_interp_no_extrap_kernel_type

   type(arg_type) :: meta_args(6) = (/                                         &
        arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), &
        arg_type(GH_FIELD,  GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_2), &
        arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                              &
        arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                              &
        arg_type(GH_FIELD,  GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1),      &
        arg_type(GH_FIELD,  GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_2)       &
        /)
  integer :: operates_on = CELL_COLUMN
end type lfric2lfric_vert_lin_interp_no_extrap_kernel_type

public :: lfric2lfric_vert_lin_interp_no_extrap_kernel_code, lfric2lfric_vert_lin_interp_no_extrap_code_r_single!, lfric2lfric_vert_lin_interp_no_extrap_code_r_double

  ! Generic interface for real32 and real64 types
  interface lfric2lfric_vert_lin_interp_no_extrap_kernel_code
    module procedure  &
      lfric2lfric_vert_lin_interp_no_extrap_code_r_single !,   &
      !lfric2lfric_vert_lin_interp_no_extrap_code_r_double
  end interface

contains

!> @brief Performs the weighted injection-prolongation for scalar fields
!> @param[in]     nlayers                  Number of layers in a model column (in the destination field)
!> @param[in,out] destination_field        The destination field to interpolate onto
!> @param[in]     source_field             Source field to be interpolated
!> @param[in]     source_layers            Number of layers in the source field
!> @param[in]     ncell                    Number of cells in the partition
!!                                         of both grids
!> @param[in]     source_heights           Heights of vertical layers in the source field
!> @param[in]     dest_heights             Heights of vertical layers in the destination field
!> @param[in]     ndf_dest                 Num of DoFs per cell on the destination grid
!> @param[in]     undf_dest                Total num of DoFs on the destination grid
!!                                         for this mesh partition
!> @param[in]     map_dest                 DoFmap of cells on the destination grid
!> @param[in]     ndf_source               Num of DoFs per cell on the source grid
!> @param[in]     undf_source              Total num of DoFs on the source
!!                                         grid for this mesh partition
!> @param[in]     map_source               DoFmap of cells on the source grid

! R_SINGLE PRECISION
! ==================
subroutine lfric2lfric_vert_lin_interp_no_extrap_code_r_single(  &
                                          nlayers,                 &
                                          destination_field,       &
                                          source_field,            &
                                          source_layers,           &
                                          ncell,                   &
                                          dest_heights,            &
                                          source_heights,          &
                                          ndf_dest,                &
                                          undf_dest,               &
                                          map_dest,                &
                                          ndf_source,              &
                                          undf_source,             &
                                          map_source)

  implicit none

  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ncell
  integer(kind=i_def), intent(in)    :: source_layers
  integer(kind=i_def), intent(in)    :: ndf_dest, ndf_source
  integer(kind=i_def), intent(in)    :: undf_dest, undf_source
  integer(kind=i_def), intent(in), dimension(ndf_dest)    :: map_dest(ndf_dest)
  integer(kind=i_def), intent(in), dimension(ndf_source)  :: map_source(ndf_source)
  real(kind=r_single), intent(inout), dimension(ndf_dest) :: destination_field(undf_dest)
  real(kind=r_single), intent(in), dimension(ndf_source)  :: source_field(undf_source)
  real(kind=r_single), intent(in), dimension(ndf_dest) :: dest_heights(undf_dest)
  real(kind=r_single), intent(in), dimension(ndf_source)  :: source_heights(undf_source)

  integer(kind=i_def) :: multidata, df, k, m, kk, level_below(nlayers), source_top_df, dest_top_df
  integer(kind=i_def) :: d_h, s_h_top, s_h_bottom, s_h_below

  ! Assume lowest order W3 or Wtheta space
  df = 1
  ! Loop is 0 -> nlayers-1 for W3 fields, but 0 -> nlayers for Wtheta fields
  dest_top_df = nlayers - 2 + ndf_dest
  source_top_df = source_layers - 2 + ndf_source
  ! Number of multidata values per grid cell
  multidata = undf_dest/((dest_top_df+1)*ncell) - 1

  do kk = 1, dest_top_df
    level_below(kk) = source_layers
    do k = 1, source_top_df
      if ( (source_heights(k) > dest_heights(kk)) .AND.            &
           (level_below(kk) == source_layers) ) THEN
        level_below(kk) = k-1
          ! potential optimisation: start from level_below(kk-1)
      end if
    end do
  end do

  do m = 0, multidata
     do kk = 1, dest_top_df

      ! EXTRAPOLATION METHOD - ! No linear extrapolation at top or bottom

      d_h = map_dest(df) + m*(dest_top_df+1) + kk
      s_h_top = map_source(df) + m*(source_top_df+1) + source_top_df
      s_h_bottom = map_source(df) + m*(source_top_df+1)
      s_h_below = map_source(df) + m*(source_top_df+1) + level_below(kk)
      
      if (dest_heights(d_h) >= source_heights(s_h_top)) then

        ! Top: Set to top input data
        destination_field(d_h) = source_field(s_h_top)

      else if (dest_heights(d_h) <= source_heights(s_h_bottom)) then

        ! Bottom: Set to bottom input data
        destination_field(d_h) = source_field(s_h_bottom)
        
      else

      ! Linearly interpolate

      ! dk(kk) =  ( (dh(kk) - sh(lb(kk))) * sf(lb(kk)+1) - (dh(kk) - sh(lb(kk)+1)) * sf(lb(kk)) )
         !          / (sh(lb(kk)+1) - sh(lb(kk)))

          destination_field(d_h) =           &
                  ( (dest_heights(d_h) - source_heights(s_h_below )) * source_field(s_h_below + 1)    &
                  - (dest_heights(d_h) - source_heights(s_h_below + 1)) * source_field(s_h_below) )   &
                  / ( source_heights(s_h_below + 1)  - source_heights(s_h_below) )

      end if
    end do
  end do

  end subroutine lfric2lfric_vert_lin_interp_no_extrap_code_r_single

  ! R_DOUBLE PRECISION
  ! ==================
  !subroutine lfric2lfric_vert_lin_interp_no_extrap_code_r_double

!!!!!!! SHARKS COPY R SINGLE TO HERE, change real fields from r_single to r_double

  !end subroutine lfric2lfric_vert_lin_interp_no_extrap_code_r_double

end module lfric2lfric_vert_lin_interp_no_extrap_kernel_mod
