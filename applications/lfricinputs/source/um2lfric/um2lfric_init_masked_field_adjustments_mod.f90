! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module um2lfric_init_masked_field_adjustments_mod

! Intrinsic modules
use, intrinsic :: iso_fortran_env, only: int32, int64, real64

implicit none

private

public :: um2lfric_init_masked_field_adjustments

contains

!> @brief     Initialise the masked field adjustments ready for post-processing
!> @details   Read land/sea masks from the UM file and LFRic ancils and use this
!!            information to determine which points are not completely on one
!!            or the other and so need post-processing after regridding.
!!            Nearest neighbour algorithms are used to determine which data
!!            should be used instead.
subroutine um2lfric_init_masked_field_adjustments()

use lfricinp_nearest_neighbour_mod,    only: find_nn_on_um_grid
use lfricinp_get_latlon_mod,           only: get_lfric_mesh_coords
use lfricinp_masks_mod,                only: lfricinp_init_masks,              &
                                             lfricinp_finalise_masks,          &
                                             um_land_mask,                     &
                                             um_maritime_mask,                 &
                                             lfric_land_mask,                  &
                                             lfric_maritime_mask
use lfricinp_stashmaster_mod,          only: stashcode_land_frac

use um2lfric_regrid_weights_mod,           only: get_weights
use um2lfric_masked_field_adjustments_mod, only: land_field_adjustments,       &
                                                 maritime_field_adjustments

implicit none

! Local variables
character(len=1),    parameter :: um_land_mask_grid_type = 'p'
integer(kind=int32)            :: cell_lid, idx_lon_nn, idx_lat_nn, i
real(kind=real64)              :: lon, lat

! Initialise masks
call lfricinp_init_masks(stashcode_land_frac)

!
! Initialise land field adjustments
!
if (allocated(um_land_mask) .and. allocated(lfric_land_mask)) then

  ! Find indices of land points on LFRic mesh that will require adjustment to
  ! UM NN land point values
  call land_field_adjustments%find_adjusted_points_src_2d_dst_1d(              &
                                      src_mask=um_land_mask,                   &
                                      dst_mask=lfric_land_mask,                &
                                      weights=get_weights(stashcode_land_frac))
  !
  ! Set map from LFRic adjusted land point indices to UM NN land point indices
  allocate(land_field_adjustments%adjusted_dst_to_src_map_2D(                  &
                                  land_field_adjustments%num_adjusted_points,2 &
                                                            ))
  do i = 1, land_field_adjustments%num_adjusted_points
    cell_lid = land_field_adjustments%adjusted_dst_indices_1D(i)
    call get_lfric_mesh_coords(cell_lid, lon, lat)
    call find_nn_on_um_grid(um_mask=um_land_mask,                              &
                            um_mask_grid_type=um_land_mask_grid_type,          &
                            lon_ref=lon,                                       &
                            lat_ref=lat,                                       &
                            idx_lon_nn=idx_lon_nn,                             &
                            idx_lat_nn=idx_lat_nn)
    land_field_adjustments%adjusted_dst_to_src_map_2D(i,1) = idx_lon_nn
    land_field_adjustments%adjusted_dst_to_src_map_2D(i,2) = idx_lat_nn
  end do
  !
  ! Set initialisation flag of land field adjustment
  land_field_adjustments%initialised = .true.

end if

!
! Initialise maritime field adjustments
!
if (allocated(um_maritime_mask) .and. allocated(lfric_maritime_mask)) then

  ! Find indices of maritime points on LFRic mesh that will require adjustment
  ! to UM NN maritime point values
  call maritime_field_adjustments%find_adjusted_points_src_2d_dst_1d(          &
                                      src_mask=um_maritime_mask,               &
                                      dst_mask=lfric_maritime_mask,            &
                                      weights=get_weights(stashcode_land_frac))
  !
  ! Set map from LFRic adjusted maritime point indices to UM NN maritime point
  ! indices
  allocate(maritime_field_adjustments%adjusted_dst_to_src_map_2D(              &
                              maritime_field_adjustments%num_adjusted_points,2 &
                                                            ))
  do i = 1, maritime_field_adjustments%num_adjusted_points
    cell_lid = maritime_field_adjustments%adjusted_dst_indices_1D(i)
    call get_lfric_mesh_coords(cell_lid, lon, lat)
    call find_nn_on_um_grid(um_mask=um_maritime_mask,                          &
                            um_mask_grid_type=um_land_mask_grid_type,          &
                            lon_ref=lon,                                       &
                            lat_ref=lat,                                       &
                            idx_lon_nn=idx_lon_nn,                             &
                            idx_lat_nn=idx_lat_nn)
    maritime_field_adjustments%adjusted_dst_to_src_map_2D(i,1) = idx_lon_nn
    maritime_field_adjustments%adjusted_dst_to_src_map_2D(i,2) = idx_lat_nn
  end do
  !
  ! Set initialisation flag of maritime field adjustment
  maritime_field_adjustments%initialised = .true.

end if

! Finalise masks
call lfricinp_finalise_masks

end subroutine um2lfric_init_masked_field_adjustments

end module um2lfric_init_masked_field_adjustments_mod
