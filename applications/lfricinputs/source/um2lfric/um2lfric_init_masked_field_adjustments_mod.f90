! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE um2lfric_init_masked_field_adjustments_mod

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY: int32, int64, real64

IMPLICIT NONE

PRIVATE

PUBLIC :: um2lfric_init_masked_field_adjustments

CONTAINS

SUBROUTINE um2lfric_init_masked_field_adjustments()

USE lfricinp_nearest_neighbour_mod,    ONLY: find_nn_on_um_grid
USE lfricinp_get_latlon_mod,           ONLY: get_lfric_mesh_coords
USE lfricinp_masks_mod,                ONLY: lfricinp_init_masks,              &
                                             lfricinp_finalise_masks,          &
                                             um_land_mask,                     &
                                             um_maritime_mask,                 &
                                             lfric_land_mask,                  &
                                             lfric_maritime_mask

USE um2lfric_regrid_weights_mod,           ONLY: get_weights
USE um2lfric_masked_field_adjustments_mod, ONLY: land_field_adjustments,       &
                                                 maritime_field_adjustments

IMPLICIT NONE

! Local variables
CHARACTER(LEN=1),    PARAMETER :: um_land_mask_grid_type = 'p'
INTEGER(KIND=int64), PARAMETER :: stashcode_land_mask = 505
INTEGER(KIND=int32)            :: cell_lid, idx_lon_nn, idx_lat_nn, i
REAL(KIND=real64)              :: lon, lat

! Initialise masks
CALL lfricinp_init_masks(stashcode_land_mask)

!
! Initialise land field adjustments
!
IF (ALLOCATED(um_land_mask) .AND. ALLOCATED(lfric_land_mask)) THEN

  ! Find indices of land points on LFRic mesh that will require adjustment to
  ! UM NN land point values
  CALL land_field_adjustments%find_adjusted_points_src_2d_dst_1d(              &
                                      src_mask=um_land_mask,                   &
                                      dst_mask=lfric_land_mask,                &
                                      weights=get_weights(stashcode_land_mask))
  !
  ! Set map from LFRic adjusted land point indices to UM NN land point indices
  ALLOCATE(land_field_adjustments%adjusted_dst_to_src_map_2D(                  &
                                  land_field_adjustments%num_adjusted_points,2 &
                                                            ))
  DO i = 1, land_field_adjustments%num_adjusted_points
    cell_lid = land_field_adjustments%adjusted_dst_indices_1D(i)
    CALL get_lfric_mesh_coords(cell_lid, lon, lat)
    CALL find_nn_on_um_grid(um_mask=um_land_mask,                              &
                            um_mask_grid_type=um_land_mask_grid_type,          &
                            lon_ref=lon,                                       &
                            lat_ref=lat,                                       &
                            idx_lon_nn=idx_lon_nn,                             &
                            idx_lat_nn=idx_lat_nn)
    land_field_adjustments%adjusted_dst_to_src_map_2D(i,1) = idx_lon_nn
    land_field_adjustments%adjusted_dst_to_src_map_2D(i,2) = idx_lat_nn
  END DO
  !
  ! Set initialisation flag of land field adjustment
  land_field_adjustments%initialised = .TRUE.

END IF

!
! Initialise maritime field adjustments
!
IF (ALLOCATED(um_maritime_mask) .AND. ALLOCATED(lfric_maritime_mask)) THEN

  ! Find indices of maritime points on LFRic mesh that will require adjustment
  ! to UM NN maritime point values
  CALL maritime_field_adjustments%find_adjusted_points_src_2d_dst_1d(          &
                                      src_mask=um_maritime_mask,               &
                                      dst_mask=lfric_maritime_mask,            &
                                      weights=get_weights(stashcode_land_mask))
  !
  ! Set map from LFRic adjusted maritime point indices to UM NN maritime point
  ! indices
  ALLOCATE(maritime_field_adjustments%adjusted_dst_to_src_map_2D(              &
                              maritime_field_adjustments%num_adjusted_points,2 &
                                                            ))
  DO i = 1, maritime_field_adjustments%num_adjusted_points
    cell_lid = maritime_field_adjustments%adjusted_dst_indices_1D(i)
    CALL get_lfric_mesh_coords(cell_lid, lon, lat)
    CALL find_nn_on_um_grid(um_mask=um_maritime_mask,                          &
                            um_mask_grid_type=um_land_mask_grid_type,          &
                            lon_ref=lon,                                       &
                            lat_ref=lat,                                       &
                            idx_lon_nn=idx_lon_nn,                             &
                            idx_lat_nn=idx_lat_nn)
    maritime_field_adjustments%adjusted_dst_to_src_map_2D(i,1) = idx_lon_nn
    maritime_field_adjustments%adjusted_dst_to_src_map_2D(i,2) = idx_lat_nn
  END DO
  !
  ! Set initialisation flag of maritime field adjustment
  maritime_field_adjustments%initialised = .TRUE.

END IF

! Finalise masks
CALL lfricinp_finalise_masks

END SUBROUTINE um2lfric_init_masked_field_adjustments

END MODULE um2lfric_init_masked_field_adjustments_mod
