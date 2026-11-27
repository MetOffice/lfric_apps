! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Some content in this file was generated or refactored with assistance from
! - Met Office Github Copilot Enterprise (Claude Opus 5.5)
module lfric2um_regrid_weights_mod

use lfricinp_regrid_weights_type_mod, only: lfricinp_regrid_weights_type

! Intrinsic modules
use, intrinsic :: iso_fortran_env, only : int32, int64, real64

implicit none

private

! Regridding Weights
type(lfricinp_regrid_weights_type), public, target ::                          &
                         mesh_face_centre_to_grid_p_bilinear,                  &
                         mesh_face_centre_to_grid_p_neareststod,               &
                         mesh_face_centre_to_grid_u_bilinear,                  &
                         mesh_face_centre_to_grid_v_bilinear,                  &
                         mesh_face_centre_to_grid_p_land_mask_bilinear,        &
                         mesh_face_centre_to_grid_p_land_mask_neareststod,     &
                         mesh_face_centre_to_grid_p_sea_mask_bilinear,         &
                         mesh_face_centre_to_grid_p_sea_mask_neareststod

public ::  lfric2um_regrid_weightsfile_ctl, get_weights

contains

!---------------------------------------------------------

subroutine lfric2um_regrid_weightsfile_ctl()
! Description:
!  Control routine to handle weights file reading and some processing
!  of weights to convert to 2D array indices
use lfric2um_namelists_mod, only: lfric2um_config
use lfricinp_um_grid_mod, only: um_grid
implicit none

! P points to face centre interpolation
call mesh_face_centre_to_grid_p_bilinear % load(                               &
     lfric2um_config%weights_file_face_centre_to_p_bilinear)
call mesh_face_centre_to_grid_p_bilinear % populate_dst_address_2D(            &
     int(um_grid % num_p_points_x, kind=int32))

! P points to face centre nearest neighbour
call mesh_face_centre_to_grid_p_neareststod % load(                            &
     lfric2um_config%weights_file_face_centre_to_p_neareststod)
call mesh_face_centre_to_grid_p_neareststod % populate_dst_address_2D(         &
     int(um_grid % num_p_points_x, kind=int32))

! U points to face centre interpolation
call mesh_face_centre_to_grid_u_bilinear % load(                               &
     lfric2um_config%weights_file_face_centre_to_u_bilinear)
call mesh_face_centre_to_grid_u_bilinear % populate_dst_address_2D(            &
     int(um_grid % num_u_points_x, kind=int32))

! V points to face centre interpolation
call mesh_face_centre_to_grid_v_bilinear % load(                               &
     lfric2um_config%weights_file_face_centre_to_v_bilinear)
call mesh_face_centre_to_grid_v_bilinear % populate_dst_address_2D(            &
     int(um_grid % num_v_points_x, kind=int32))

! Optional land-masked P points weights, used for land fields
call load_optional_p_weights(                                                  &
     mesh_face_centre_to_grid_p_land_mask_bilinear,                            &
     lfric2um_config%weights_file_face_centre_to_p_land_mask_bilinear)
call load_optional_p_weights(                                                  &
     mesh_face_centre_to_grid_p_land_mask_neareststod,                         &
     lfric2um_config%weights_file_face_centre_to_p_land_mask_neareststod)

! Optional sea-masked P points weights, used for sea fields
call load_optional_p_weights(                                                  &
     mesh_face_centre_to_grid_p_sea_mask_bilinear,                             &
     lfric2um_config%weights_file_face_centre_to_p_sea_mask_bilinear)
call load_optional_p_weights(                                                  &
     mesh_face_centre_to_grid_p_sea_mask_neareststod,                          &
     lfric2um_config%weights_file_face_centre_to_p_sea_mask_neareststod)

end subroutine lfric2um_regrid_weightsfile_ctl

!---------------------------------------------------------

subroutine load_optional_p_weights(weights, filename)
! Description:
!  Load P points weights from filename, unless it is 'unset', in which
!  case the weights are left unallocated
use lfricinp_um_grid_mod, only: um_grid
implicit none

type(lfricinp_regrid_weights_type), intent(in out) :: weights
character(len=*),                   intent(in)     :: filename

if (trim(filename) == 'unset') return

call weights % load(filename)
call weights % populate_dst_address_2D(                                        &
     int(um_grid % num_p_points_x, kind=int32))

end subroutine load_optional_p_weights

!---------------------------------------------------------

function get_weights(stashcode) result (weights)

! Description:
!  Takes stashcode as argument, interogates stashmaster grid
!  code to determine which grid location the field points sit
!  on and returns pointer to the appropriate weights file

! Intrinsic modules
use, intrinsic :: iso_fortran_env, only : int64

! lfricinputs modules
use lfricinp_stashmaster_mod,    only: get_stashmaster_item, grid,    &
                                       land_compressed, ozone_points, &
                                       p_points, u_points, v_points,  &
                                       p_points_values_over_sea

use lfricinp_regrid_options_mod, only: interp_method, nn_fields, &
                                       specify_nearest_neighbour


! LFRic modules
use log_mod, only: log_event, log_scratch_space, LOG_LEVEL_DEBUG, &
                   LOG_LEVEL_INFO,  LOG_LEVEL_ERROR

implicit none

! Arguments
integer(kind=int64), intent(in) :: stashcode
! Result
type(lfricinp_regrid_weights_type), pointer :: weights

! Local variables
integer(kind=int64) :: horiz_grid_code = 0
type(lfricinp_regrid_weights_type), pointer :: land_mask_weights
type(lfricinp_regrid_weights_type), pointer :: sea_mask_weights
character(len=32) :: method

! Check if interpolation method is supported.
if (trim(interp_method) /= 'bilinear') then
  write(log_scratch_space, '(A)') 'Unsupported interpolation method: '  &
                                // trim(interp_method)
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)
ENDIF

! Get grid type code from STASHmaster entry
horiz_grid_code = get_stashmaster_item(stashcode, grid)

select case(horiz_grid_code)
case( u_points )
  weights => mesh_face_centre_to_grid_u_bilinear
  write(log_scratch_space, '((A,I4,A))')                                &
     "Will use bilinear interpolation for stashcode: ", stashcode, " on grid_u"
  call log_event(log_scratch_space, LOG_LEVEL_DEBUG)

case( v_points )
  weights => mesh_face_centre_to_grid_v_bilinear
  write(log_scratch_space, '((A,I4,A))')                                &
     "Will use bilinear interpolation for stashcode: ", stashcode, " on grid_v"
  call log_event(log_scratch_space, LOG_LEVEL_DEBUG)

case( p_points, ozone_points, land_compressed, p_points_values_over_sea )

  ! Use nearest neighbour weights if requested for this stashcode,
  ! otherwise use the default (bilinear) weights
  if (any(specify_nearest_neighbour(1:nn_fields) == stashcode)) then
    weights => mesh_face_centre_to_grid_p_neareststod
    land_mask_weights => mesh_face_centre_to_grid_p_land_mask_neareststod
    sea_mask_weights => mesh_face_centre_to_grid_p_sea_mask_neareststod
    method = 'nearest neighbour'
  else
    weights => mesh_face_centre_to_grid_p_bilinear
    land_mask_weights => mesh_face_centre_to_grid_p_land_mask_bilinear
    sea_mask_weights => mesh_face_centre_to_grid_p_sea_mask_bilinear
    method = 'bilinear'
  end if

  ! Land and sea fields use the land- or sea-masked version of those
  ! weights, if loaded
  if (horiz_grid_code == land_compressed .and.                          &
      allocated(land_mask_weights % remap_matrix)) then
    weights => land_mask_weights
    method = 'land-masked ' // trim(method)
  else if (horiz_grid_code == p_points_values_over_sea .and.            &
           allocated(sea_mask_weights % remap_matrix)) then
    weights => sea_mask_weights
    method = 'sea-masked ' // trim(method)
  end if

  write(log_scratch_space, '((A,I4,A))') "Will use " // trim(method) // &
     " interpolation for stashcode: ", stashcode, " on grid_p"
  call log_event(log_scratch_space, LOG_LEVEL_DEBUG)

case DEFAULT
  write(log_scratch_space, '(2(A,I0))')                                 &
        "Unsupported horizontal grid type code: ", horiz_grid_code,     &
        " encountered during regrid of stashcode", stashcode
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)

end select

if (.not. allocated(weights%remap_matrix)) then
  call log_event("Attempted to select unallocated weights matrix",      &
                  LOG_LEVEL_ERROR)
end if

end function get_weights

end module lfric2um_regrid_weights_mod
