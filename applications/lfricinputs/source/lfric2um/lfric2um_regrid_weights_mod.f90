! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module lfric2um_regrid_weights_mod

use lfricinp_regrid_weights_type_mod, only: lfricinp_regrid_weights_type

! Intrinsic modules
use, intrinsic :: iso_fortran_env, only : int32, int64, real64

implicit none

private

! Regridding Weights
type(lfricinp_regrid_weights_type), public, target ::                          &
                                    mesh_face_centre_to_grid_p_bilinear,       &
                                    mesh_face_centre_to_grid_u_bilinear,       &
                                    mesh_face_centre_to_grid_v_bilinear

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

end subroutine lfric2um_regrid_weightsfile_ctl

!---------------------------------------------------------

function get_weights(stashcode) result (weights)

! Description:
!  Takes stashcode as argument, interogates stashmaster grid
!  code to determine which grid location the field points sit
!  on and returns pointer to the appropriate weights file

! Intrinsic modules
use, intrinsic :: iso_fortran_env, only : int64
! lfricinputs modules
use lfricinp_stashmaster_mod, only: get_stashmaster_item, p_points, &
                                    u_points, v_points, ozone_points, grid
use lfricinp_regrid_options_mod, only: interp_method

! LFRic modules
use log_mod, only: log_event, log_scratch_space, LOG_LEVEL_ERROR

implicit none

! Arguments
integer(kind=int64), intent(in) :: stashcode
! Result
type(lfricinp_regrid_weights_type), pointer :: weights

! Local variables
integer(kind=int64) :: horiz_grid_code = 0

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
case( v_points )
  weights => mesh_face_centre_to_grid_v_bilinear
case( p_points, ozone_points )
   weights => mesh_face_centre_to_grid_p_bilinear
case DEFAULT
  write(log_scratch_space, '(2(A,I0))')                                        &
        "Unsupported horizontal grid type code: ",                             &
        horiz_grid_code, " encountered during regrid of stashcode", stashcode
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)
end select

if (.not. allocated(weights%remap_matrix)) then
  call log_event("Attempted to select unallocated weights matrix",             &
                  LOG_LEVEL_ERROR)
end if

end function get_weights

end module lfric2um_regrid_weights_mod
