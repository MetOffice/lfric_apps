! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lfric2um_regrid_weights_mod

USE lfricinp_regrid_weights_type_mod, ONLY: lfricinp_regrid_weights_type

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY : int32, int64, real64

IMPLICIT NONE

PRIVATE

! Regridding Weights
TYPE(lfricinp_regrid_weights_type), PUBLIC, TARGET ::                          &
                                    mesh_face_centre_to_grid_p_bilinear,       &
                                    mesh_face_centre_to_grid_u_bilinear,       &
                                    mesh_face_centre_to_grid_v_bilinear

PUBLIC ::  lfric2um_regrid_weightsfile_ctl, get_weights

CONTAINS

!---------------------------------------------------------

SUBROUTINE lfric2um_regrid_weightsfile_ctl()
! Description:
!  Control routine to handle weights file reading and some processing
!  of weights to convert to 2D array indices
USE lfric2um_namelists_mod, ONLY: lfric2um_config
USE lfricinp_um_grid_mod, ONLY: um_grid
IMPLICIT NONE

! P points to face centre interpolation
CALL mesh_face_centre_to_grid_p_bilinear % load(                               &
     lfric2um_config%weights_file_face_centre_to_p_bilinear)
CALL mesh_face_centre_to_grid_p_bilinear % populate_dst_address_2D(            &
     INT(um_grid % num_p_points_x, KIND=int32))

! U points to face centre interpolation
CALL mesh_face_centre_to_grid_u_bilinear % load(                               &
     lfric2um_config%weights_file_face_centre_to_u_bilinear)
CALL mesh_face_centre_to_grid_u_bilinear % populate_dst_address_2D(            &
     INT(um_grid % num_u_points_x, KIND=int32))

! V points to face centre interpolation
CALL mesh_face_centre_to_grid_v_bilinear % load(                               &
     lfric2um_config%weights_file_face_centre_to_v_bilinear)
CALL mesh_face_centre_to_grid_v_bilinear % populate_dst_address_2D(            &
     INT(um_grid % num_v_points_x, KIND=int32))

END SUBROUTINE lfric2um_regrid_weightsfile_ctl

!---------------------------------------------------------

FUNCTION get_weights(stashcode) RESULT (weights)

! Description:
!  Takes stashcode as argument, interogates stashmaster grid
!  code to determine which grid location the field points sit
!  on and returns pointer to the appropriate weights file

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY : int64
! lfricinputs modules
USE lfricinp_stashmaster_mod, ONLY: get_stashmaster_item, p_points, &
   u_points, v_points, ozone_points, grid
USE lfricinp_regrid_options_mod, ONLY: interp_method

! LFRic modules
USE log_mod, ONLY: log_event, log_scratch_space, LOG_LEVEL_ERROR

IMPLICIT NONE

! Arguments
INTEGER(KIND=int64), INTENT(IN) :: stashcode
! Result
TYPE(lfricinp_regrid_weights_type), POINTER :: weights

! Local variables
INTEGER(KIND=int64) :: horiz_grid_code = 0

! Check if interpolation method is supported.
IF (TRIM(interp_method) /= 'bilinear') THEN
  WRITE(log_scratch_space, '(A)') 'Unsupported interpolation method: '  &
                                // TRIM(interp_method)
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
ENDIF

! Get grid type code from STASHmaster entry
horiz_grid_code = get_stashmaster_item(stashcode, grid)

SELECT CASE(horiz_grid_code)
CASE( u_points )
  weights => mesh_face_centre_to_grid_u_bilinear
CASE( v_points )
  weights => mesh_face_centre_to_grid_v_bilinear
CASE( p_points, ozone_points )
   weights => mesh_face_centre_to_grid_p_bilinear
CASE DEFAULT
  WRITE(log_scratch_space, '(2(A,I0))')                                        &
        "Unsupported horizontal grid type code: ",                             &
        horiz_grid_code, " encountered during regrid of stashcode", stashcode
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END SELECT

IF (.NOT. ALLOCATED(weights%remap_matrix)) THEN
  CALL log_event("Attempted to select unallocated weights matrix",             &
                  LOG_LEVEL_ERROR)
END IF

END FUNCTION get_weights

END MODULE lfric2um_regrid_weights_mod
