! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE um2lfric_post_process_fields_mod

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY : real64, int64

IMPLICIT NONE

PRIVATE

PUBLIC :: um2lfric_post_process_fields

CONTAINS

SUBROUTINE um2lfric_post_process_fields(field, stashcode)
! Description: Applies transformations to fields. Should be applied
! post regridding but before the fields are transferred into the
! LFRic field types.

! lfric modules
USE log_mod,         ONLY: log_event, LOG_LEVEL_INFO, log_scratch_space
! um2lfric modules
USE lfricinp_um_grid_mod, ONLY: um_grid
! lfricinputs modules
USE lfricinp_stashmaster_mod, ONLY: stashcode_area_cf, stashcode_v, &
    get_stashmaster_item, pseudL, snow_layers_and_tiles
USE lfricinp_reorder_snow_field_mod, ONLY: lfricinp_reorder_snow_field
USE lfricinp_add_bottom_level_mod,   ONLY: lfricinp_add_bottom_level
USE lfricinp_regrid_options_mod,     ONLY: winds_on_w3
USE lfricinp_um_parameters_mod,      ONLY: um_rmdi

IMPLICIT NONE

! Arguments
REAL(KIND=real64), ALLOCATABLE, INTENT(IN OUT) :: field(:,:)
INTEGER(KIND=int64), INTENT(IN) :: stashcode

! Any post procesing steps determined by STASHcode

SELECT CASE (stashcode)
CASE (stashcode_area_cf)
  WRITE(log_scratch_space, '(A,I0)') "Adding bottom level for stashcode: ",    &
       stashcode
  CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
  CALL lfricinp_add_bottom_level(field)

CASE (stashcode_v)
  IF (.NOT. winds_on_w3) THEN
    WRITE(log_scratch_space,'(A)') "Change sign of 'V' component of " //       &
                                   "wind on W2H function space"
    CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
    ! Switch sign of data, but leave RMDI values unchanged
    field(:,:) = field(:,:) * MERGE(-1.0, 1.0, field /= um_rmdi)
  END IF
END SELECT

! Post processing steps determined by pseudo level code

SELECT CASE (get_stashmaster_item(stashcode, pseudL))
CASE(snow_layers_and_tiles)
  WRITE(log_scratch_space, '(A,I0)') "Reordering snow layers for stashcode: ", &
       stashcode
  CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
  CALL lfricinp_reorder_snow_field(field, um_grid)
END SELECT

END SUBROUTINE um2lfric_post_process_fields

END MODULE um2lfric_post_process_fields_mod
