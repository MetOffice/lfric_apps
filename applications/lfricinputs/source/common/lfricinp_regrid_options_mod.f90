! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lfricinp_regrid_options_mod

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY : int64

! UM2LFRic modules
USE lfricinp_um_parameters_mod,     ONLY: fnamelen, um_imdi

IMPLICIT NONE
PRIVATE
PUBLIC :: lfricinp_init_regrid_options, interp_method, winds_on_w3, &
          regrid_type, specify_nearest_neighbour, nn_fields

CHARACTER(LEN=fnamelen) :: interp_method = 'bilinear'
CHARACTER(LEN=fnamelen) :: regrid_type = 'global_to_global'
LOGICAL :: winds_on_w3 = .TRUE.

INTEGER(KIND=int64), PARAMETER :: max_stash_list = 999
INTEGER(KIND=int64)     :: specify_nearest_neighbour(max_stash_list)
INTEGER :: nn_fields ! Number of nearest neighbour interpolation fields

CONTAINS

SUBROUTINE lfricinp_init_regrid_options(fname)

USE lfricinp_unit_handler_mod, ONLY: get_free_unit
USE log_mod,                   ONLY: log_event, log_scratch_space,             &
                                     LOG_LEVEL_ERROR

IMPLICIT NONE

CHARACTER(LEN=fnamelen) :: fname
INTEGER                 :: status = -1
CHARACTER(LEN=512)      :: message = 'No namelist read'
INTEGER                 :: unit_number
INTEGER                 :: i_stash

NAMELIST /regrid_options/ interp_method, regrid_type, winds_on_w3, &
     specify_nearest_neighbour

specify_nearest_neighbour(:) = um_imdi

CALL get_free_unit(unit_number)

OPEN(UNIT=unit_number, FILE=fname, IOSTAT=status, IOMSG=message)
IF (status /= 0) CALL log_event(message, LOG_LEVEL_ERROR)

READ(unit_number, NML=regrid_options, IOSTAT=status, IOMSG=message)
IF (status /= 0) CALL log_event(message, LOG_LEVEL_ERROR)

IF ( (TRIM(interp_method) /= 'bilinear') .AND.                                 &
     (TRIM(interp_method) /= 'copy') ) THEN
  WRITE(log_scratch_space, '(A)') 'Interpolation method: '                     &
                                  // TRIM(interp_method) // ' not recognised.'
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END IF

IF ( (TRIM(regrid_type) /= 'global_to_global') .AND.                           &
     (TRIM(regrid_type) /= 'lam_to_lam') .AND.                                 &
     (TRIM(regrid_type) /= 'lbc_to_lbc') ) THEN
  WRITE(log_scratch_space, '(A)') 'Regridding type: ' // TRIM(regrid_type) //  &
                                  ' not recognised.'
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END IF

IF ( (TRIM(interp_method) == 'copy') .AND.                                     &
     (TRIM(regrid_type) /= 'lam_to_lam' .AND.                                  &
      TRIM(regrid_type) /= 'lbc_to_lbc') ) THEN
  WRITE(log_scratch_space, '(A)') 'The copy interpolation method is only ' //  &
                                  'valid for LAM to LAM or LBC to LBC ' //     &
                                  'regridding'
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END IF

IF ( (TRIM(interp_method) == 'bilinear') .AND.                                 &
     (.NOT. winds_on_w3) ) THEN
  WRITE(log_scratch_space, '(A)') 'For bilinear interpolation method winds ' //&
                                  'must be placed on W3'
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END IF

IF ( (TRIM(interp_method) == 'copy') .AND. winds_on_w3 ) THEN
  WRITE(log_scratch_space, '(A)')                                              &
                            'For copying method winds must not be placed on W3'
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END IF

nn_fields=0
! Count how many fields in specified lists
DO i_stash = 1, max_stash_list
  IF (specify_nearest_neighbour(i_stash) == um_imdi) THEN
    EXIT
  ELSE
    nn_fields = nn_fields + 1
  END IF
END DO

CLOSE(unit_number)

END SUBROUTINE lfricinp_init_regrid_options

END MODULE lfricinp_regrid_options_mod
