! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lfricinp_check_shumlib_status_mod

IMPLICIT NONE

PRIVATE
PUBLIC :: shumlib

CONTAINS

SUBROUTINE shumlib(routinename, status,                                        &
                   print_on_success, ignore_warning, errorstatus)

! Each routine in the shumlib fieldsfile API returns a status object.
! This routine checks the status of such an object. It can act
! as a wrapper routine if the shumlib function is called directly
! in the argument list, for example:
!
! CALL shumlib("my_shumlib_func", my_shumlib_func(var1, var2))
!
! Here my_shumlib_func is called first, the arguments var1 and
! var2 are passed to the routine. They remain in scope in the calling
! routine and are available to be set if their intent allows.
! my_shumlib_func then returns a status object that is passed into this
! routine for checking.

! Intrinsic modules
USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_BOOL

! LFRic modules
USE log_mod, ONLY: log_event, LOG_LEVEL_ERROR, LOG_LEVEL_INFO

! Shumlib modules
USE f_shum_ff_status_mod, ONLY: shum_ff_status_type, OPERATOR(==),            &
                                OPERATOR(<), SHUMLIB_SUCCESS

IMPLICIT NONE

! Arguments
TYPE(shum_ff_status_type), INTENT(IN) :: status
CHARACTER(LEN=*),          INTENT(IN) :: routinename
LOGICAL(KIND=C_BOOL), OPTIONAL,  INTENT(IN) :: print_on_success, ignore_warning
INTEGER,              OPTIONAL, INTENT(OUT) :: errorstatus
! Internal variables

! Message - set to be the maximum SHUMlib message length (1024) plus a
! reasonable size for a routine name (128)
CHARACTER(LEN=1152) :: message
! Set message
WRITE(message, '(A,A,A,A)') '[', TRIM(routinename), '] ',TRIM(status%message)

IF (status < SHUMLIB_SUCCESS) THEN
  ! This is potentially a warning
  IF (PRESENT(errorstatus)) errorstatus = -1
  IF (PRESENT(ignore_warning)) THEN
    IF (ignore_warning) THEN
      ! Ignore warning is true, print and carry on
      CALL log_event(message, LOG_LEVEL_INFO)
    ELSE
      ! Ignore warning is false, print and abort
      CALL log_event(message, LOG_LEVEL_ERROR)
    END IF
  ELSE
    ! Ignore warning is not set, print and abort
    CALL log_event(message, LOG_LEVEL_ERROR)
  END IF
ELSE IF (status == SHUMLIB_SUCCESS) THEN
  ! This is a success; if print_on_success is provided, and true, print message
  ! otherwise carry on as normal
  IF (PRESENT(errorstatus)) errorstatus = 0
  IF (PRESENT(print_on_success)) THEN
    IF (print_on_success) THEN
      CALL log_event(message, LOG_LEVEL_INFO)
    END IF
  END IF
ELSE
  IF (PRESENT(errorstatus)) errorstatus = 1
  ! This is a definite failure, abort whatever else is provided
  CALL log_event(message, LOG_LEVEL_ERROR)
END IF

END SUBROUTINE shumlib

END MODULE lfricinp_check_shumlib_status_mod
