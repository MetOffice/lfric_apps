! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE field_list_mod
!
! This module defines and provides access to the internal global field list of
! the API. It also provides routines for initialising the field list and
! creating a pointer to the field in the list.
!

USE field_mod,         ONLY: field_type
USE log_mod,           ONLY: log_event, LOG_LEVEL_ERROR
USE constants_def_mod, ONLY: field_name_len, max_no_fields

IMPLICIT NONE

! Field array containing all fields
TYPE(field_type), TARGET :: field_list(max_no_fields)

! Field write id array
CHARACTER(LEN=field_name_len) :: field_io_name_list(max_no_fields)

! Actual number of fields in field list
INTEGER :: no_fields

CONTAINS

SUBROUTINE init_field_list()
!
! Routine that initialises the global field list
!

IMPLICIT NONE

!
! Local variables
!
! Iterable
INTEGER :: l

! Set number of field to zero
no_fields = 0

! Set associated field array that contains the field dump write names to blank
! strings.
DO l = 1, max_no_fields
  field_io_name_list(l) = REPEAT(' ', field_name_len)
END DO

END SUBROUTINE init_field_list


FUNCTION get_field_pointer(field_id)
!
! Function that produces a pointer to a field with a given field id
!

IMPLICIT NONE

!
! Arguments
!
! Field identifier
CHARACTER(LEN=*), INTENT(IN) :: field_id

!
! Local variables
!
! Field index
INTEGER :: field_index

! pointer to field in the in the field list
TYPE(field_type), POINTER :: get_field_pointer

! Iterable
INTEGER :: l

! Loop over field list items to find field id and store its position index in
! the global field list.
field_index = 0
DO l = 1, no_fields
  IF (TRIM(field_id) == TRIM(field_list(l)%get_name())) THEN
    field_index = l
    EXIT
  END IF
END DO

! Create field pointer to the field in the global field list with the given id.
! If no field in global field list exist with the given id, report issue to
! user and abort the API.
IF (field_index == 0) THEN ! No field in global field list has the given id

  get_field_pointer => NULL()
  CALL log_event('Field ' // TRIM(field_id) // ' is not in global field list', &
                 LOG_LEVEL_ERROR)

 ELSE ! Field with given id has been found in the global field list

  get_field_pointer => field_list(field_index)

END IF

END FUNCTION get_field_pointer


FUNCTION get_field_index(field_id)
!
! Function that produces the index in the global field list a given field_id
! corresponds to.
!

IMPLICIT NONE

!
! Arguments
!
! Field identifier
CHARACTER(LEN=*), INTENT(IN) :: field_id

!
! Local variables
!
! Field index
INTEGER :: get_field_index

! Iterable
INTEGER :: l

! Loop over field list items to find field id and store its position index in
! the global field list.
get_field_index = 0
DO l = 1, no_fields
  IF (TRIM(field_id) == TRIM(field_list(l)%get_name())) THEN
    get_field_index = l
    EXIT
  END IF
END DO

! If no field in global field list exist with the given id, report issue to
! user and abort the API.
IF (get_field_index == 0) THEN ! No field in global field list has the given id

  CALL log_event('Field ' // TRIM(field_id) // ' is not in global field list', &
                 LOG_LEVEL_ERROR)

END IF

END FUNCTION get_field_index

END MODULE field_list_mod
