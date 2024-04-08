! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE gen_io_check_mod
!
! This module provides a routine to verify the data in a dependency graph is
! consistent with what the dependency graph generator expects.
!

USE dependency_graph_mod, ONLY: dependency_graph
USE log_mod,              ONLY: log_event, log_scratch_space,                  &
                                LOG_LEVEL_ERROR
USE constants_def_mod,    ONLY: genpar_len

IMPLICIT NONE

CONTAINS

SUBROUTINE gen_io_check(dep_graph, input_field_no, input_field_fs,             &
                        output_field_no, output_field_fs, parameter_no)
!
! This is a service routine to be called from within a generator. It verifies
! whether the data (field definitions and parameter list) in the dependency
! graph is consistent with what the dependency graph generator expects.
! The actual set of checks performed is determined by the arguments supplied
! when calling this routine from within the generator.
!

IMPLICIT NONE

!
! Argument definitions:
!
! Dependency graph to be processed
CLASS(dependency_graph), INTENT(IN) :: dep_graph
!
! Number of in/out fields expected
INTEGER, OPTIONAL, INTENT(IN) :: input_field_no, output_field_no
!
! Input and output lists of function spaces
INTEGER, OPTIONAL, INTENT(IN) :: input_field_fs(:), output_field_fs(:)
!
! Number of parameters expected in parameter list
INTEGER, OPTIONAL, INTENT(IN) :: parameter_no

!
! Local variables
!
! Actual number of input fields
INTEGER :: actual_number_of_input_fields
!
! Actual number of input fields
INTEGER :: actual_number_of_output_fields
!
! Dummy function space variable
INTEGER :: fspace
!
! Logical for checking if there is a function space mismatch
LOGICAL :: l_func_space_mismatch
!
! Local copy of parameter list
CHARACTER(LEN=genpar_len) :: parlist
!
! Logical flag to check if currently scanning a parameter in parmeter list
LOGICAL :: l_scanning_parameter
!
! Actual number of parameters
INTEGER :: actual_number_of_parameters
!
! Iterable(s)
INTEGER :: i

!
! Find actual number of input and output fields
!
! Find number of input fields in the dependency graph
IF (ALLOCATED(dep_graph % input_field)) THEN
  actual_number_of_input_fields = SIZE(dep_graph % input_field)
ELSE
  actual_number_of_input_fields = 0
END IF

! Find number of output fields in the dependency graph
IF (ALLOCATED(dep_graph % output_field)) THEN
  actual_number_of_output_fields = SIZE(dep_graph % output_field)
ELSE
  actual_number_of_output_fields = 0
END IF

!
! Checks on number of input and output fields
!
! Check number of input fields in dependency graph to the specified number. If
! they fail to compare report the issue and abort.
!
IF (PRESENT(input_field_no)) THEN

  IF (actual_number_of_input_fields /= input_field_no) THEN
    WRITE(log_scratch_space,'(A)') 'GEN_IO_CHECK: Incorrect number of ' //     &
                                   'input fields'
    CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
  END IF

END IF

! Check number of output fields in dependency graph to the specified number. If
! they fail to compare report the issue and abort.
IF (PRESENT(output_field_no)) THEN

  IF (actual_number_of_output_fields /= output_field_no) THEN
    WRITE(log_scratch_space,'(A)') 'GEN_IO_CHECK: Incorrect number of ' //     &
                                   'output fields'
    CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
  END IF

END IF
!
! Done checking number of input and output fields
!

!
! Compare actual field function spaces with the user supplied lists
!
! First compare sizes
IF (PRESENT(input_field_fs)) THEN

  IF (SIZE(input_field_fs) /= actual_number_of_input_fields) THEN
    WRITE(log_scratch_space,'(A)') 'GEN_IO_CHECK: Incorrect number of ' //     &
                                   'input function spaces provided!'
    CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
  END IF

END IF
!
IF (PRESENT(output_field_fs)) THEN

  IF (SIZE(output_field_fs) /= actual_number_of_output_fields) THEN
    WRITE(log_scratch_space,'(A)') 'GEN_IO_CHECK: Incorrect number of ' //     &
                                   'output function spaces provided!'
    CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
  END IF

END IF
!
! Now compare function spaces
l_func_space_mismatch = .false.
!
IF (PRESENT(input_field_fs)) THEN

  DO i = 1, actual_number_of_input_fields
    fspace = dep_graph % input_field(i) % field_ptr % which_function_space()
    IF (input_field_fs(i) /= fspace) l_func_space_mismatch = .true.
  END DO

END IF
!
IF (PRESENT(output_field_fs)) THEN

  DO i = 1, actual_number_of_output_fields
    fspace = dep_graph % output_field(i) % field_ptr % which_function_space()
    IF (output_field_fs(i) /= fspace) l_func_space_mismatch = .true.
  END DO

END IF
!
IF (l_func_space_mismatch) THEN

  WRITE(log_scratch_space,'(A)') 'GEN_IO_CHECK: Function space mismatch '//    &
                                 'found. Check your field definitions and ' // &
                                 'generator used are commensurate.'
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)

END IF
!
! Done checking function spaces
!

!
! Check if number of parameters in parameter list is consisent with whats
! expected
!
IF (PRESENT(parameter_no)) THEN

  ! Iterate over generator parameter list character by character to determine
  ! number of parameters.

  parlist = ADJUSTL(dep_graph % genpar)
  actual_number_of_parameters = 0
  l_scanning_parameter = .false.

  DO i = 1, LEN(TRIM(parlist))             ! Loop over each character

    IF(parlist(i:i) /= ' ') THEN           ! Non-blank character detected ...

      IF (.NOT. l_scanning_parameter) THEN ! ... not currently scanning a
                                           ! parameter, so must a new parameter
        actual_number_of_parameters = actual_number_of_parameters + 1
        l_scanning_parameter = .true.      ! ... scanning a parameter, so set to
                                           ! true

      END IF

    ELSE                                   ! Blank character detected, so not
                                           ! scanning a current parameter in
      l_scanning_parameter = .false.       ! list

    END IF

  END DO

  IF (actual_number_of_parameters /= parameter_no) THEN

    WRITE(log_scratch_space,'(A)') 'GEN_IO_CHECK: The expected ' //            &
                                   'and actual number of parameters do ' //    &
                                   'not match!'
    CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)

  END IF

END IF
!
! Done checking number of parameters
!

END SUBROUTINE gen_io_check

END MODULE gen_io_check_mod
