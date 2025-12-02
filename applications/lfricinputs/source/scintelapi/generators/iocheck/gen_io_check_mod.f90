! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module gen_io_check_mod
!
! This module provides a routine to verify the data in a dependency graph is
! consistent with what the dependency graph generator expects.
!

use dependency_graph_mod, only: dependency_graph
use log_mod,              only: log_event, log_scratch_space,                  &
                                LOG_LEVEL_ERROR
use constants_def_mod,    only: genpar_len

implicit none

contains

subroutine gen_io_check(dep_graph, input_field_no, input_field_fs,             &
                        output_field_no, output_field_fs, parameter_no)
!
! This is a service routine to be called from within a generator. It verifies
! whether the data (field definitions and parameter list) in the dependency
! graph is consistent with what the dependency graph generator expects.
! The actual set of checks performed is determined by the arguments supplied
! when calling this routine from within the generator.
!

implicit none

!
! Argument definitions:
!
! Dependency graph to be processed
class(dependency_graph), intent(in) :: dep_graph
!
! Number of in/out fields expected
integer, optional, intent(in) :: input_field_no, output_field_no
!
! Input and output lists of function spaces
integer, optional, intent(in) :: input_field_fs(:), output_field_fs(:)
!
! Number of parameters expected in parameter list
integer, optional, intent(in) :: parameter_no

!
! Local variables
!
! Actual number of input fields
integer :: actual_number_of_input_fields
!
! Actual number of input fields
integer :: actual_number_of_output_fields
!
! Dummy function space variable
integer :: fspace
!
! Logical for checking if there is a function space mismatch
logical :: l_func_space_mismatch
!
! Local copy of parameter list
character(len=genpar_len) :: parlist
!
! Logical flag to check if currently scanning a parameter in parmeter list
logical :: l_scanning_parameter
!
! Actual number of parameters
integer :: actual_number_of_parameters
!
! Iterable(s)
integer :: i

!
! Find actual number of input and output fields
!
! Find number of input fields in the dependency graph
if (allocated(dep_graph % input_field)) then
  actual_number_of_input_fields = size(dep_graph % input_field)
else
  actual_number_of_input_fields = 0
end if

! Find number of output fields in the dependency graph
if (allocated(dep_graph % output_field)) then
  actual_number_of_output_fields = size(dep_graph % output_field)
else
  actual_number_of_output_fields = 0
end if

!
! Checks on number of input and output fields
!
! Check number of input fields in dependency graph to the specified number. If
! they fail to compare report the issue and abort.
!
if (present(input_field_no)) then

  if (actual_number_of_input_fields /= input_field_no) then
    write(log_scratch_space,'(A)') 'GEN_IO_CHECK: Incorrect number of ' //     &
                                   'input fields'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
  end if

end if

! Check number of output fields in dependency graph to the specified number. If
! they fail to compare report the issue and abort.
if (present(output_field_no)) then

  if (actual_number_of_output_fields /= output_field_no) then
    write(log_scratch_space,'(A)') 'GEN_IO_CHECK: Incorrect number of ' //     &
                                   'output fields'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
  end if

end if
!
! Done checking number of input and output fields
!

!
! Compare actual field function spaces with the user supplied lists
!
! First compare sizes
if (present(input_field_fs)) then

  if (size(input_field_fs) /= actual_number_of_input_fields) then
    write(log_scratch_space,'(A)') 'GEN_IO_CHECK: Incorrect number of ' //     &
                                   'input function spaces provided!'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
  end if

end if
!
if (present(output_field_fs)) then

  if (size(output_field_fs) /= actual_number_of_output_fields) then
    write(log_scratch_space,'(A)') 'GEN_IO_CHECK: Incorrect number of ' //     &
                                   'output function spaces provided!'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
  end if

end if
!
! Now compare function spaces
l_func_space_mismatch = .false.
!
if (present(input_field_fs)) then

  do i = 1, actual_number_of_input_fields
    fspace = dep_graph % input_field(i) % field_ptr % which_function_space()
    if (input_field_fs(i) /= fspace) l_func_space_mismatch = .true.
  end do

end if
!
if (present(output_field_fs)) then

  do i = 1, actual_number_of_output_fields
    fspace = dep_graph % output_field(i) % field_ptr % which_function_space()
    if (output_field_fs(i) /= fspace) l_func_space_mismatch = .true.
  end do

end if
!
if (l_func_space_mismatch) then

  write(log_scratch_space,'(A)') 'GEN_IO_CHECK: Function space mismatch '//    &
                                 'found. Check your field definitions and ' // &
                                 'generator used are commensurate.'
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)

end if
!
! Done checking function spaces
!

!
! Check if number of parameters in parameter list is consisent with whats
! expected
!
if (present(parameter_no)) then

  ! Iterate over generator parameter list character by character to determine
  ! number of parameters.

  parlist = adjustl(dep_graph % genpar)
  actual_number_of_parameters = 0
  l_scanning_parameter = .false.

  do i = 1, len(trim(parlist))             ! Loop over each character

    if(parlist(i:i) /= ' ') then           ! Non-blank character detected ...

      if (.not. l_scanning_parameter) then ! ... not currently scanning a
                                           ! parameter, so must a new parameter
        actual_number_of_parameters = actual_number_of_parameters + 1
        l_scanning_parameter = .true.      ! ... scanning a parameter, so set to
                                           ! true

      end if

    else                                   ! Blank character detected, so not
                                           ! scanning a current parameter in
      l_scanning_parameter = .false.       ! list

    end if

  end do

  if (actual_number_of_parameters /= parameter_no) then

    write(log_scratch_space,'(A)') 'GEN_IO_CHECK: The expected ' //            &
                                   'and actual number of parameters do ' //    &
                                   'not match!'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)

  end if

end if
!
! Done checking number of parameters
!

end subroutine gen_io_check

end module gen_io_check_mod
