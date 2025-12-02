! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module dump_generator_mod
!
! This module provides routine that generates the fields in the global field
! list, (i.e fill them with data as per user defined API configuration) and
! write the fields to dump.
!

use field_list_mod,    only: field_list, field_io_name_list, no_fields,        &
                             get_field_index
use log_mod,           only: log_event, log_scratch_space, LOG_LEVEL_INFO

implicit none

contains

subroutine dump_generator()
!
! This routine generates all the fields in the internally stored global field
! list, and write the fields to dump.
!

use constants_def_mod,         only: empty_string, field_name_len
use dependency_graph_list_mod, only: no_dependency_graphs,                     &
                                     dependency_graph_list,                    &
                                     generation_order

implicit none

!
! Local variables
!
! Arrays of input and output field names in a given dependency graph
character(len=field_name_len), allocatable :: input_field_name(:),             &
                                              output_field_name(:)
!
! Number of input and output fields in a given dependency graph
integer :: input_field_no, output_field_no
!
! Index of field in global field list
integer :: global_field_index
!
! Field dependency matrix showing dependencies between fields
integer :: field_dependency_matrix(no_fields, no_fields)
!
! Logical flag array to indicate which fields in global field list has been
! generated
logical :: l_field_generated(no_fields)
!
! Iterable(s)
integer :: i, j, k, l, ii, jj

!
! Set-up the field dependency matrix. This is used exclusively to determine when
! its safe to finalise a field after it has been generated.
!
! First initialise the matrix ...
field_dependency_matrix = 0
!
! ... then set dependencies ...
do l = 1, no_dependency_graphs

  if (allocated(dependency_graph_list(l)%input_field)) then
    do j = 1, size(dependency_graph_list(l)%output_field)

      jj = get_field_index(                                                    &
           dependency_graph_list(l)%output_field(j)%field_ptr%get_name()       &
                        )
      do i = 1, size(dependency_graph_list(l)%input_field)

        ii = get_field_index(                                                  &
             dependency_graph_list(l)%input_field(i)%field_ptr%get_name()      &
                        )
        field_dependency_matrix(ii,jj) = 1

      end do

    end do
  end if

end do
!
! Done setting up field dependency matrix
!

!
! Initialise logical array that flags up if a field in the global field list has
! already been generated
!
l_field_generated(:) = .false.
!
! Done initialising field generation logical array
!

!
! Apply field generators to dependency graphs in correct order. Also write
! fields to dump as soon as they are generated. Any fields which are not
! required in the generation of the remaining fields will be finalised
! immediately.
!
do i = 1, no_dependency_graphs ! Loop over all dependency graphs.

  ! Find index of next dependency graph to process
  k = generation_order(i)

  ! Report to user which generator will be run
  write(log_scratch_space,'(A)') 'About to run generator ' //                  &
                                 trim(dependency_graph_list(k)%gen%identifier)
  call log_event(log_scratch_space, LOG_LEVEL_INFO)

  ! Get input field name list for later use
  if (allocated(dependency_graph_list(k)%input_field)) then

    input_field_no = size(dependency_graph_list(k)%input_field)
    allocate(input_field_name(input_field_no))
    do l = 1, input_field_no
      input_field_name(l) =                                                    &
           trim(dependency_graph_list(k)%input_field(l)%field_ptr%get_name())
    end do

  ENDIF

  ! Get output field name list for later use
  output_field_no = size(dependency_graph_list(k)%output_field)
  allocate(output_field_name(output_field_no))
  do l = 1, output_field_no
    output_field_name(l) =                                                     &
          trim(dependency_graph_list(k)%output_field(l)%field_ptr%get_name())
  end do

  ! Report to user which fields will be generated and what fields they will
  ! generated from.
  do l = 1, output_field_no
    if (allocated(dependency_graph_list(k)%input_field)) then
      write(log_scratch_space,'(100(A,1X))') 'Dependencies for ' //            &
                                              output_field_name(l) // ' are:', &
                                              (input_field_name(j),            &
                                               j = 1, input_field_no)
    else
      write(log_scratch_space,'(A)') 'There are no dependencies for ' //       &
                                      output_field_name(l)
    end if
    call log_event(log_scratch_space, LOG_LEVEL_INFO)
  end do

  ! Generate the fields in output field list of the dependency graph.
  call dependency_graph_list(k)%generate_fields

  ! Report to user that running the generator has completed.
  write(log_scratch_space,'(A)') 'Done running generator ' //                  &
                                 trim(dependency_graph_list(k)%gen%identifier)
  call log_event(log_scratch_space, LOG_LEVEL_INFO)

  ! Write the newly generated fields to dump if needed. Since a field has
  ! now been generated, its dependencies on other fields in the field dependency
  ! can now be removed, and finally also flag up that its been generated.
  do l = 1, output_field_no

    global_field_index =  get_field_index(output_field_name(l))
    if (trim(field_io_name_list(global_field_index)) /= empty_string) then
      call dump_write(global_field_index)
    end if

    l_field_generated(global_field_index) = .true.
    field_dependency_matrix(:,global_field_index) = 0

  end do

  ! Finalise any fields in the global field list which are no longer needed,
  ! i.e. no other remaining fields depend on it and its been generated already.
  do l = 1, no_fields
    if (all(field_dependency_matrix(l,:) == 0) .and. l_field_generated(l)) then
      call field_list(l) % field_final()
    end if
  end do

  ! Deallocate field name lists
  if (allocated(input_field_name)) deallocate(input_field_name)
  if (allocated(output_field_name)) deallocate(output_field_name)

end do ! End of loop of dependency graphs.

end subroutine dump_generator


subroutine dump_write(global_field_index)
!
! This routine writes a field in the global field list, which has their
! write_name parameter set, to the target dump.
!

use lfric_xios_write_mod, only: write_field_generic
use field_parent_mod,     only: write_interface
use field_mod,            only: field_proxy_type

implicit none

!
! Argument(s)
!
! Index of field in global field list to be written to dump
integer :: global_field_index

!
! Local variables
!
! Number of layers in a field
integer :: no_layers
!
! Field proxy
type(field_proxy_type) :: field_proxy
!
! IO procedure pointers
procedure(write_interface), pointer :: tmp_write_ptr

! Define IO procedure pointers for 2D and 3D fields
tmp_write_ptr => write_field_generic

! Report which field is about to be written to dump
write(log_scratch_space,'(A)')                                                 &
     'Writing field ' // trim(field_list(global_field_index)%get_name()) //    &
     ' to target dump'
call log_event(log_scratch_space, LOG_LEVEL_INFO)

! Find number of vertical layers in the field
field_proxy = field_list(global_field_index) % get_proxy()
no_layers = field_proxy % vspace % get_nlayers()

! Set field write behaviour
call field_list(global_field_index) % set_write_behaviour(tmp_write_ptr)

! Write the field to dump
call field_list(global_field_index) %                                          &
           write_field('write_' // trim(field_io_name_list(global_field_index)))

! Nullify IO procedure pointers
nullify(tmp_write_ptr)

end subroutine dump_write

end module dump_generator_mod
