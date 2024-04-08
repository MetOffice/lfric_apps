! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE dump_generator_mod
!
! This module provides routine that generates the fields in the global field
! list, (i.e fill them with data as per user defined API configuration) and
! write the fields to dump.
!

USE field_list_mod,    ONLY: field_list, field_io_name_list, no_fields,        &
                             get_field_index
USE log_mod,           ONLY: log_event, log_scratch_space, LOG_LEVEL_INFO

IMPLICIT NONE

CONTAINS

SUBROUTINE dump_generator()
!
! This routine generates all the fields in the internally stored global field
! list, and write the fields to dump.
!

USE constants_def_mod,         ONLY: empty_string, field_name_len
USE dependency_graph_list_mod, ONLY: no_dependency_graphs,                     &
                                     dependency_graph_list,                    &
                                     generation_order

IMPLICIT NONE

!
! Local variables
!
! Arrays of input and output field names in a given dependency graph
CHARACTER(LEN=field_name_len), ALLOCATABLE :: input_field_name(:),             &
                                              output_field_name(:)
!
! Number of input and output fields in a given dependency graph
INTEGER :: input_field_no, output_field_no
!
! Index of field in global field list
INTEGER :: global_field_index
!
! Field dependency matrix showing dependencies between fields
INTEGER :: field_dependency_matrix(no_fields, no_fields)
!
! Logical flag array to indicate which fields in global field list has been
! generated
LOGICAL :: l_field_generated(no_fields)
!
! Iterable(s)
INTEGER :: i, j, k, l, ii, jj

!
! Set-up the field dependency matrix. This is used exclusively to determine when
! its safe to finalise a field after it has been generated.
!
! First initialise the matrix ...
field_dependency_matrix = 0
!
! ... then set dependencies ...
DO l = 1, no_dependency_graphs

  IF (ALLOCATED(dependency_graph_list(l)%input_field)) THEN
    DO j = 1, SIZE(dependency_graph_list(l)%output_field)

      jj = get_field_index(                                                    &
           dependency_graph_list(l)%output_field(j)%field_ptr%get_name()       &
                        )
      DO i = 1, SIZE(dependency_graph_list(l)%input_field)

        ii = get_field_index(                                                  &
             dependency_graph_list(l)%input_field(i)%field_ptr%get_name()      &
                        )
        field_dependency_matrix(ii,jj) = 1

      END DO

    END DO
  END IF

END DO
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
DO i = 1, no_dependency_graphs ! Loop over all dependency graphs.

  ! Find index of next dependency graph to process
  k = generation_order(i)

  ! Report to user which generator will be run
  WRITE(log_scratch_space,'(A)') 'About to run generator ' //                  &
                                 TRIM(dependency_graph_list(k)%gen%identifier)
  CALL log_event(log_scratch_space, LOG_LEVEL_INFO)

  ! Get input field name list for later use
  IF (ALLOCATED(dependency_graph_list(k)%input_field)) THEN

    input_field_no = SIZE(dependency_graph_list(k)%input_field)
    ALLOCATE(input_field_name(input_field_no))
    DO l = 1, input_field_no
      input_field_name(l) =                                                    &
           TRIM(dependency_graph_list(k)%input_field(l)%field_ptr%get_name())
    END DO

  ENDIF

  ! Get output field name list for later use
  output_field_no = SIZE(dependency_graph_list(k)%output_field)
  ALLOCATE(output_field_name(output_field_no))
  DO l = 1, output_field_no
    output_field_name(l) =                                                     &
          TRIM(dependency_graph_list(k)%output_field(l)%field_ptr%get_name())
  END DO

  ! Report to user which fields will be generated and what fields they will
  ! generated from.
  DO l = 1, output_field_no
    IF (ALLOCATED(dependency_graph_list(k)%input_field)) THEN
      WRITE(log_scratch_space,'(100(A,1X))') 'Dependencies for ' //            &
                                              output_field_name(l) // ' are:', &
                                              (input_field_name(j),            &
                                               j = 1, input_field_no)
    ELSE
      WRITE(log_scratch_space,'(A)') 'There are no dependencies for ' //       &
                                      output_field_name(l)
    END IF
    CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
  END DO

  ! Generate the fields in output field list of the dependency graph.
  CALL dependency_graph_list(k)%generate_fields

  ! Report to user that running the generator has completed.
  WRITE(log_scratch_space,'(A)') 'Done running generator ' //                  &
                                 TRIM(dependency_graph_list(k)%gen%identifier)
  CALL log_event(log_scratch_space, LOG_LEVEL_INFO)

  ! Write the newly generated fields to dump if needed. Since a field has
  ! now been generated, its dependencies on other fields in the field dependency
  ! can now be removed, and finally also flag up that its been generated.
  DO l = 1, output_field_no

    global_field_index =  get_field_index(output_field_name(l))
    IF (TRIM(field_io_name_list(global_field_index)) /= empty_string) THEN
      CALL dump_write(global_field_index)
    END IF

    l_field_generated(global_field_index) = .true.
    field_dependency_matrix(:,global_field_index) = 0

  END DO

  ! Finalise any fields in the global field list which are no longer needed,
  ! i.e. no other remaining fields depend on it and its been generated already.
  DO l = 1, no_fields
    IF (ALL(field_dependency_matrix(l,:) == 0) .AND. l_field_generated(l)) THEN
      CALL field_list(l) % field_final()
    END IF
  END DO

  ! Deallocate field name lists
  IF (ALLOCATED(input_field_name)) DEALLOCATE(input_field_name)
  IF (ALLOCATED(output_field_name)) DEALLOCATE(output_field_name)

END DO ! End of loop of dependency graphs.

END SUBROUTINE dump_generator


SUBROUTINE dump_write(global_field_index)
!
! This routine writes a field in the global field list, which has their
! write_name parameter set, to the target dump.
!

USE lfric_xios_write_mod, ONLY: write_field_generic
USE field_parent_mod,     ONLY: write_interface
USE field_mod,            ONLY: field_proxy_type

IMPLICIT NONE

!
! Argument(s)
!
! Index of field in global field list to be written to dump
INTEGER :: global_field_index

!
! Local variables
!
! Number of layers in a field
INTEGER :: no_layers
!
! Field proxy
TYPE(field_proxy_type) :: field_proxy
!
! IO procedure pointers
PROCEDURE(write_interface), POINTER :: tmp_write_ptr

! Define IO procedure pointers for 2D and 3D fields
tmp_write_ptr => write_field_generic

! Report which field is about to be written to dump
WRITE(log_scratch_space,'(A)')                                                 &
     'Writing field ' // TRIM(field_list(global_field_index)%get_name()) //    &
     ' to target dump'
CALL log_event(log_scratch_space, LOG_LEVEL_INFO)

! Find number of vertical layers in the field
field_proxy = field_list(global_field_index) % get_proxy()
no_layers = field_proxy % vspace % get_nlayers()

! Set field write behaviour
CALL field_list(global_field_index) % set_write_behaviour(tmp_write_ptr)

! Write the field to dump
CALL field_list(global_field_index) %                                          &
           write_field('write_' // TRIM(field_io_name_list(global_field_index)))

! Nullify IO procedure pointers
NULLIFY(tmp_write_ptr)

END SUBROUTINE dump_write

END MODULE dump_generator_mod
