! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE dependency_analyser_mod
!
! This module provides routines for performing an analysis of the
! dependencies between dependency graphs, and based on this analysis is find the
! correct order is which the output fields in the dependency graphs should be
! generated.
!
USE log_mod, ONLY: log_event, LOG_LEVEL_ERROR, LOG_LEVEL_INFO,                 &
                   log_scratch_space

IMPLICIT NONE

CONTAINS

SUBROUTINE dependency_analyser()

! This subroutine performs a number tasks and checks to validate the output
! fields, as contained in a list of dependency graphs, can be successfully
! generated:
!
! 1) It checks no dependency graph object contains the same field in its input
! and output lists.
!
! 2) It checks that the output fields across dependency graphs are unique, i.e.
! each output field should belong to one, and only one, dependency graph.
!
! 3) It checks that each dependency graph input field, has a defined method for
! generation, i.e. it in the list of all dependency graph output fields.
!
! 4) It constructs a dependency matrix expressing the dependencies between
! dependency graphs.
!
! 5) Using the information in the dependency matrix, it attemps to construct
! the order in which the dependency graph generators should be run.

USE dependency_graph_list_mod, ONLY: dependency_graph_list,                    &
                                     no_dependency_graphs,                     &
                                     generation_order
USE constants_mod,             ONLY: str_def

IMPLICIT NONE

!
! Local variables:
!
! Iterable
INTEGER :: i, j, l, n

! No of output fields across all dependency graphs
INTEGER :: no_output_fields

! Array of output field names in all dependency graphs
CHARACTER(LEN=str_def), ALLOCATABLE :: out_field_ids(:)

! Dependency graph matrix
INTEGER :: dep_graph_matrix(no_dependency_graphs, no_dependency_graphs)

! Dummy field name variables
CHARACTER(LEN=str_def) :: fid, ofid, ifid

! Logical flag to check dependency between two dependency graphs
LOGICAL :: dependency

! Logical flag to check input and output fields are distinct for any given
! dependency graph
LOGICAL :: l_inout_fields_not_unique

! Logical flag to check for output overlap across dependency graphs
LOGICAL :: l_output_not_unique

! Logical flag to check input fields has a definied method for generation
LOGICAL :: l_input_field_not_generated

CALL log_event('Running DEPENDENCY ANALYSER', LOG_LEVEL_INFO)

!
! Perform checks 1, 2 and 3.
!
! 1) Check the input fields and output fields for each dependency graph are
! distinct, i.e. the lists has no field ids in common.
!
l_inout_fields_not_unique = .false. ! Assume input and output list are unique.

DO l = 1, no_dependency_graphs ! Loop over all dependency graphs

  IF (ALLOCATED(dependency_graph_list(l)%input_field)) THEN ! Check input field
                                                            ! list exist

    DO i = 1, SIZE(dependency_graph_list(l)%input_field) ! Loop over input
                                                         !fields

      ! Set input field name for later compararison
      ifid = dependency_graph_list(l)%input_field(i)%field_ptr%get_name()

      DO j = 1, SIZE(dependency_graph_list(l)%output_field) ! Loop over ouput
                                                            ! field list
        ! Set output field name for comparison
        ofid = dependency_graph_list(l)%output_field(j)%field_ptr%get_name()

        ! Check if input and output field list names are the same. If they are
        ! then report the issue and the offending field name to user and set
        ! l_inout_fields_not_unique = .true.
        IF (TRIM(ifid) == TRIM(ofid)) THEN

          l_inout_fields_not_unique = .true.
          WRITE(log_scratch_space, '(A)') 'Field ' // TRIM(ifid) //            &
                                   ' occurs in both the input and output' //   &
                                   ' lists of a dependency graph'
          CALL log_event(log_scratch_space, LOG_LEVEL_INFO)

        END IF

      END DO ! End of loop over output fields

    END DO ! End of loop over input fields

  END IF ! End of check if input field list exist

END DO ! End of loop over dependency graphs

! If in any of the dependency graphs the input and output fields contain the
! same field, report this to the user and abort the API.
IF (l_inout_fields_not_unique) THEN
  WRITE(log_scratch_space, '(A)') 'At least one of the dependency graphs ' //  &
                                  'has at least one field common to the '  //  &
                                  'its input and output field lists.'
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END IF

! Interlude: Create list of output field names. This list will come in handy in
! further checks later on.
!
! Find total number of output fields accross all dependendency graphs
no_output_fields = 0
DO l = 1, no_dependency_graphs
  no_output_fields = no_output_fields                                          &
                     + SIZE(dependency_graph_list(l)%output_field)
END DO

! Allocate array to hold output field names
ALLOCATE(out_field_ids(no_output_fields))

! Fill output field name array with field names
n = 0
DO l = 1, no_dependency_graphs
  DO i = 1, SIZE(dependency_graph_list(l)%output_field)
    n = n + 1
    fid = dependency_graph_list(l)%output_field(i)%field_ptr%get_name()
    out_field_ids(n) = fid
  END DO
END DO
!
! End of interlude.

! 2) Check there are no output field overlap across dependency graphs.
l_output_not_unique = .false. ! Assume no overlap in output field lists

DO j = 1, no_output_fields ! Loop over all output fields

  ! Set output field name for comparison
  fid = out_field_ids(j)

  ! Check if the field id occurs multiple times in the output field name array.
  ! If it does, report issue to the user and set l_output_not_unique = .true.
  IF (COUNT(out_field_ids == fid) > 1) THEN

    l_output_not_unique = .true.
    WRITE(log_scratch_space, '(A)') 'Field ' // TRIM(fid) // ' occurs in ' //  &
                                    'output field list of multiple ' //        &
                                    'dependency graphs'
    CALL log_event(log_scratch_space, LOG_LEVEL_INFO)

  END IF

END DO ! End of loop over all output fields

! If there are any overlap between output field lists in the dependency graphs,
! report this to the user and abort the API.
IF (l_output_not_unique) THEN
  WRITE(log_scratch_space, '(A)')                                              &
          'Output fields not unique across dependency graphs'
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END IF

! 3) Check every dependency graph input field has a defined method of
!    generation, i.e its in the list of all dependency graph output fields.
l_input_field_not_generated = .false. ! Assume all input fields has a method of
                                      ! generation

DO l = 1, no_dependency_graphs ! Loop over all dependency graphs

  IF (ALLOCATED(dependency_graph_list(l)%input_field)) THEN ! Check if input
                                                            ! field list exists

    DO i = 1, SIZE(dependency_graph_list(l)%input_field) ! Loop over all input
                                                         ! fields

      ! Set name of input field for comparison
      fid = dependency_graph_list(l)%input_field(i)%field_ptr%get_name()

      ! Check whether input field name is amongst output field names across all
      ! dependency graphs. If not the report issue to user and set the logical
      ! flag l_input_field_not_generated = .true., i.e. the input field has no
      ! method of generation defined.
      IF (COUNT(out_field_ids == fid) == 0) THEN
        l_input_field_not_generated = .true.
        WRITE(log_scratch_space,'(A)') 'Field ' // TRIM(fid) // ' not in ' //  &
                                       'output list of any dependency graph.'
        CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
      END IF

    END DO ! End of loop over all input fields

  END IF ! End of check on input field existence

END DO  ! End of loop over dependency graphs

! If any dependency graph input field has not method of generation defined,
! report this to the user and abort the API.
IF (l_input_field_not_generated) THEN
  WRITE(log_scratch_space, '(A)')                                              &
        'Some dependency graph input fields have no defined generation method.'
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END IF
!
! Done with checks 1, 2 and 3
!

!
! Perform tasks 4 and 5
!
! Set-up dependency graph matrix. The dependency matrix encapsulates the
! dependencies between dependency graphs. If some of the input fields in
! dependency graph with index J, are in the output field list of dependency
! graph with index I, then this implies that dependency graph J depends on
! dependency graph I. Or in other words the output fields in dependency graph I
! has to be generated BEFORE the output fields in dependency graph J. In such a
! case the (i,j) entry in dependency matrix is set to 1, other wise it is 0.
!
DO j = 1, no_dependency_graphs ! Loop over all dependency graphs J

  DO i = 1, no_dependency_graphs ! Loop over all dependency graphs I

    dependency = .false. ! Assume dependency graph J do not dependency graph I

    DO l = 1, SIZE(dependency_graph_list(i)%output_field) ! Loop over all output
                                                          ! fields in dependency
                                                          !graph I

      ! Set field name of output field in dependency graph I for comparison
      fid = dependency_graph_list(i)%output_field(l)%field_ptr%get_name()

      dependency = .false. ! Explicit reset of flag for each output field to be
                           ! processed. Really unnecessary though.

      IF (ALLOCATED(dependency_graph_list(j)%input_field)) THEN ! Check if input
                                                                ! field list for
                                                                ! dependency
                                                                ! graph J exists


        DO n = 1, SIZE(dependency_graph_list(j)%input_field) ! Loop over input
                                                             ! fields for
                                                             ! dependency
                                                             ! graph J

          ! Evaluate logical that checks if output field in dependency graph I
          ! is same field as input field in dependency graph J
          dependency =                                                         &
          (dependency_graph_list(j)%input_field(n)%field_ptr%get_name() == fid)

        END DO ! End of loop over input fields for dependency graph J

      END IF ! End of check if input field list for dependency graph J exists

      IF (dependency) EXIT ! If ouput field in dependency graph I exists in
                           ! input field list of dependency graph J exit the
                           ! loop

    END DO ! End of loop over output fields in dependency graph I

    IF (dependency) THEN   ! dep graph j has a field dependency on dep graph i
      dep_graph_matrix(i,j) = 1
    ELSE                   ! dep graph j has no field dependency on dep graph i
      dep_graph_matrix(i,j) = 0
    END IF

  END DO ! Loop over all dependency graphs I

END DO ! Loop over all dependency graphs J

! Find/set order in which to create fields
CALL resolve_dependency(no_dependency_graphs, dep_graph_matrix,                &
                        generation_order)
!
! Done with tasks 4 and 5
!

END SUBROUTINE dependency_analyser


SUBROUTINE resolve_dependency(n, graph_matrix, order)

IMPLICIT NONE

!
! Arguments
!
! Number of nodes in full graph
INTEGER, INTENT(IN) :: n

! Dependency graph matrix
INTEGER, INTENT(IN) :: graph_matrix(n,n)

! Order in which nodes should be evaluated (the stack)
INTEGER, INTENT(OUT) :: order(n)

!
! Local variables
!
! Iterable
INTEGER :: i, j

! Logical flag for indicating whether dependencies can be resolved
LOGICAL :: cant_resolve_dependency

! Local copy of the dependency matrix
INTEGER :: graph_matrix_local(n,n)

! Make local copy of the dependency graph matrix
DO j = 1, n
  DO i = 1, n
    graph_matrix_local(i,j) = graph_matrix(i,j)
  END DO
END DO

! Initialise generation order stack to "empty"
DO i = 1, n
 order(i) = 0
END DO

! Start filling up the generation order stack.
i = 0
DO

  cant_resolve_dependency = .true.

  ! Search over all dependency graphs not yet on the generation order stack
  DO j = 1, n
    IF (ALL(order /= j)) THEN

      ! If node j has no remaining dependencies then add to next slot on
      ! the generation order stack and remove all dependencies on node j from
      ! dependency graph  matrix
      IF (ALL(graph_matrix_local(:,j) == 0)) THEN  ! No remaining dependencies
                                                   ! dependency graph j.
        cant_resolve_dependency = .false.

        i = i + 1                                  ! Pop dependency graph j on
        order(i) = j                               ! generation order stack.

        graph_matrix_local(j,:) = 0                ! Remove dependency graph j
                                                   ! from all dependencies in
                                                   ! dependency matrix.
      END IF

    END IF
  END DO

  ! If dependencies cannot be resolved stop further graph analysis
  IF (cant_resolve_dependency) EXIT

  ! Exit if all dependencies has been processed
  IF ((ALL(graph_matrix_local == 0)).AND.(ALL(order > 0))) EXIT

END DO

! Log an event
IF (cant_resolve_dependency) THEN
  CALL log_event('Cannot resolve graph dependencies', LOG_LEVEL_ERROR)
ELSE
  CALL log_event('Graph dependencies RESOLVED', LOG_LEVEL_INFO)
END IF

END SUBROUTINE resolve_dependency

END MODULE dependency_analyser_mod
