! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE read_from_file_mod
!
! This module contains a generator to read a field from file.
!

USE dependency_graph_mod, ONLY: dependency_graph

IMPLICIT NONE

CONTAINS

SUBROUTINE read_from_file(dep_graph)
!
! This generator reads the output field in the dependency graph from file
!

USE gen_io_check_mod,       ONLY: gen_io_check
use lfric_xios_read_mod,    ONLY: read_field_generic
USE field_mod,              ONLY: field_proxy_type
USE field_parent_mod,       ONLY: read_interface
USE log_mod,                ONLY: log_event, log_scratch_space, LOG_LEVEL_ERROR
USE constants_def_mod,      ONLY: genpar_len, field_name_len

IMPLICIT NONE

!
! Argument definitions:
!
! Dependency graph to be processed
CLASS(dependency_graph), INTENT(IN OUT) :: dep_graph

!
! Local variables
!
! IO procedure pointers
PROCEDURE(read_interface), POINTER :: tmp_read_ptr
!
! Parameter list
CHARACTER(LEN=genpar_len) :: parlist
!
! Id of field to read from file
CHARACTER(LEN=field_name_len) field_io_name

! Error code for reading parameter list
INTEGER :: ioerr

!
! Perform some initial input checks
!
CALL gen_io_check(                                                             &
                  dep_graph=dep_graph,                                         &
                  input_field_no=0,                                            &
                  output_field_no=1,                                           &
                  parameter_no=1                                               &
                 )
!
! Done with initial input checks
!

! Get id of field to read from file
parlist = dep_graph % genpar
READ(parlist,'(A)',IOSTAT=ioerr) field_io_name
IF (ioerr /= 0 ) THEN
  WRITE(log_scratch_space,'(A)') 'Error occured when parsing parameter ' //    &
                                 'list in routine READ_FROM_FILE. ' //         &
                                 'Input should be a single character string.'
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END IF

tmp_read_ptr => read_field_generic

! Read field from file
CALL dep_graph % output_field(1) % field_ptr %                                 &
                                   set_read_behaviour(tmp_read_ptr)
CALL dep_graph % output_field(1) % field_ptr %                                 &
                                   read_field('read_' // TRIM(field_io_name))

! Nullify read pointer
NULLIFY(tmp_read_ptr)

END SUBROUTINE read_from_file

END MODULE read_from_file_mod
