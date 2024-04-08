! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE dependency_graph_mod
!
! This module defines and provides access to the field generator and dependency
! graph derived type objects.
!

USE field_ptr_mod,     ONLY: field_ptr_object
USE constants_def_mod, ONLY: gen_id_len, genpar_len

IMPLICIT NONE

! Define the field generator derived type. It simply contains a identifier of
! intrinsic character type and a pointer to a dependency_graph_operator
! procedure
TYPE field_generator
  CHARACTER(len=gen_id_len) :: identifier
  PROCEDURE(dependency_graph_operator), POINTER, NOPASS :: generator => NULL()
END TYPE field_generator

! Define the dependency graph derived type. The input_field and output_field
! arrays will effectively contain pointers to the internal global field list of
! the API. The field generator gen is intended to be run on the dependency graph
! itself, taking the input field data and generate the output field data. The
! internal file (i.e. string) variable genpar is intended as a mechanism to pass
! parameters to the generator operating on the dependency graph. The procedure
! generate_fields should be called to run the generator gen on the the
! dependency graph (see definition of the routine run_generator below)
TYPE dependency_graph
  TYPE(field_ptr_object), ALLOCATABLE :: input_field(:)
  TYPE(field_ptr_object), ALLOCATABLE :: output_field(:)
  TYPE(field_generator) :: gen
  CHARACTER(len=genpar_len) :: genpar
  CONTAINS
    PROCEDURE :: generate_fields => run_generator
END TYPE dependency_graph

! Define the interface of a dependency graph operator procedure, which simply
! has a single dependency graph object as input and output. Note the IMPORT
! declaration is needed as the dependency graph derived type is defined in the
! same module.
ABSTRACT INTERFACE
  SUBROUTINE dependency_graph_operator(dep_graph)
    IMPORT :: dependency_graph
    CLASS(dependency_graph), intent(inout) :: dep_graph
  END SUBROUTINE dependency_graph_operator
END INTERFACE

PRIVATE :: run_generator

CONTAINS

SUBROUTINE run_generator(dep_graph)
!
! This routine takes a dependency graph object as input and run the generator,
! contained in dependency graph, on itself.
!

IMPLICIT NONE

!
! Argument definitions:
!
! Field object to run generator on
CLASS(dependency_graph), intent(inout) :: dep_graph

! Call the generator on the dependency graph
CALL dep_graph % gen % generator(dep_graph)

END SUBROUTINE run_generator

END MODULE dependency_graph_mod
