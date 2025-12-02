! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module dependency_graph_mod
!
! This module defines and provides access to the field generator and dependency
! graph derived type objects.
!

use field_ptr_mod,     only: field_ptr_object
use constants_def_mod, only: gen_id_len, genpar_len

implicit none

! Define the field generator derived type. It simply contains a identifier of
! intrinsic character type and a pointer to a dependency_graph_operator
! procedure
type field_generator
  character(len=gen_id_len) :: identifier
  procedure(dependency_graph_operator), pointer, nopass :: generator => null()
end type field_generator

! Define the dependency graph derived type. The input_field and output_field
! arrays will effectively contain pointers to the internal global field list of
! the API. The field generator gen is intended to be run on the dependency graph
! itself, taking the input field data and generate the output field data. The
! internal file (i.e. string) variable genpar is intended as a mechanism to pass
! parameters to the generator operating on the dependency graph. The procedure
! generate_fields should be called to run the generator gen on the the
! dependency graph (see definition of the routine run_generator below)
type dependency_graph
  type(field_ptr_object), allocatable :: input_field(:)
  type(field_ptr_object), allocatable :: output_field(:)
  type(field_generator) :: gen
  character(len=genpar_len) :: genpar
  contains
    procedure :: generate_fields => run_generator
end type dependency_graph

! Define the interface of a dependency graph operator procedure, which simply
! has a single dependency graph object as input and output. Note the import
! declaration is needed as the dependency graph derived type is defined in the
! same module.
abstract interface
  subroutine dependency_graph_operator(dep_graph)
    import :: dependency_graph
    class(dependency_graph), intent(inout) :: dep_graph
  end subroutine dependency_graph_operator
end interface

private :: run_generator

contains

subroutine run_generator(dep_graph)
!
! This routine takes a dependency graph object as input and run the generator,
! contained in dependency graph, on itself.
!

implicit none

!
! Argument definitions:
!
! Field object to run generator on
class(dependency_graph), intent(inout) :: dep_graph

! Call the generator on the dependency graph
call dep_graph % gen % generator(dep_graph)

end subroutine run_generator

end module dependency_graph_mod
