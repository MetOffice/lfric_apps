! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module dependency_graph_list_mod
!
! This module defines and provides access to the global dependency graph list
! and a routine to initialise the list
!

use dependency_graph_mod, only: dependency_graph
use constants_def_mod,    only: max_no_dependency_graphs

implicit none

! 1d array of dependency graphs
type(dependency_graph) :: dependency_graph_list(max_no_dependency_graphs)

! Array of indices showing order in which dependency graphs should be processed
integer :: generation_order(max_no_dependency_graphs)

! Number of requested fields
integer :: no_dependency_graphs

contains

subroutine init_dependency_graph_list()
!
! This routine initialises the dependency graph list
!

implicit none

! Iterable
integer :: l

no_dependency_graphs = 0

do l = 1, max_no_dependency_graphs
  generation_order(l) = 0
end do

end subroutine init_dependency_graph_list

end module dependency_graph_list_mod
