! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE dependency_graph_list_mod
!
! This module defines and provides access to the global dependency graph list
! and a routine to initialise the list
!

USE dependency_graph_mod, ONLY: dependency_graph
USE constants_def_mod,    ONLY: max_no_dependency_graphs

IMPLICIT NONE

! 1d array of dependency graphs
TYPE(dependency_graph) :: dependency_graph_list(max_no_dependency_graphs)

! Array of indices showing order in which dependency graphs should be processed
INTEGER :: generation_order(max_no_dependency_graphs)

! Number of requested fields
INTEGER :: no_dependency_graphs

CONTAINS

SUBROUTINE init_dependency_graph_list()
!
! This routine initialises the dependency graph list
!

IMPLICIT NONE

! Iterable
INTEGER :: l

no_dependency_graphs = 0

DO l = 1, max_no_dependency_graphs
  generation_order(l) = 0
END DO

END SUBROUTINE init_dependency_graph_list

END MODULE dependency_graph_list_mod
