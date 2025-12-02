! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
program scintelapi
!
! This is the top level driver program. All input to the API (i.e. namelist file
! names/paths) are parsed with the call to the routine scintelapi_initialise().
!

use scintelapi_interface_mod, only: scintelapi_add_dependency_graphs_from_nl,  &
                                    scintelapi_add_fields_from_nl,             &
                                    scintelapi_initialise, scintelapi_finalise
use dependency_analyser_mod,  only: dependency_analyser
use dump_generator_mod,       only: dump_generator

implicit none

! Initialise the LFRic, XIOS, API, etc. infrastructure
call scintelapi_initialise()

! Read namelist file for field definitions, and add said fields to internal
! field list
call scintelapi_add_fields_from_nl()

! Read namelist file for dependency graph definitions, and add said dependency
! graphs to internal list
call scintelapi_add_dependency_graphs_from_nl()

! Do dependency analysis
call dependency_analyser()

! Generate the dump
call dump_generator()

! Finalise the API
call scintelapi_finalise()

end program scintelapi
