! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
PROGRAM scintelapi
!
! This is the top level driver program. All input to the API (i.e. namelist file
! names/paths) are parsed with the call to the routine scintelapi_initialise().
!

USE scintelapi_interface_mod, ONLY: scintelapi_add_dependency_graphs_from_nl,  &
                                    scintelapi_add_fields_from_nl,             &
                                    scintelapi_initialise, scintelapi_finalise
USE dependency_analyser_mod,  ONLY: dependency_analyser
USE dump_generator_mod,       ONLY: dump_generator

IMPLICIT NONE

! Initialise the LFRic, XIOS, API, etc. infrastructure
CALL scintelapi_initialise()

! Read namelist file for field definitions, and add said fields to internal
! field list
CALL scintelapi_add_fields_from_nl()

! Read namelist file for dependency graph definitions, and add said dependency
! graphs to internal list
CALL scintelapi_add_dependency_graphs_from_nl()

! Do dependency analysis
CALL dependency_analyser()

! Generate the dump
CALL dump_generator()

! Finalise the API
CALL scintelapi_finalise()

END PROGRAM scintelapi
