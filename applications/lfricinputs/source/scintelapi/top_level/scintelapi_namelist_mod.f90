! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE scintelapi_namelist_mod
!
! This module holds the paths to the LFRic infracstructure namelist file and the
! field definition and dependency graph namelist file. The former is used to
! initialise the LFRic infrastructure and the latter to configure the API.
!
! It also includes a routine to read the namelist paths from the command line
! argument list.
!

USE constants_def_mod, ONLY: file_name_len

IMPLICIT NONE

! Science Intelligence input namelist file
CHARACTER(LEN=file_name_len), PUBLIC :: scintelapi_nl

! Array containing required LFRic configuration namelists
CHARACTER(*), PARAMETER  :: required_lfric_namelists(6) = ['logging       ', &
                                                           'finite_element', &
                                                           'base_mesh     ', &
                                                           'planet        ', &
                                                           'extrusion     ', &
                                                           'io            ']


END MODULE scintelapi_namelist_mod
