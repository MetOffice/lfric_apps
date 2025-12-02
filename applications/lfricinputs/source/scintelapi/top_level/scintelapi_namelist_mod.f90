! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module scintelapi_namelist_mod
!
! This module holds the paths to the LFRic infracstructure namelist file and the
! field definition and dependency graph namelist file. The former is used to
! initialise the LFRic infrastructure and the latter to configure the API.
!
! It also includes a routine to read the namelist paths from the command line
! argument list.
!

use constants_def_mod, only: file_name_len

implicit none

! Science Intelligence input namelist file
character(len=file_name_len), public :: scintelapi_nl

! Array containing required LFRic configuration namelists
character(*), parameter  :: required_lfric_namelists(6) = ['logging       ', &
                                                           'finite_element', &
                                                           'base_mesh     ', &
                                                           'planet        ', &
                                                           'extrusion     ', &
                                                           'io            ']


end module scintelapi_namelist_mod
