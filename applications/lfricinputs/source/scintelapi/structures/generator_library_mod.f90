! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE generator_library_mod
!
! This module provides access to the generator library. It also contains
! routines that creates/defines the library and finds the position index in the
! library of a generator with a given identifier.
!

USE dependency_graph_mod, ONLY: field_generator
USE log_mod,              ONLY: log_event, LOG_LEVEL_ERROR

IMPLICIT NONE

! Number of generators
INTEGER, PARAMETER :: no_generators = 14

! Generator array containing all generators
TYPE(field_generator), TARGET :: generator_list(no_generators)

CONTAINS

SUBROUTINE init_generator_lib()
!
! This routine creates/defines the generator library
!
! Basic tranforms
USE init_field_mod,           ONLY: init_field
USE read_from_file_mod,       ONLY: read_from_file
USE copy_field_data_mod,      ONLY: copy_field_data
USE a_times_X_dep_alg_mod,    ONLY: a_times_X_dep_alg
USE X_plus_Y_dep_alg_mod,     ONLY: X_plus_Y_dep_alg
USE aX_plus_bY_dep_alg_mod,   ONLY: aX_plus_bY_dep_alg
USE X_minus_Y_dep_alg_mod,    ONLY: X_minus_Y_dep_alg
USE X_times_Y_dep_alg_mod,    ONLY: X_times_Y_dep_alg
USE X_divideby_Y_dep_alg_mod, ONLY: X_divideby_Y_dep_alg
USE X_powint_n_mod,           ONLY: X_powint_n
USE X_powreal_a_mod,          ONLY: X_powreal_a
! Science transforms
USE soil_moist_content_to_soil_moist_stress_mod, ONLY:                         &
                                     soil_moist_content_to_soil_moist_stress
USE soil_moist_stress_to_soil_moist_content_mod, ONLY:                         &
                                     soil_moist_stress_to_soil_moist_content
USE specific_humidity_to_mixing_ratio_mod, ONLY: specific_humidity_to_mixing_ratio

IMPLICIT NONE

! Iterable
INTEGER :: l

!
! Basic transforms
!

! init field operator
l = 1
generator_list(l)%identifier = 'init_field                              '
generator_list(l)%generator => init_field

! read from dump
l = 2
generator_list(l)%identifier = 'read_from_file                          '
generator_list(l)%generator => read_from_file

! Copy data between fields
l = 3
generator_list(l)%identifier = 'copy_field_data                         '
generator_list(l)%generator => copy_field_data

! scale field operator
l = 4
generator_list(l)%identifier = 'a_times_X_dep_alg                       '
generator_list(l)%generator => a_times_X_dep_alg

! Add two fields
l = 5
generator_list(l)%identifier = 'X_plus_Y_dep_alg                        '
generator_list(l)%generator => X_plus_Y_dep_alg

! Linearly combine two fields
l = 6
generator_list(l)%identifier = 'aX_plus_bY_dep_alg                      '
generator_list(l)%generator => aX_plus_bY_dep_alg

! Subtract two fields
l = 7
generator_list(l)%identifier = 'X_minus_Y_dep_alg                       '
generator_list(l)%generator => X_minus_Y_dep_alg

! Multiply two fields
l = 8
generator_list(l)%identifier = 'X_times_Y_dep_alg                       '
generator_list(l)%generator => X_times_Y_dep_alg

! Divide two fields
l = 9
generator_list(l)%identifier = 'X_divideby_Y_dep_alg                    '
generator_list(l)%generator => X_divideby_Y_dep_alg

! Raise field to power n, where n is an integer)
l = 10
generator_list(l)%identifier = 'X_powint_n                              '
generator_list(l)%generator => X_powint_n

! Raise field to power n, where n is a real
l = 11
generator_list(l)%identifier = 'X_powreal_a                             '
generator_list(l)%generator => X_powreal_a

!
! Science transforms
!

! Soil moisture content to soil moisture stress conversion
l = 12
generator_list(l)%identifier = 'soil_moist_content_to_soil_moist_stress '
generator_list(l)%generator => soil_moist_content_to_soil_moist_stress

! Soil moisture content to soil moisture stress conversion
l = 13
generator_list(l)%identifier = 'soil_moist_stress_to_soil_moist_content '
generator_list(l)%generator => soil_moist_stress_to_soil_moist_content

! Convert specific humidities to equivalent mixing ratios
l = 14
generator_list(l)%identifier = 'spec_humidity_to_mixing_ratios'
generator_list(l)%generator => specific_humidity_to_mixing_ratio

END SUBROUTINE init_generator_lib


FUNCTION generator_index(generator_id)
!
! This function returns the library position index of a generator with a given
! id.
!

IMPLICIT NONE

!
! Arguments
!
! Generator identifier
CHARACTER(LEN=*), INTENT(IN) :: generator_id

!
! Local variables
!
! Generator index
INTEGER :: generator_index

! Iterable
INTEGER :: l

! Loop over generator library items to find generator index corresponding to id
generator_index = 0
DO l = 1, no_generators
  IF (trim(generator_id) == trim(generator_list(l)%identifier)) THEN
    generator_index = l
    EXIT
  END IF
END DO

! Check if generator has been found, if not raise error.
IF (generator_index == 0) THEN
  CALL log_event('Generator not found in generator library', LOG_LEVEL_ERROR)
END IF

END FUNCTION generator_index

END MODULE generator_library_mod
