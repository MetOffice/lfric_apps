! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module generator_library_mod
!
! This module provides access to the generator library. It also contains
! routines that creates/defines the library and finds the position index in the
! library of a generator with a given identifier.
!

use dependency_graph_mod, only: field_generator
use log_mod,              only: log_event, LOG_LEVEL_ERROR

implicit none

! Number of generators
integer, parameter :: no_generators = 14

! Generator array containing all generators
type(field_generator), target :: generator_list(no_generators)

contains

subroutine init_generator_lib()
!
! This routine creates/defines the generator library
!
! Basic tranforms
use init_field_mod,           only: init_field
use read_from_file_mod,       only: read_from_file
use copy_field_data_mod,      only: copy_field_data
use a_times_X_dep_alg_mod,    only: a_times_X_dep_alg
use X_plus_Y_dep_alg_mod,     only: X_plus_Y_dep_alg
use aX_plus_bY_dep_alg_mod,   only: aX_plus_bY_dep_alg
use X_minus_Y_dep_alg_mod,    only: X_minus_Y_dep_alg
use X_times_Y_dep_alg_mod,    only: X_times_Y_dep_alg
use X_divideby_Y_dep_alg_mod, only: X_divideby_Y_dep_alg
use X_powint_n_mod,           only: X_powint_n
use X_powreal_a_mod,          only: X_powreal_a
! Science transforms
use soil_moist_content_to_soil_moist_stress_mod, &
                              only: soil_moist_content_to_soil_moist_stress
use soil_moist_stress_to_soil_moist_content_mod, &
                              only: soil_moist_stress_to_soil_moist_content
use specific_humidity_to_mixing_ratio_mod, &
                              only: specific_humidity_to_mixing_ratio

implicit none

! Iterable
integer :: l

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

end subroutine init_generator_lib


function generator_index(generator_id)
!
! This function returns the library position index of a generator with a given
! id.
!

implicit none

!
! Arguments
!
! Generator identifier
character(len=*), intent(in) :: generator_id

!
! Local variables
!
! Generator index
integer :: generator_index

! Iterable
integer :: l

! Loop over generator library items to find generator index corresponding to id
generator_index = 0
do l = 1, no_generators
  if (trim(generator_id) == trim(generator_list(l)%identifier)) then
    generator_index = l
    exit
  end if
end do

! Check if generator has been found, if not raise error.
if (generator_index == 0) then
  call log_event('Generator not found in generator library', LOG_LEVEL_ERROR)
end if

end function generator_index

end module generator_library_mod
