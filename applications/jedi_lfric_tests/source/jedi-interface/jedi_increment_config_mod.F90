!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module providing a configuration for the JEDI increment emulator.
!>
!> @details A class is defined to hold the configuration data required to
!>          construct a JEDI increment emulator. An initialiser is included and
!>          this is currently hard coded. In JEDI this information would be
!>          stored in a yaml file and eckit is used to parse and store a
!>          configuration object.
!>
module jedi_increment_config_mod

  use constants_mod,             only : i_def, str_def, l_def
  use jedi_lfric_field_meta_mod, only : jedi_lfric_field_meta_type
  use fs_continuity_mod,         only : W3, Wtheta
  use jedi_lfric_datetime_mod,   only : jedi_datetime_type

  implicit none

  private

type, public :: jedi_increment_config_type

  !> The field meta data
  type( jedi_lfric_field_meta_type ) :: field_meta_data

  ! date: '2018-04-14T21:00:00Z'
  ! mpas define the formating for this as:
  ! dateTimeString = '$Y-$M-$D_$h:$m:$s'
  type( jedi_datetime_type )         :: inc_time
  !> File prefix for read
  character( len=str_def )           :: read_file_prefix

contains

  !> Initialiser.
  procedure :: initialise

  !> jedi_increment_config finalizer
  final     :: jedi_increment_config_destructor

end type jedi_increment_config_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains


!> @brief    Initialiser for jedi_increment_config_type
!>
subroutine initialise( self )

  implicit none

  class( jedi_increment_config_type ), intent(inout) :: self

  ! Local
  integer( kind=i_def ), parameter :: nvars = 10
  character( len=str_def )         :: variable_names(nvars)
  integer( kind=i_def )            :: variable_function_spaces(nvars)
  logical( kind=l_def )            :: variable_is_2d(nvars)

  ! Configuration inputs
  call self%inc_time%init( '2018-04-14T21:00:00' )
  self%read_file_prefix = "read_"

  ! Setup arrays required for field_meta_data
  ! Variable names
  variable_names(1) = "theta"
  variable_names(2) = "rho"
  variable_names(3) = "exner"
  variable_names(4) = "u_in_w3"
  variable_names(5) = "v_in_w3"
  variable_names(6) = "w_in_wth"
  variable_names(7) = "m_v"
  variable_names(8) = "m_cl"
  variable_names(9) = "m_r"
  variable_names(10) = "m_ci"

  ! Variable function spaces
  variable_function_spaces(1) = Wtheta
  variable_function_spaces(2) = W3
  variable_function_spaces(3) = W3
  variable_function_spaces(4) = W3
  variable_function_spaces(5) = W3
  variable_function_spaces(6) = Wtheta
  variable_function_spaces(7) = Wtheta
  variable_function_spaces(8) = Wtheta
  variable_function_spaces(9) = Wtheta
  variable_function_spaces(10) = Wtheta

  ! Variable is_2d
  variable_is_2d = .false.

  call self%field_meta_data%initialise( variable_names,           &
                                        variable_function_spaces, &
                                        variable_is_2d )

end subroutine initialise

!> @brief    Finalizer for jedi_increment_config_type
!>
subroutine jedi_increment_config_destructor( self )

  implicit none

  type( jedi_increment_config_type ), intent(inout)    :: self

end subroutine jedi_increment_config_destructor

end module jedi_increment_config_mod
