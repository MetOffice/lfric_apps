!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module providing a configuration for the JEDI geometry emulator.
!>
!> @details A class is defined to hold the configuration data required to
!>          construct a JEDI geometry emulator. An initialiser is included and
!>          this is currently hard coded. In JEDI this information would be
!>          stored in a yaml file and eckit is used to parse and store a
!>          configuration object.
!>
module jedi_geometry_config_mod

  use constants_mod,               only : i_def, str_def
  use jedi_lfric_file_meta_mod,    only : jedi_lfric_file_meta_type
  use jedi_lfric_tests_config_mod, only : test_trajectory_path
  use jedi_lfric_duration_mod,     only : jedi_duration_type
  use jedi_lfric_datetime_mod,     only : jedi_datetime_type

  implicit none

  private

type, public :: jedi_geometry_config_type

  !> The file meta data
  type(jedi_lfric_file_meta_type) :: file_meta_data(2)
  !> The name of the context to setup
  character( len=str_def )        :: context_name
  !> The name of the namelist file that contains the configuration
  character(len=str_def)          :: filename
  !> The io clock time-step
  type( jedi_duration_type )      :: io_time_step
  !> The io calender start date
  type( jedi_datetime_type )      :: io_calender_start

contains

  !> Initialiser.
  procedure :: initialise

  !> jedi_geometry_config finalizer
  final     :: jedi_geometry_config_destructor

end type jedi_geometry_config_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains


!> @brief    Initialiser for jedi_geometry_config_type
!>
subroutine initialise( self, filename )

  implicit none

  class( jedi_geometry_config_type ), intent(inout) :: self
  character(len=*),                      intent(in) :: filename

  ! Local
  character( len=str_def ) :: xios_id
  character( len=str_def ) :: io_mode_str
  character( len=str_def ) :: field_group_id
  character( len=str_def ) :: file_name
  integer( kind=i_def )    :: freq

  ! The name of the namelist file that contains the configuration
  self%filename=filename

  ! Set the IO time-step
  call self%io_time_step%init( 'P0DT1H0M0S' )
  call self%io_calender_start%init( '2018-04-14T21:00:00' )

  ! The name of the context to use
  self%context_name="jedi_state"
  ! Setup two files - one read and one write
  ! Note, xios_id is arbitrary but has to be unique within a context
  ! Read
  file_name = test_trajectory_path
  xios_id = "read_model_data"
  io_mode_str = "read"
  field_group_id = "read_fields"
  freq=1_i_def
  call self%file_meta_data(1)%initialise( file_name,   &
                                          xios_id,     &
                                          io_mode_str, &
                                          freq,        &
                                          field_group_id )
  ! Write
  file_name = "write_model_data"
  xios_id = "write_model_data"
  io_mode_str = "write"
  field_group_id = "write_fields"
  freq=1_i_def
  call self%file_meta_data(2)%initialise( file_name,   &
                                          xios_id,     &
                                          io_mode_str, &
                                          freq,        &
                                          field_group_id )

end subroutine initialise

!> @brief    Finalizer for jedi_geometry_config_type
!>
subroutine jedi_geometry_config_destructor( self )

  implicit none

  type( jedi_geometry_config_type ), intent(inout)    :: self

end subroutine jedi_geometry_config_destructor

end module jedi_geometry_config_mod
