!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module providing a configuration for the pseudo model emulator.
!>
!> @details A class is defined to hold the configuration data required to
!>          construct the pseudo model emulator. An initialiser is included and
!>          this is currently hard coded. In JEDI this information would be
!>          stored in a yaml file and eckit is used to parse and store a
!>          configuration object.
!>
module jedi_pseudo_model_config_mod

  use constants_mod,           only : i_def, str_def
  use jedi_lfric_datetime_mod, only : jedi_datetime_type
  use jedi_lfric_duration_mod, only : jedi_duration_type

  implicit none

  private

type, public :: jedi_pseudo_model_config_type

  !> List of the dates to be read
  type( jedi_datetime_type ), allocatable :: state_times(:)

  !> File prefix for read
  character(len=str_def)                  :: read_file_prefix

contains

  !> jedi_pseudo_model initialiser.
  procedure :: initialise

  !> Finalizer
  final     :: jedi_pseudo_model_config_destructor

end type jedi_pseudo_model_config_type

!------------------------------------------------------------------------------
! Contained functions/subroutines
!------------------------------------------------------------------------------
contains

!> @brief    Initialiser for jedi_pseudo_model_config_type
!>
subroutine initialise( self )

  implicit none

  class( jedi_pseudo_model_config_type ), intent(inout) :: self

  ! Local
  integer( kind=i_def )      :: i, datetime_entries
  type( jedi_datetime_type ) :: next_datetime
  type( jedi_duration_type ) :: time_step

  call next_datetime%init( '2018-04-14T21:00:00' )
  call time_step%init( 'P0DT1H0M0S' )

  datetime_entries = 9_i_def
  allocate(self%state_times(datetime_entries))

  ! initialise datetime states 1-9 time steps after
  ! lfric calendar_start namelist variable time
  do i = 1, datetime_entries
    next_datetime = next_datetime + time_step
    self%state_times(i) = next_datetime
  end do

  self%read_file_prefix="read_"

end subroutine initialise

!> @brief    Finalizer for jedi_pseudo_model_config_type
!>
subroutine jedi_pseudo_model_config_destructor(self)!

  implicit none

  type(jedi_pseudo_model_config_type), intent(inout) :: self

  if ( allocated(self%state_times) ) deallocate(self%state_times)

end subroutine jedi_pseudo_model_config_destructor

end module jedi_pseudo_model_config_mod
