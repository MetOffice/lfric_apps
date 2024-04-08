!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module providing a configuration for the JEDI linear model
!>        emulator.
!>
!> @details A class is defined to hold the configuration data required to
!>          construct a JEDI linear model emulator. An initialiser is included
!>          and this is currently hard coded. In JEDI this information would be
!>          stored in a yaml file and eckit is used to parse and store a
!>          configuration object.
!>
module jedi_linear_model_config_mod

  use constants_mod,             only : str_def
  use jedi_lfric_datetime_mod,   only : jedi_datetime_type
  use jedi_lfric_duration_mod,   only : jedi_duration_type

  implicit none

  private

type, public :: jedi_linear_model_config_type

  !> The TLM forecast_length
  type( jedi_duration_type ) :: forecast_length
  !> The configuration filename
  character( len=str_def )   :: config_filename

contains

  !> Initialiser.
  procedure :: initialise

  !> jedi_state_config finalizer
  final     :: jedi_state_config_destructor

end type jedi_linear_model_config_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

!> @brief    Initialiser for jedi_linear_model_config_type
!>
!> @param [in] config_filename The name of the configuration file
subroutine initialise( self, config_filename )

  implicit none

  class( jedi_linear_model_config_type ), intent(inout) :: self
  character(len=*),                       intent(in)    :: config_filename

  ! Configuration inputs

  ! Set initial model forecast length and configuration file name
  call self%forecast_length%init( 'P0DT6H0M0S' )
  self%config_filename = config_filename

end subroutine initialise

!> @brief    Finalizer for jedi_linear_model_config_type
!>
subroutine jedi_state_config_destructor( self )

  implicit none

  type( jedi_linear_model_config_type ), intent(inout)    :: self

end subroutine jedi_state_config_destructor

end module jedi_linear_model_config_mod
