!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief   A module providing a method to get the time step duration from a
!>          configuration.
!>
module jedi_lfric_timestep_mod

  use config_mod,              only : config_type
  use constants_mod,           only : i_timestep, r_second
  use jedi_lfric_duration_mod, only : jedi_duration_type

  implicit none

  private
  public :: get_configuration_timestep

contains

  !> @brief  Get time step from the configuration object
  !>
  !> @param [in] config The configuration to extract timestep from
  function get_configuration_timestep( config ) result(timestep)

    implicit none

    type( config_type ), intent(in) :: config

    type( jedi_duration_type ) :: timestep

    ! Local
    real( kind=r_second ) :: dt

    ! Get configuration time-step
    dt = config%timestepping%dt()
    call timestep%init( int(dt, i_timestep) )

  end function get_configuration_timestep

end module jedi_lfric_timestep_mod
