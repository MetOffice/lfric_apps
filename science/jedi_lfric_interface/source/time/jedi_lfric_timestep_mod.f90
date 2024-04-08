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

  use constants_mod,           only : i_timestep, r_second
  use jedi_lfric_duration_mod, only : jedi_duration_type
  use namelist_collection_mod, only : namelist_collection_type
  use namelist_mod,            only : namelist_type

  implicit none

  private
  public :: get_configuration_timestep

contains

  !> @brief  Get time step from the configuration object
  !>
  !> @param [in] configuration The configuration to extract timestep from
  function get_configuration_timestep( configuration ) result(timestep)

    implicit none

    type( namelist_collection_type ), intent(in) :: configuration

    type( jedi_duration_type ) :: timestep

    ! Local
    type( namelist_type ), pointer :: timestepping_nml
    real( kind=r_second )          :: dt

    ! Get configuration time-step
    timestepping_nml => configuration%get_namelist('timestepping')
    call timestepping_nml%get_value( 'dt', dt )
    call timestep%init( int(dt, i_timestep) )

  end function get_configuration_timestep

end module jedi_lfric_timestep_mod
