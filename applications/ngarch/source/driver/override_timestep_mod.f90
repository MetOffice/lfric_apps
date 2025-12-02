!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Module for overriding the GungHo timestep with an ngarch timestep
!> @details The structure of GungHo's timestep selection does not allow for
!>          new timesteps to be selected from without altering the source code.
!>          The override_timestep routine is triggered after GungHo's setup
!>          and replaces the GungHo timestep with an ngarch timestep, selected
!>          by the ngarch namelist in configuration.
!>          New timestep modules have to be included in the case statement if
!>          they are to be selectable, and the appropriate option needs to be
!>          set up in ngarch/rose-meta/lfric-ngarch/HEAD/rose-meta.conf.
module override_timestep_mod

  use timestep_method_mod, only: timestep_method_type
  use no_timestep_alg_mod, only: no_timestep_type
  use driver_modeldb_mod,  only: modeldb_type
  use ngarch_config_mod,   only: method, method_lfric_atm, method_casim, method_bl
  use log_mod,             only: log_event, LOG_LEVEL_ERROR

  use casim_timestep_mod,          only: casim_timestep_type
  use boundary_layer_timestep_mod, only: boundary_layer_timestep_type

  implicit none

contains

  !> @brief Replaces the GungHo timestep with the configured ngarch timestep.
  !> @param [in]  modeldb  The structure that holds model state
  subroutine override_timestep( modeldb )
    implicit none

    class( modeldb_type ), intent(inout) :: modeldb

    class( timestep_method_type ), pointer :: timestep_method => null()

    select case( method )
      case( method_lfric_atm )
        ! Do not overwrite the normal timestep

      case( method_casim )
        allocate( timestep_method, source=casim_timestep_type() )

      case( method_bl )
        allocate( timestep_method, source=boundary_layer_timestep_type() )

      case default
        call log_event( "ngarch: Invalid time stepping option chosen, "// &
                        "stopping program! ",LOG_LEVEL_ERROR )

    end select

    if ( associated( timestep_method ) ) then
      call modeldb%values%remove_key_value( "timestep_method" )
      call modeldb%values%add_key_value( "timestep_method", timestep_method )
    end if

  end subroutine override_timestep

end module override_timestep_mod
