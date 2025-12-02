!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Steps the shallow water miniapp through one timestep.
!!
!> @details Handles the stepping (for a single timestep) of the
!!          shallow water app, using either a semi-implicit or an
!!          explicit time stepping scheme.

module shallow_water_step_mod

  use constants_mod,                  only: i_def
  use field_mod,                      only: field_type
  use field_collection_mod,           only: field_collection_type
  use driver_modeldb_mod,             only: modeldb_type
  use log_mod,                        only: log_event,         &
                                            log_scratch_space, &
                                            LOG_LEVEL_INFO,    &
                                            LOG_LEVEL_ERROR
  use swe_timestep_alg_mod,           only: swe_timestep_alg_si,     &
                                            swe_timestep_alg_ssprk3, &
                                            swe_timestep_alg_rk4
  use shallow_water_settings_config_mod,                               &
                                      only: time_scheme,               &
                                            time_scheme_semi_implicit, &
                                            time_scheme_ssprk3,        &
                                            time_scheme_rk4

  implicit none

  private
  public shallow_water_step

  contains

  !> @brief Steps the shallow water miniapp through one timestep.
  !> @details Extracts prognostic fields from model_data then calls
  !!          swe_timestep_alg to step the shallow water miniapp
  !!          through one timestep. An iterated semi-implicit or an
  !!          explicit Runge-Kutta timestepping scheme is used.
  !> @param [in,out] modeldb      The working data set for the model run
  !>
  subroutine shallow_water_step( modeldb )

    implicit none

    type( modeldb_type ), target, intent(inout) :: modeldb

    type( field_collection_type ), pointer :: prognostic_fields => null()

    type(field_type), pointer :: wind => null()
    type(field_type), pointer :: geopot => null()
    type(field_type), pointer :: buoyancy => null()
    type(field_type), pointer :: q => null()
    type(field_type), pointer :: tracer_const => null()
    type(field_type), pointer :: tracer_pv => null()
    type(field_type), pointer :: tracer_step => null()
    type(field_type), pointer :: s_geopot => null()

    prognostic_fields => modeldb%fields%get_field_collection("prognostics")

    call prognostic_fields%get_field('wind', wind)
    call prognostic_fields%get_field('geopot', geopot)
    call prognostic_fields%get_field('buoyancy', buoyancy)
    call prognostic_fields%get_field('q', q)
    call prognostic_fields%get_field('s_geopot', s_geopot)
    call prognostic_fields%get_field('tracer_const', tracer_const)
    call prognostic_fields%get_field('tracer_pv', tracer_pv)
    call prognostic_fields%get_field('tracer_step', tracer_step)

    write( log_scratch_space, &
           '(A,I0)' ) 'Start of timestep ', modeldb%clock%get_step()
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

    select case( time_scheme )

    case( time_scheme_semi_implicit )
      call swe_timestep_alg_si( modeldb%clock,           &
                                wind,                    &
                                geopot, buoyancy, q,     &
                                tracer_const, tracer_pv, &
                                tracer_step, s_geopot )
    case ( time_scheme_ssprk3 )
      call swe_timestep_alg_ssprk3( modeldb%clock,       &
                                    wind,                &
                                    geopot, buoyancy, q, &
                                    s_geopot )
    case ( time_scheme_rk4 )
      call swe_timestep_alg_rk4( modeldb%clock,       &
                                 wind,                &
                                 geopot, buoyancy, q, &
                                 s_geopot )
    case default
      call log_event("No valid time stepping scheme selected", LOG_LEVEL_ERROR)

    end select

  end subroutine shallow_water_step

end module shallow_water_step_mod
