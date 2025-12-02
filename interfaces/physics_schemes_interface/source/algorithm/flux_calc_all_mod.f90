!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!>@brief A module to specify the surface latent and sensible heat fluxes
module flux_calc_all_mod
  use sci_time_interp_mod,          only: time_interpolate_list, &
                                          sinusoidal_function,   &
                                          diurnal_function
  use constants_mod,                only: i_def, r_def, rmdi
  use driver_modeldb_mod,           only: modeldb_type
  use specified_surface_config_mod, only: &
                                 time_data, profile_size,                &
                                 function_name_fluxes,                   &
                                 function_name_fluxes_sinusoidal,        &
                                 function_name_fluxes_diurnal,           &
                                 function_name_fluxes_constant,          &
                                 function_name_fluxes_time_interpolated, &
                                 function_amplitude_e,                   &
                                 function_period_e,                      &
                                 function_phase_e,                       &
                                 function_amplitude_h,                   &
                                 function_period_h,                      &
                                 function_phase_h,                       &
                                 specified_flux_e, specified_flux_h,     &
                                 time_units,                             &
                                 time_units_seconds, time_units_minutes, &
                                 time_units_hours, time_units_days,      &
                                 length_of_day, time_of_max_flux
   use log_mod,            only: log_event,                              &
                                 LOG_LEVEL_ERROR,                        &
                                 LOG_LEVEL_INFO

  implicit none

  public :: flux_calc_init, &
            flux_calc_step
  contains

!> @brief   Adds the fluxes to the modeldb
!> @param[in] modeldb The structure that holds the model state
subroutine flux_calc_init(modeldb)

  implicit none

  type(modeldb_type), intent(inout), target  :: modeldb

  real(r_def) :: flux_e
  real(r_def) :: flux_h

  flux_e = rmdi
  flux_h = rmdi

  call modeldb%values%add_key_value('flux_e', flux_e)
  call modeldb%values%add_key_value('flux_h', flux_h)

end subroutine flux_calc_init

!> @brief   Specify the surface sensible and latent heat fluxes
!> @details Specify the fluxes from constant values, time interpolated values
!!          or from analytical functions.
!> @param[in] modeldb The structure that holds the model state
subroutine flux_calc_step(modeldb)

  implicit none

  type(modeldb_type), intent(inout), target  :: modeldb

  real(r_def) :: time_now
  real(r_def), pointer :: flux_e
  real(r_def), pointer :: flux_h

  call modeldb%values%get_value('flux_e', flux_e)
  call modeldb%values%get_value('flux_h', flux_h)

  if ( function_name_fluxes /= function_name_fluxes_constant ) then

    ! Calculate current time in seconds since start date

    time_now = real( modeldb%clock%get_step(), r_def ) * &
               real( modeldb%clock%get_seconds_per_step(), r_def )

    ! Convert into appropriate units
    select case ( time_units )
      case( time_units_seconds )
        time_now = time_now
      case( time_units_minutes )
        time_now = time_now / 60
      case( time_units_hours )
        time_now = time_now / 3600
      case( time_units_days )
        time_now = time_now / 86400
      case default
        call log_event( 'Unknown time units', LOG_LEVEL_ERROR )
    end select

  end if

  select case ( function_name_fluxes )

    !-- Constant fluxes --
    case( function_name_fluxes_constant )
      flux_e = specified_flux_e(1)
      flux_h = specified_flux_h(1)

    !-- Time interpolated fluxes --
    ! This procedure uses the current time to interpolate
    ! a given list of values and times. If the current time falls
    ! outside of the range of the provided times list, then the max/min
    ! value is used instead accordingly.
    case( function_name_fluxes_time_interpolated )

      call time_interpolate_list( flux_e, time_data, &
                                  specified_flux_e, time_now )
      call time_interpolate_list( flux_h, time_data, &
                                  specified_flux_h, time_now )

    !-- Analytical function fluxes --
    ! Specify fluxes from a prescribed sine wave.
    case( function_name_fluxes_sinusoidal )

      call sinusoidal_function( flux_e, time_now,     &
                                function_amplitude_e, &
                                function_period_e,    &
                                function_phase_e )
      call sinusoidal_function( flux_h, time_now,     &
                                function_amplitude_h, &
                                function_period_h,    &
                                function_phase_h )

    !-- Analytical function fluxes --
    ! Specify fluxes from a prescribed diurnal cycle.
    case( function_name_fluxes_diurnal )
      call diurnal_function( flux_e, time_now,     &
                             function_amplitude_e, &
                             time_of_max_flux,     &
                             length_of_day )
      call diurnal_function( flux_h, time_now,     &
                             function_amplitude_h, &
                             time_of_max_flux,     &
                             length_of_day )

    case default
      call log_event( 'Unknown flux function name', LOG_LEVEL_ERROR )

  end select

end subroutine flux_calc_step

end module flux_calc_all_mod
