!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> Module controlling the initialisation and finalisation of time related
!> functionality for JEDI IO. This module is a copy of "driver_time_mod.F90"
!> included in components/driver. The method "init_time" has been modified to
!> allow for a time-step to be specified. Currently the time-step is included
!> in a configuration module with global scope. Ongoing work to enable specific
!> configurations make it possible to specify the time-step via a configuration
!> object.This module is a stop gap until that functionality is available.
!>
module jedi_lfric_driver_time_mod

  use calendar_mod,      only: calendar_type
  use clock_mod,         only: clock_type
  use constants_mod,     only: i_def, r_second, i_timestep, str_def
  use log_mod,           only: log_event, log_level_error
  use model_clock_mod,   only: model_clock_type
  use step_calendar_mod, only: step_calendar_type

  implicit none

  private
  public :: jedi_lfric_init_time, jedi_lfric_final_time

contains

  !> Initialise clock and calendar from configuration
  !>
  !> @param[in]  time_step      The time step to use in the clock
  !> @param[in]  calendar_start The start date to use in the calendar
  !> @param[out] clock          The clock to create
  !> @param[out] calendar       The calendar to create
  subroutine jedi_lfric_init_time(time_step, calendar_start, clock, calendar)

    implicit none

    real(r_second),                      intent(in)  :: time_step
    character(str_def),                  intent(in)  :: calendar_start
    type(model_clock_type), allocatable, intent(out) :: clock
    class(calendar_type),   allocatable, intent(out) :: calendar

    ! Local
    integer(i_def)      :: rc
    integer(i_timestep) :: timestep_start
    integer(i_timestep) :: timestep_end
    real(r_second)      :: spinup_period
    character(str_def)  :: calendar_origin

    ! Create the clock - these values are hard coded for JEDI
    calendar_origin = calendar_start
    ! Choice of calendar here
    if (.not. allocated(calendar)) then
      allocate( calendar, source=step_calendar_type(calendar_origin, &
                                                    calendar_start), stat=rc )
      if (rc /= 0) then
        call log_event( "Unable to allocate calendar", log_level_error )
      end if
    end if

    ! Create the clock - these values are hard coded for JEDI
    timestep_start = 1_i_timestep
    timestep_end   = huge(1_i_timestep)
    spinup_period  = 0_r_second
    allocate( clock, source=model_clock_type( timestep_start, &
                                              timestep_end,   &
                                              time_step,      &
                                              max(spinup_period, 0.0_r_second) &
                                            ), stat=rc )

    if (rc /= 0) then
      call log_event( "Unable to allocate model clock", log_level_error )
    end if

  end subroutine jedi_lfric_init_time

  !> Finalise the clock and calendar
  !>
  !> @param[out] clock    The model clock
  !> @param[out] calendar The model calendar
  subroutine jedi_lfric_final_time(clock, calendar)

    implicit none

    type(model_clock_type),   allocatable, intent(out) :: clock
    class(calendar_type),     allocatable, intent(out) :: calendar

     if (allocated(clock))    deallocate(clock)
     if (allocated(calendar)) deallocate(calendar)

  end subroutine jedi_lfric_final_time

end module jedi_lfric_driver_time_mod
