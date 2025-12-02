!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the term_s
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Initialisation of iau fields. Set up the iau time control.
!>        Reads the iau_inc field from file.


module iau_time_control_mod

  use constants_mod,   only: r_def, i_def, l_def
  use model_clock_mod, only: model_clock_type
  use log_mod,         only: log_event,             &
                             log_scratch_space,     &
                             log_level_error
#ifdef UM_PHYSICS
  use iau_config_mod,  only: iau_ts_start,          &
                             iau_window_length,     &
                             iau_mode,              &
                             iau_mode_instantaneous
#endif

  implicit none

  public  :: calc_iau_ts_end, &
             calc_iau_ts_num, &
             calc_iau_weight

  contains

  !> @brief   Calculates the iau_ts_num variable.
  !> @details Returns the timestep index for a specified time
  !>          for use by the IAU.
  !> @param[in]  model_clock   The gungho clock
  !> @param[in]  iau_time      A time in seconds to convert into a timestep
  !>                           index. Can be specified as seconds from the
  !>                           beginning of the model run or as the IAU time
  !>                           window.
  !> @param[out] iau_ts_num    Timestep index for given time
  function calc_iau_ts_num( model_clock, iau_time ) result( timestep_index )

    implicit none

    class(model_clock_type), intent(in)    :: model_clock
    integer(kind=i_def), intent(in)        :: iau_time
    integer(kind=i_def)                    :: timestep_index
    real(kind=r_def)                       :: dt

    dt = 0.0_r_def
    timestep_index = 0.0_i_def
#ifdef UM_PHYSICS
    dt  = real(model_clock%get_seconds_per_step(), r_def)

    if ( dt <= 0.0_r_def ) then
      write(log_scratch_space,*) "dt cannot be leq 0, found dt = ", dt
      call log_event( log_scratch_space, log_level_error )
    end if

    timestep_index = nint( iau_time / dt ) + 1.0_i_def
#endif

  end function calc_iau_ts_num

  !> @brief   Calculates the iau_ts_end variable.
  !> @details Returns the end of the IAU application window.
  !> This is an integer valued index.
  !> @param[in]  model_clock   The gungho clock
  !> @param[out] iau_ts_end    End index of IAU time window
  function calc_iau_ts_end( model_clock ) result( iau_ts_end )

    implicit none

    class(model_clock_type), intent(in)    :: model_clock
#ifdef UM_PHYSICS
    integer(kind=i_def)                    :: iau_time
#endif
    integer(kind=i_def)                    :: iau_ts_end
    integer(kind=i_def)                    :: iau_ts_num

    iau_ts_num = 0.0_i_def
    iau_ts_end = 0.0_i_def
#ifdef UM_PHYSICS
    iau_time = iau_window_length
    iau_ts_num = calc_iau_ts_num( model_clock, iau_time )
    iau_ts_end = iau_ts_start + iau_ts_num - 1.0_r_def
#endif

  end function calc_iau_ts_end

  !> @brief   Calculates the iau_weight variable.
  !> @details Returns the weight of the IAU increments.
  !> Calculates the weight to apply to the total DA increment in order to add a
  !> uniform fraction at each timestep the iau is active for, such that the total
  !> increment is added over the iau time window. Supports weight calculation for
  !> standard increments and tendencies (quantity per second to add)
  !> @param[in]  model_clock   The gungho clock
  !> @param[in]  iau_tendency  Logical for increment being in the form of
  !>                           a tendency (quantity per second)
  !> @param[out] iau_weight    Weight of IAU increments
  function calc_iau_weight( model_clock, iau_tendency ) result( iau_weight )

    implicit none
    class(model_clock_type), intent(in)    :: model_clock
    logical(kind=l_def), intent(in)        :: iau_tendency
#ifdef UM_PHYSICS
    integer(kind=i_def)                    :: iau_time
#endif
    real(kind=r_def)                       :: iau_weight
    integer(kind=i_def)                    :: iau_ts_num

    iau_ts_num = 0.0_i_def
    iau_weight = 0.0_r_def
#ifdef UM_PHYSICS
    if (iau_mode == iau_mode_instantaneous) then
      iau_weight = 1.0_r_def
    else if ( iau_tendency ) then
      iau_weight = real( model_clock%get_seconds_per_step(), r_def )
    else
      iau_time = iau_window_length
      iau_ts_num = calc_iau_ts_num( model_clock, iau_time )
      iau_weight = 1.0_r_def / real( iau_ts_num, r_def )
    end if
#endif

  end function calc_iau_weight

end module iau_time_control_mod
