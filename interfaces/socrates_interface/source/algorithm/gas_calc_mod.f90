module gas_calc_mod

  use constants_mod, only: i_def, r_def, l_def
  use xios,          only: xios_date, xios_get_current_date
  use lfric_xios_utils_mod, only: seconds_from_date
  use radiative_gases_config_mod, only: l_cts_fcg_rates
  use log_mod,              only: log_event,         &
                                  log_scratch_space, &
                                  LOG_LEVEL_ERROR
  implicit none

contains

  subroutine gas_calc(gas_now, gas_index_max, gas_year, gas_conc, gas_rate)

    implicit none

    real(kind=r_def), intent(out)   :: gas_now ! gas concentration at time step
    integer(kind=i_def), intent(in) :: gas_index_max
    integer(kind=i_def), intent(in) :: gas_year(:)
    real(kind=r_def), intent(in)    :: gas_conc(:)
    real(kind=r_def), intent(in)    :: gas_rate(:)

    ! Date manipulation variables
    type(xios_date) :: datetime, start_date, end_date, start_date_p1y, time_origin
    real(r_def)     :: fraction

    ! other local variables
    integer(i_def) :: i, i_year, gas_index

    ! Get date and time
    call xios_get_current_date(datetime)
    i_year = int(datetime%year, i_def)

    ! Check that namelist data exists for this year
    if ( i_year < gas_year(1)) then
      call log_event( "gas_calc: no gas data for this year", LOG_LEVEL_ERROR )
    end if

    ! loop over years to find current index
    gas_index = 0
    DO i=1, gas_index_max
      IF ( i_year >= gas_year(i) ) gas_index = gas_index+1
    END DO

    start_date = xios_date(gas_year(gas_index), 1, 1, 0, 0, 0)
    end_date   = xios_date(gas_year(gas_index+1), 1, 1, 0, 0, 0)
    start_date_p1y = xios_date(gas_year(gas_index)+1, 1, 1, 0, 0, 0)

    if ( gas_rate(gas_index) > 0.0_r_def ) then

      ! Use input rates to calculate current concentration

      gas_now = gas_conc(1)
      ! Find last concentration, denoted by last year with negative gas_rate
      ! Then apply gas_rates prior to current year
      do i=1, gas_index-1
        if ( gas_rate(i) < 0.0_r_def ) then
          gas_now = gas_conc(i+1)
        else
          gas_now = gas_now * gas_rate(i) ** ( gas_year(i+1) - gas_year(i) )
        end if
      end do

      if (l_cts_fcg_rates) then
        fraction = fraction_of_period(start_date, start_date_p1y, datetime, continuous=.true.)
      else
        fraction = fraction_of_period(start_date, start_date_p1y, datetime, continuous=.false.)
      end if

      gas_now = gas_now * (gas_rate(gas_index) ** fraction)

    else
      ! Interpolate input concentrations to calculate current concentration
      gas_now = gas_conc(gas_index) +                               &
                ( gas_conc(gas_index+1) - gas_conc(gas_index) )     &
                * fraction_of_period(start_date, end_date, datetime, continuous=.true.)

    end if

  end subroutine gas_calc

  function fraction_of_period(date_start, date_end, current_date, continuous)
    ! There is a know limitation with this routine such that any date before the
    ! XIOS calendar origin will return -30931200s.  So this routine
    ! will only be correct for periods starting after or at the XIOS calendar origin.

    implicit none

    type(xios_date), intent(in) :: date_start
    type(xios_date), intent(in) :: date_end
    type(xios_date), intent(in) :: current_date
    logical(l_def),  intent(in) :: continuous

    real(r_def) :: fraction_of_period
    integer(i_def)  :: start_in_seconds
    integer(i_def)  :: period_in_seconds
    integer(i_def)  :: now_in_seconds

    start_in_seconds = int(seconds_from_date(date_start), i_def)
    period_in_seconds = int(seconds_from_date(date_end), i_def) - start_in_seconds
    now_in_seconds = int(seconds_from_date(current_date), i_def) - start_in_seconds
    if (continuous)then
      fraction_of_period = real( now_in_seconds, r_def) &
                              / real( period_in_seconds, r_def)
    else
      fraction_of_period = real( now_in_seconds / period_in_seconds, r_def)
    end if

  end function fraction_of_period

end module gas_calc_mod
