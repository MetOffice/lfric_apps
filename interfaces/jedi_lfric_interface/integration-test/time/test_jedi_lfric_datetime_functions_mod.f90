!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!>@brief Tests the jedi_datetime_functions_mod routines that log an ERROR
module test_jedi_lfric_datetime_functions_mod

  use constants_mod,                 only : i_timestep

  implicit none

  private
  public test_YYYYMMDD_to_JDN,         &
         test_JDN_to_YYYYMMDD_invalid, &
         test_hhmmss_to_seconds,       &
         test_seconds_to_hhmmss_large, &
         test_seconds_to_hhmmss_neg

contains

  !> @brief Test logging an ERROR when an invalid date is passed
  subroutine test_YYYYMMDD_to_JDN()

    use jedi_lfric_datetime_functions_mod, only : YYYYMMDD_to_JDN

    implicit none

    integer(i_timestep) :: year  = -3000
    integer(i_timestep) :: month = 12
    integer(i_timestep) :: day   = 30

    integer(i_timestep) :: julian_day_number

    call YYYYMMDD_to_JDN( year, month, day, julian_day_number )

  end subroutine test_YYYYMMDD_to_JDN

  !> @brief Test logging an ERROR when an invalid JDN is passed
  subroutine test_JDN_to_YYYYMMDD_invalid()

    use jedi_lfric_datetime_functions_mod, only : JDN_to_YYYYMMDD

    implicit none

    integer(i_timestep) :: julian_day_number = -20000

    integer(i_timestep) :: year
    integer(i_timestep) :: month
    integer(i_timestep) :: day

    call JDN_to_YYYYMMDD( julian_day_number, year, month, day )

  end subroutine test_JDN_to_YYYYMMDD_invalid

  !> @brief Test logging an ERROR when an invalid time is passed
  subroutine test_hhmmss_to_seconds()

    use jedi_lfric_datetime_functions_mod, only : hhmmss_to_seconds

    implicit none

    integer(i_timestep) :: hour   = 34
    integer(i_timestep) :: minute = 12
    integer(i_timestep) :: second = 36

    integer(i_timestep) :: time_seconds

    call hhmmss_to_seconds( hour, minute, second, time_seconds )

  end subroutine test_hhmmss_to_seconds

  !> @brief Test logging an ERROR when the number of seconds
  !!        passed is too large
  subroutine test_seconds_to_hhmmss_large()

    use jedi_lfric_datetime_functions_mod, only : seconds_to_hhmmss

    implicit none

    integer(i_timestep) :: time_seconds = 86400

    integer(i_timestep) :: hour
    integer(i_timestep) :: minute
    integer(i_timestep) :: second

    call seconds_to_hhmmss( time_seconds, hour, minute, second )

  end subroutine test_seconds_to_hhmmss_large

  !> @brief Test logging an ERROR when the number of seconds
  !!        passed is negative
  subroutine test_seconds_to_hhmmss_neg()

    use jedi_lfric_datetime_functions_mod, only : seconds_to_hhmmss

    implicit none

    integer(i_timestep) :: time_seconds = -1

    integer(i_timestep) :: hour
    integer(i_timestep) :: minute
    integer(i_timestep) :: second

    call seconds_to_hhmmss( time_seconds, hour, minute, second )

  end subroutine test_seconds_to_hhmmss_neg

end module test_jedi_lfric_datetime_functions_mod
