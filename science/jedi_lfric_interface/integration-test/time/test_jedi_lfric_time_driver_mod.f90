!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!>@brief Drives the execution of the jedi interface tests
module test_jedi_lfric_time_driver_mod

  use test_jedi_lfric_datetime_functions_mod, only : test_YYYYMMDD_to_JDN,         &
                                                     test_JDN_to_YYYYMMDD_invalid, &
                                                     test_hhmmss_to_seconds,       &
                                                     test_seconds_to_hhmmss_large, &
                                                     test_seconds_to_hhmmss_neg
  use test_jedi_lfric_datetime_mod,           only : test_init_string_err,             &
                                                     test_copy_from_jedi_datetime_err, &
                                                     test_add_duration_to_datetime,    &
                                                     test_duration_from_datetimes
  use test_jedi_lfric_duration_mod,           only : test_duration_init_bad_string_err,     &
                                                     test_duration_divide_zero_err,         &
                                                     test_duration_divide_remainder_err,    &
                                                     test_duration_divide_int_zero_err,     &
                                                     test_duration_divide_int_remainder_err

  implicit none

  private
  public test_jedi_interface_init,             &
         test_jedi_interface_final,            &
         run_init_string_err,                  &
         run_copy_from_jedi_datetime_err,      &
         run_add_duration_to_datetime,         &
         run_duration_from_datetimes,          &
         run_YYYYMMDD_to_JDN,                  &
         run_JDN_to_YYYYMMDD_invalid,          &
         run_hhmmss_to_seconds,                &
         run_seconds_to_hhmmss_large,          &
         run_seconds_to_hhmmss_neg,            &
         run_duration_init_bad_string_err,     &
         run_duration_divide_zero_err,         &
         run_duration_divide_remainder_err,    &
         run_duration_divide_int_zero_err,     &
         run_duration_divide_int_remainder_err

contains

  !> @brief Initialise testing for the jedi-interface.
  subroutine test_jedi_interface_init()

    implicit none

  end subroutine test_jedi_interface_init

  !> @brief Finalises testing for the jedi-interface
  subroutine test_jedi_interface_final()

    implicit none

  end subroutine test_jedi_interface_final

  !> @brief Tests logging an error when initialising a datetime
  !!        with a bad string
  subroutine run_init_string_err()

    implicit none

    call test_init_string_err()

  end subroutine run_init_string_err

  !> @brief Tests logging an error when attempting to initialise
  !!        a datetime with another uninitialised datetime
  subroutine run_copy_from_jedi_datetime_err()

    implicit none

    call test_copy_from_jedi_datetime_err()

  end subroutine run_copy_from_jedi_datetime_err

  !> @brief Tests adding a jedi_duration instance to a
  !!        jedi_datetime instance
  subroutine run_add_duration_to_datetime()

    implicit none

    call test_add_duration_to_datetime()

  end subroutine run_add_duration_to_datetime

  !> @brief Test getting a duration instance by subtracting
  !!        one jedi datetime from another
  subroutine run_duration_from_datetimes()

    implicit none

    call test_duration_from_datetimes()

  end subroutine run_duration_from_datetimes

  !> @brief Runs the YYYYMMDD_to_JDN test
  subroutine run_YYYYMMDD_to_JDN()

    implicit none

    call test_YYYYMMDD_to_JDN()

  end subroutine run_YYYYMMDD_to_JDN

  !> @brief Runs the JDN_to_YYYYMMDD test with an invalid JDN
  subroutine run_JDN_to_YYYYMMDD_invalid()

    implicit none

    call test_JDN_to_YYYYMMDD_invalid()

  end subroutine run_JDN_to_YYYYMMDD_invalid

  !> @brief Runs the hhmmss_to_seconds test
  subroutine run_hhmmss_to_seconds()

    implicit none

    call test_hhmmss_to_seconds()

  end subroutine run_hhmmss_to_seconds

  !> @brief Runs the hhmmss_to_seconds test with a too large time
  subroutine run_seconds_to_hhmmss_large()

    implicit none

    call test_seconds_to_hhmmss_large()

  end subroutine run_seconds_to_hhmmss_large

  !> @brief Runs the hhmmss_to_seconds test with a negative time
  subroutine run_seconds_to_hhmmss_neg()

    implicit none

    call test_seconds_to_hhmmss_neg()

  end subroutine run_seconds_to_hhmmss_neg

  !> @brief Runs the hhmmss_to_seconds test with a negative time
  subroutine run_duration_init_bad_string_err( bad_string )

    implicit none

    character(*), intent(inout) :: bad_string

    call test_duration_init_bad_string_err( bad_string )

  end subroutine run_duration_init_bad_string_err

  !> @brief Test logging an error when trying to divide
  !!        a jedi duration by zero
  subroutine run_duration_divide_zero_err()

    implicit none

    call test_duration_divide_zero_err()

  end subroutine run_duration_divide_zero_err

  !> @brief Test logging an error when two durations don't
  !!        divide evenly into each other
  subroutine run_duration_divide_remainder_err()

    implicit none

    call test_duration_divide_remainder_err()

  end subroutine run_duration_divide_remainder_err

  !> @brief Test logging an error when trying to divide
  !!        a jedi duration by an integer zero
  subroutine run_duration_divide_int_zero_err()

    implicit none

    call test_duration_divide_int_zero_err()

  end subroutine run_duration_divide_int_zero_err

  !> @brief Test logging an error when a durations does
  !!        not divide evenly by an integer
  subroutine run_duration_divide_int_remainder_err()

    implicit none

    call test_duration_divide_int_remainder_err()

  end subroutine run_duration_divide_int_remainder_err

end module test_jedi_lfric_time_driver_mod
