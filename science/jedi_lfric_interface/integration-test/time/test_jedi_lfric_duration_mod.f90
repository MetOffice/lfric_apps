!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!>@brief Tests for jedi_lfric_duration_mod
module test_jedi_lfric_duration_mod

  use jedi_lfric_duration_mod, only : jedi_duration_type

  implicit none

  private
  public test_duration_init_bad_string_err,     &
         test_duration_divide_zero_err,         &
         test_duration_divide_remainder_err,    &
         test_duration_divide_int_zero_err,     &
         test_duration_divide_int_remainder_err

contains

  !> @brief Test logging an error when trying to initialise
  !!        a jedi_duration with a bad string
  subroutine test_duration_init_bad_string_err( bad_string )

    implicit none

    character(*), intent(inout) :: bad_string

    type(jedi_duration_type)    :: jedi_duration

    call jedi_duration%init( bad_string )

  end subroutine test_duration_init_bad_string_err

  !> @brief Test logging an error when trying to divide
  !!        a jedi duration by a duration of zero
  subroutine test_duration_divide_zero_err()

    implicit none

    type(jedi_duration_type)    :: jedi_duration
    type(jedi_duration_type)    :: jedi_duration_2
    type(jedi_duration_type)    :: jedi_duration_3

    call jedi_duration%init() ! 0 seconds
    call jedi_duration_2%init( 4 )

    jedi_duration_3 = jedi_duration_2 / jedi_duration

  end subroutine test_duration_divide_zero_err

  !> @brief Test logging an error when two durations don't
  !!        divide evenly into each other
  subroutine test_duration_divide_remainder_err()

    implicit none

    type(jedi_duration_type)    :: jedi_duration
    type(jedi_duration_type)    :: jedi_duration_2
    type(jedi_duration_type)    :: jedi_duration_3

    call jedi_duration%init( 3 )
    call jedi_duration_2%init( 4 )

    jedi_duration_3 = jedi_duration_2 / jedi_duration

  end subroutine test_duration_divide_remainder_err

  !> @brief Test logging an error when trying to divide
  !!        a jedi duration by zero
  subroutine test_duration_divide_int_zero_err()

    implicit none

    type(jedi_duration_type)    :: jedi_duration
    type(jedi_duration_type)    :: jedi_duration_2

    call jedi_duration%init( 4 )

    jedi_duration_2 = jedi_duration / 0

  end subroutine test_duration_divide_int_zero_err

  !> @brief Test logging an error when a duration does
  !!        not divide evenly by the integer
  subroutine test_duration_divide_int_remainder_err()

    implicit none

    type(jedi_duration_type)    :: jedi_duration
    type(jedi_duration_type)    :: jedi_duration_2

    call jedi_duration%init( 4 )

    jedi_duration_2 = jedi_duration / 3

  end subroutine test_duration_divide_int_remainder_err

end module test_jedi_lfric_duration_mod
