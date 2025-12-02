!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
module test_cli_mod

  use, intrinsic :: iso_fortran_env, only : error_unit

  implicit none

  private
  public get_cli_argument

contains

  function get_cli_argument( idx ) result(argument)

    implicit none

    integer, intent(in), optional :: idx
    character(:), allocatable :: argument

    integer :: argument_index
    integer :: argument_length
    integer :: status

    if (present(idx)) then
      argument_index = idx
    else
      argument_index = 1
    end if

    call get_command_argument( argument_index,         &
                               length=argument_length, &
                               status=status )
    if (status /= 0) then
        write( error_unit, &
               '("Unable to get command line argument (", I0, ")")' ) status
        error stop 1
    end if

    allocate( character(argument_length) :: argument )
    call get_command_argument( argument_index, value=argument )

  end function get_cli_argument

end module test_cli_mod
