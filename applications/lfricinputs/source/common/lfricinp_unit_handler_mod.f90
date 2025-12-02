! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module lfricinp_unit_handler_mod

use, intrinsic :: iso_fortran_env, only : int32, int64

implicit none

private
public :: get_free_unit

interface get_free_unit
module procedure :: get_free_unit_32, get_free_unit_64
end interface get_free_unit

contains
!-------------------------------------------------------------------------------

subroutine get_free_unit_32(new_unit)
! Description:
!  Searches for a free unit number and returns first one it finds
implicit none
integer(kind=int32), intent(out) :: new_unit
logical :: unit_already_open
real :: rand

unit_already_open = .true.

do while (unit_already_open)
  call random_number(rand)
  new_unit = 24000 + floor(rand*1000)
  ! Check if unit already open, if true while loop will continue
  inquire(unit=new_unit, opened=unit_already_open)
end do

end subroutine get_free_unit_32

!-------------------------------------------------------------------------------

subroutine get_free_unit_64(new_unit)
implicit none
! Description:
!  Searches for a free unit number and returns first one it finds
integer(kind=int64), intent(out) :: new_unit
logical :: unit_already_open
real :: rand

unit_already_open = .true.

do while (unit_already_open)
  call random_number(rand)
  new_unit = 24000 + floor(rand*1000)
  ! Check if unit already open, if true while loop will continue
  inquire(unit=new_unit, opened=unit_already_open)
end do

end subroutine get_free_unit_64
!-------------------------------------------------------------------------------

end module lfricinp_unit_handler_mod

