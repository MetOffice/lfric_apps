! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module lfricinp_add_bottom_level_mod

! Intrinsic modules
use, intrinsic :: iso_fortran_env, only : real64, int64

implicit none

private

public :: lfricinp_add_bottom_level

contains

subroutine lfricinp_add_bottom_level(field)
! Add an extra bottom level to field

! Arguments
real(kind=real64), allocatable, intent(in out) :: field(:,:)
! Note 1st index is 2D field, 2nd index is level number

! Local variables
real(kind=real64), allocatable :: field_temp(:,:)
integer(kind=int64) :: num_levels
integer(kind=int64) :: field_size_2d
integer(kind=int64) :: lev, i

num_levels = size(field, 2)
field_size_2d = size(field, 1)
! Allocate a temp field with an extra level
allocate(field_temp(field_size_2d, num_levels+1))
do lev = 2, num_levels+1
  do i = 1, field_size_2d
    field_temp(i, lev) = field(i, lev-1)
  end do
end do
! Take a copy of level
field_temp(:,1) = field_temp(:,2)
! Move across from temp field to proper field
call move_alloc(field_temp, field)

end subroutine lfricinp_add_bottom_level

end module lfricinp_add_bottom_level_mod
