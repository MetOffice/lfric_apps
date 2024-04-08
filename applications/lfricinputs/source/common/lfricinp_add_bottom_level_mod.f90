! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lfricinp_add_bottom_level_mod

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY : real64, int64

IMPLICIT NONE

PRIVATE

PUBLIC :: lfricinp_add_bottom_level

CONTAINS

SUBROUTINE lfricinp_add_bottom_level(field)
! Add an extra bottom level to field

! Arguments
REAL(KIND=real64), ALLOCATABLE, INTENT(IN OUT) :: field(:,:)
! Note 1st index is 2D field, 2nd index is level number

! Local variables
REAL(KIND=real64), ALLOCATABLE :: field_temp(:,:)
INTEGER(KIND=int64) :: num_levels
INTEGER(KIND=int64) :: field_size_2d
INTEGER(KIND=int64) :: lev, i

num_levels = SIZE(field, 2)
field_size_2d = SIZE(field, 1)
! Allocate a temp field with an extra level
ALLOCATE(field_temp(field_size_2d, num_levels+1))
DO lev = 2, num_levels+1
  DO i = 1, field_size_2d
    field_temp(i, lev) = field(i, lev-1)
  END DO
END DO
! Take a copy of level
field_temp(:,1) = field_temp(:,2)
! Move across from temp field to proper field
CALL MOVE_ALLOC(field_temp, field)

END SUBROUTINE lfricinp_add_bottom_level

END MODULE lfricinp_add_bottom_level_mod
