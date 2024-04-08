! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lfricinp_unit_handler_mod

USE, INTRINSIC :: iso_fortran_env, ONLY : int32, int64

IMPLICIT NONE

PRIVATE
PUBLIC :: get_free_unit

INTERFACE get_free_unit
MODULE PROCEDURE :: get_free_unit_32, get_free_unit_64
END INTERFACE get_free_unit

CONTAINS
!-------------------------------------------------------------------------------

SUBROUTINE get_free_unit_32(new_unit)
! Description:
!  Searches for a free unit number and returns first one it finds
IMPLICIT NONE
INTEGER(KIND=int32), INTENT(OUT) :: new_unit
LOGICAL :: unit_already_open
REAL :: rand

unit_already_open = .TRUE.

DO WHILE (unit_already_open)
  CALL RANDOM_NUMBER(rand)
  new_unit = 24000 + FLOOR(rand*1000)
  ! Check if unit already open, if true while loop will continue
  INQUIRE(UNIT=new_unit, OPENED=unit_already_open)
END DO

END SUBROUTINE get_free_unit_32

!-------------------------------------------------------------------------------

SUBROUTINE get_free_unit_64(new_unit)
IMPLICIT NONE
! Description:
!  Searches for a free unit number and returns first one it finds
INTEGER(KIND=int64), INTENT(OUT) :: new_unit
LOGICAL :: unit_already_open
REAL :: rand

unit_already_open = .TRUE.

DO WHILE (unit_already_open)
  CALL RANDOM_NUMBER(rand)
  new_unit = 24000 + FLOOR(rand*1000)
  ! Check if unit already open, if true while loop will continue
  INQUIRE(UNIT=new_unit, OPENED=unit_already_open)
END DO

END SUBROUTINE get_free_unit_64
!-------------------------------------------------------------------------------

END MODULE lfricinp_unit_handler_mod

