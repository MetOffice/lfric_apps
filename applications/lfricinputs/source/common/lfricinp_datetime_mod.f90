! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!> @brief     Module containing datetime TYPE
!> @details   Holds details of dates and times required by input and output
!!            files
!>

MODULE lfricinp_datetime_mod

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY : int64, real64

USE constants_mod, ONLY: i_def, r_second

IMPLICIT NONE
PRIVATE
PUBLIC :: datetime, datetime_type

TYPE :: datetime_type

  INTEGER(KIND=int64) :: num_times
  REAL(KIND=real64)   :: fctimes(99)
  REAL(KIND=real64)   :: first_fctime
  CHARACTER(LEN=19)   :: validity_times(99)
  CHARACTER(LEN=19)   :: first_validity_time
  CHARACTER(LEN=10)   :: calendar

  INTEGER(KIND=i_def) :: first_step
  INTEGER(KIND=i_def) :: last_step
  REAL(r_second)      :: spinup_period
  REAL(r_second)      :: seconds_per_step

CONTAINS

  PROCEDURE :: initialise

END TYPE datetime_type

TYPE(datetime_type) :: datetime

CONTAINS

SUBROUTINE initialise(self)

IMPLICIT NONE

CLASS(datetime_type) :: self

self % fctimes(:) = -1.0
self % validity_times(:) = 'XXXX-XX-XX XX:XX:XX'

self % num_times = 1
self % fctimes(1) = 0.0

self % first_fctime = 0.0
self % validity_times(1) = '2016-01-01 15:00:00'
self % first_validity_time = '2016-01-01 15:00:00'
self % calendar = 'Gregorian '

self % first_step = 1
self % last_step = 1
self % spinup_period = 0.0
self % seconds_per_step = 1.0

END SUBROUTINE initialise

END MODULE lfricinp_datetime_mod
