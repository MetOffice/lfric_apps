! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!> @brief     Module containing datetime type
!> @details   Holds details of dates and times required by input and output
!!            files
!>

module lfricinp_datetime_mod

! Intrinsic modules
use, intrinsic :: iso_fortran_env, only : int64, real64

use constants_mod, only: i_def, r_second

implicit none
private
public :: datetime, datetime_type

type :: datetime_type

  integer(kind=int64) :: num_times
  real(kind=real64)   :: fctimes(99)
  real(kind=real64)   :: first_fctime
  character(len=19)   :: validity_times(99)
  character(len=19)   :: first_validity_time
  character(len=10)   :: calendar

  integer(kind=i_def) :: first_step
  integer(kind=i_def) :: last_step
  real(r_second)      :: spinup_period
  real(r_second)      :: seconds_per_step

contains

  procedure :: initialise

end type datetime_type

type(datetime_type) :: datetime

contains

subroutine initialise(self)

implicit none

class(datetime_type) :: self

self % fctimes(:) = -1.0_real64
self % validity_times(:) = 'XXXX-XX-XX XX:XX:XX'

self % num_times = 1
self % fctimes(1) = 0.0_real64

self % first_fctime = 0.0_real64
self % validity_times(1) = '2016-01-01 15:00:00'
self % first_validity_time = '2016-01-01 15:00:00'
self % calendar = 'Gregorian '

self % first_step = 1
self % last_step = 1
self % spinup_period = 0.0_r_second
self % seconds_per_step = 1.0_r_second

end subroutine initialise

end module lfricinp_datetime_mod
