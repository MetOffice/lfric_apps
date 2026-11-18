! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description: Derived model time/step information including start/end
!              step numbers and frequencies (in steps) of interface field
!              generation, boundary field updating, ancillary field
!              updating; and assimilation start/end times.
!              NB: Last three are set by IN_BOUND, INANCCTL, IN_ACCTL.
!              Also contains current time/date information, current
!              step number (echoed in history file) and steps-per-group.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: top_level
!
! Code Description:
!   Language: FORTRAN 90

module model_time_mod

use control_max_sizes, only: max_n_intf_a

implicit none

integer :: i_year               ! Current model time (years)
integer :: i_month              ! Current model time (months)
integer :: i_day                ! Current model time (days)
integer :: i_hour               ! Current model time (hours)
integer :: i_minute             ! Current model time (minutes)
integer :: i_second             ! Current model time (seconds)
integer :: i_day_number         ! Current model time (day no)

end module model_time_mod
