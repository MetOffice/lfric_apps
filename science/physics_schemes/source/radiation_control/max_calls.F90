! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Module to set the maximum number of calls to the radiation.
!
! Description:
!   This module declares the maximum number of calls to the
!   radiation code permitted on a single timestep.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: radiation_control
!
!- End of header

module max_calls

implicit none

!   Size allocated for arrays concerned with the SW call
integer, parameter :: npd_swcall=2

!   Size allocated for arrays concerned with the LW call
integer, parameter :: npd_lwcall=2

end module max_calls
