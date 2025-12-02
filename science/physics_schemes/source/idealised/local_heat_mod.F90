! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module local_heat_mod

use missing_data_mod, only: rmdi, imdi

implicit none
!
! Description: Local heating options and parameters
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: IDEALISED
!
! Code description:
!   Language: Fortran 90.
!   This code is written to programming standards in UMDP3.

integer, parameter :: omit      = 0   ! no local heating
integer, parameter :: analytic  = 1   ! Heating centred on a point like
                                      ! the heating used by Oliver Halliday
                                      ! in an analytic study

integer :: local_heat_option = imdi

real    :: local_heat_xoffset = rmdi
real    :: local_heat_yoffset = rmdi
real    :: local_heat_amp     = rmdi
real    :: local_heat_sigma   = rmdi
real    :: local_heat_base    = rmdi
real    :: local_heat_top     = rmdi
real    :: local_heat_period  = rmdi


end module local_heat_mod
