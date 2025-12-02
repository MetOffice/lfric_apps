! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module trelax_mod
implicit none
! Description: Module containing integer to select the
!   relaxation timescale for newtonian forcing setups
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: IDEALISED
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

integer, parameter :: tr_none=0
integer, parameter :: tr_HeldSuarez=1
integer, parameter :: tr_EL=2
integer, parameter :: tr_SHJ=3
integer, parameter :: tr_Jupiter=4
integer, parameter :: tr_HD209458b_Iro=5
integer, parameter :: tr_Y_Dwarf=6
integer, parameter :: tr_GJ1214b=7

end module trelax_mod
