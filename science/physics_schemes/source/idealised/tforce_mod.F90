! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module tforce_mod

implicit none
!
! Description: Options for temperature profiles to be used
!              with idealised forcing.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: IDEALISED
!
! Code description:
!   Language: Fortran 90.
!   This code is written to programming standards in UMDP3.


integer, parameter :: tf_none=0
integer, parameter :: tf_HeldSuarez=1
integer, parameter :: tf_TLE=2
integer, parameter :: tf_EL=3
integer, parameter :: tf_SHJ=4
integer, parameter :: tf_Jupiter=5
integer, parameter :: tf_HD209458b_Heng=6
integer, parameter :: tf_HD209458b_Heng_smooth=7
integer, parameter :: tf_HD209458b_iro=8
integer, parameter :: tf_Y_Dwarf=9
integer, parameter :: tf_GJ1214b=10
integer, parameter :: tf_GJ1214b_dT800=11
integer, parameter :: tf_file=99
integer, parameter :: tf_isothermal=100

end module tforce_mod
