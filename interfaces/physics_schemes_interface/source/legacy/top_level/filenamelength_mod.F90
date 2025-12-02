! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Module for allocation of Filename length.

module filenamelength_mod

implicit none

! Description:
!   The module sets the max length of a filename

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: top_level

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.

integer, parameter      :: filenamelength = 256 ! Max length of a filename

!-------------------------------------------------------------------

end module filenamelength_mod
