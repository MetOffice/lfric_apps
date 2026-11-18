! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

! Data module containing missing data indicies

! Migrated from include file c_mdi.h

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: top_level
!
! Code Description:
!   Language: FORTRAN 90
!

module missing_data_mod

implicit none

     ! PP missing data indicator (-1.0E+30)
real, parameter    :: rmdi_pp  = -1.0e+30

! Old real missing data indicator (-32768.0)
real, parameter    :: rmdi_old = -32768.0

! New real missing data indicator (-2**30)
real, parameter    :: rmdi     = -32768.0*32768.0

! Integer missing data indicator
integer, parameter :: imdi     = -32768

end module missing_data_mod
