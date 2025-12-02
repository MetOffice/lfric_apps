! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Module for allocation of error message length.

module errormessagelength_mod

use filenamelength_mod, only: filenamelength

implicit none

! Description:
!   The module sets the max length of several filename-related strings.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: top_level

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.

integer, parameter      :: errormessagelength = 2*filenamelength

!-------------------------------------------------------------------

end module errormessagelength_mod
