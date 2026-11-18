#if ! defined(DRHOOK)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module parkind1

!
! Description:
!   Dummy module to replace the DrHook library. Defines data types
!   which would otherwise be declared by the DrHook library.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: dummy_libraries
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.1.


implicit none

integer, parameter :: jpim = selected_int_kind(9)
integer, parameter :: jprb = selected_real_kind(13,300)

end module parkind1
#endif
