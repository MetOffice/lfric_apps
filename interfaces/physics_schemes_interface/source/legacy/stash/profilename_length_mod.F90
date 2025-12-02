! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Module for allocation of stash profile name length.

module profilename_length_mod

implicit none

! Description:
!   The module sets the max length of a stash profile name

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: stash

! Code Description:
!   Language: FORTRAN 95

integer, parameter      :: profilename_length = 16 ! Max length of a profilename
integer, parameter      :: packagename_length = 40 ! Max length of a profilename
integer, parameter      :: fileid_length      = 20 ! Max length of a fileid
!-------------------------------------------------------------------

end module profilename_length_mod


