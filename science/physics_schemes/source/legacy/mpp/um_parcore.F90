! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP.

module UM_ParCore

! Note: These vars were formally part of parvars, but do not 'switch'
! with decomposition, and hence were inappropriately homed there, and also
! cause circular dependencies when used from very basic routines.

implicit none

! Initialise these for sane defaults in serial use cases
integer :: mype=0
integer :: nproc=1      ! number of processes
integer :: nproc_max=1  ! maximum number of processes

end module UM_ParCore
