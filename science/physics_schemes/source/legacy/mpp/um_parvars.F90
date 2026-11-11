! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP
!
! Provides variables that describe the current decomposition, and a routine to
! switch between decompositions


module UM_ParVars

implicit none

private

! maximum number of spatial dimensions
integer,parameter :: Ndim_max           = 3 ! 3d data

! position of personal data in global data (in terms of standard
! Fortran array notation
integer, public :: datastart(Ndim_max)

end module UM_ParVars
