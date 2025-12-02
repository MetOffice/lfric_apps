! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top_level

module cderived_mod

implicit none

! Description:
!   This file contains declarations for derived constants and
!   grid definition coords in radians within the atmospheric model.

!
! Grid definition co-ordinates in radians
real:: Delta_lambda       ! EW (x) grid spacing in radians
real:: Delta_phi          ! NS (y) grid spacing in radians

end module cderived_mod
