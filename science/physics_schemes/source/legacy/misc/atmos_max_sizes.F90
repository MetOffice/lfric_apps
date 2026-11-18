! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Misc
!
! Description
!   This module provides parameters giving the maximum likely sizes
!   of key UM resolution variables, useful for sizing static arrays.

module atmos_max_sizes

implicit none

integer, parameter :: model_levels_max = 250 ! max no of total levels

end module atmos_max_sizes
